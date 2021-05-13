#=
This script uses Ipopt (nonlinear solver) to optimize the turbine locations in a polygon-boundary wind farm.
=#

using Ipopt
using DelimitedFiles 
using PyPlot
import ForwardDiff


# ==================================================================================
# ============================= USER INPUTS IN THIS BLOCK ==========================
# ==================================================================================

# SET THE BOUNDARY FILE:
boundary_file = "input/boundary_files/example_boundary.txt"

# INITIAL LAYOUT FILE PATH:
initial_layout_file = "input/initial_layouts/example_initial_layout.txt"

# TURBINE PARAMETERS FILE PATH:
turbine_params_file = "input/turbine_files/vestas_V80_2MW/vestas_V80_2MW_param.jl"

# WIND ROSE FILE PATH:
windrose_file = "input/wind_resource/example_windrose_single_speed.yaml"

# MODEL SET FILE PATH:
model_set_file = "input/model_sets/model_set_1_example.jl"

# FINAL LAYOUT FILE PATH:
final_layout_file = "output/final_layouts/example_final_layout.txt"

# FINAL LAYOUT FIGURE FILE PATH:
final_layout_figure_file = "output/final_layout_figures/example_final_layout_figure.png"

# ==================================================================================


# ==================================================================================
# ========================= FLOWFARM WRAPPER FUNCTIONS =============================
# ==================================================================================

# set up boundary constraint wrapper function
function boundary_wrapper(x, params)
    # include relevant globals
    params.boundary_vertices
    params.boundary_normals

    # get number of turbines
    nturbines = Int(length(x)/2)
    
    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return boundary distances
    return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y)
end

# set up spacing constraint wrapper function
function spacing_wrapper(x, params)
    # include relevant globals
    params.rotor_diameter

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # get and return spacing distances
    return 2.5*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end

# set up objective wrapper function
function aep_wrapper(x, params)
    # include relevant globals
    params.turbine_z
    params.rotor_diameter
    params.hub_height
    params.turbine_yaw
    params.ct_models
    params.generator_efficiency
    params.cut_in_speed
    params.cut_out_speed
    params.rated_speed
    params.rated_power
    params.windresource
    params.power_models
    params.model_set
    params.rotor_points_y
    params.rotor_points_z
    params.obj_scale

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines] 
    turbine_y = x[nturbines+1:end]

    # calculate AEP
    AEP = obj_scale*ff.calculate_aep(turbine_x, turbine_y, turbine_z, rotor_diameter,
                hub_height, turbine_yaw, ct_models, generator_efficiency, cut_in_speed,
                cut_out_speed, rated_speed, rated_power, windresource, power_models, model_set,
                rotor_sample_points_y=rotor_points_y,rotor_sample_points_z=rotor_points_z)
    
    # return the objective as an array
    return [AEP]
end


# ==================================================================================
# ====================== FUNCTIONS FOR IPOPT (OPTIMIZER) ===========================
# ==================================================================================

# objective function
function obj(x) 
    return -aep_wrapper(x)[1]
end

# constraint function
function con(x, g)
    # calculate spacing constraint value
    spacing_con = spacing_wrapper(x)

    # calculate boundary constraint
    boundary_con = boundary_wrapper(x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[:] = [spacing_con; boundary_con]
end

# objective gradient function
function obj_grad(x, grad_f)
    grad_f[:] = -ForwardDiff.jacobian(aep_wrapper,x)
end

# constraint gradients function
function con_grad(x, mode, rows, cols, values)
    if mode == :Structure
        # report the sparcity structure of the jacobian
        for i = 1:prob.m
            rows[(i-1)*prob.n+1:i*prob.n] = i*ones(Int,prob.n)
            cols[(i-1)*prob.n+1:i*prob.n] = 1:prob.n
        end
    else
        # calculate spacing constraint jacobian
        ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

        # calculate boundary constraint jacobian
        db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

        # combine constaint jacobians into overall constaint jacobian arrays
        for i = 1:prob.m
            for j = 1:prob.n
                values[(i-1)*prob.n+j] = [ds_dx; db_dx][i, j]
            end
        end
    end
end


# ==================================================================================
# ========================= SET UP FLOWFARM PARAMETERS =============================
# ==================================================================================

# import model set with wind farm and related details
include(model_set_file)

# scale objective to be approximately between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
boundary_vertices = readdlm(boundary_file, skipstart=1) .* 1609.0   # convert from miles to meters
include("input/boundary_files/boundary_normals_calculator.jl")
boundary_normals = boundary_normals_calculator(boundary_vertices)

# set globals for use in wrapper functions
struct params_struct2{}
    model_set
    rotor_points_y
    rotor_points_z
    turbine_z
    ambient_ti
    rotor_diameter
    boundary_vertices
    boundary_normals
    obj_scale
    hub_height
    turbine_yaw
    ct_models
    generator_efficiency
    cut_in_speed
    cut_out_speed
    rated_speed
    rated_power
    windresource
    power_models
end

params = params_struct2(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_vertices, boundary_normals, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models)


# ==================================================================================
# ============================= SET UP OPTIMIZATION ================================
# ==================================================================================

# initialize design variable array
x = [copy(turbine_x);copy(turbine_y)]
global x

# report initial objective value
println("Starting AEP: ", aep_wrapper(x, params)[1]*1e-6/obj_scale, " MWh")

# add initial turbine location to plot
figure()
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false, color="C0", linestyle="--", alpha=0.3))
end

# get number of design variables
n_designvariables = length(x)

# get number of constraints
function numberofspacingconstraints(nturb)
    # calculates number of spacing constraints needed for given number of turbines
    ncon = 0
    for i = 1:nturb-1; ncon += i; end
    return ncon
end
n_spacingconstraints = numberofspacingconstraints(nturbines)
n_boundaryconstraints = length(boundary_wrapper(x, params))
n_constraints = n_spacingconstraints + n_boundaryconstraints

# set general lower and upper bounds for design variables
lb = ones(n_designvariables) * -Inf
ub = ones(n_designvariables) * Inf

# set lower and upper bounds for constraints
lb_g = ones(n_constraints) * -Inf
ub_g = zeros(n_constraints)

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)

# create the optimization problem
prob = createProblem(n_designvariables, lb, ub, n_constraints, lb_g, ub_g, n_designvariables*n_constraints, 0,
    obj, con, obj_grad, con_grad)
addOption(prob, "hessian_approximation", "limited-memory")

# set up for WEC (wake expansion coefficient) optimization
wec_steps = 6
wec_max = 3.0
wec_end = 1.0
wec_values = collect(LinRange(wec_max, wec_end, wec_steps))

# intialize 
xopt = fill(zeros(1), wec_steps)
fopt = fill(0.0, wec_steps)
info = fill("", wec_steps)


# ==================================================================================
# ========================== RUN AND TIME OPTIMIZATION =============================
# ==================================================================================

t1t = time()    # start optimization timer

# perform an optimization for each decreasing WEC value
for i = 1:length(wec_values)

    # set the WEC value in FlowFarm
    println("Running with WEC = ", wec_values[i])
    params.model_set.wake_deficit_model.wec_factor[1] = wec_values[i]

    # warm start the optimization with the most recent turbine coordinates found
    global x
    global xopt
    prob.x = deepcopy(x)

    # optimize
    t1 = time()
    status = solveProblem(prob)     # this runs the optimization
    t2 = time()
    clk = t2-t1

    # get results from the optimizer
    xopt[i] = prob.x           # design variables (turbine coordinates)
    fopt[i] = prob.obj_val     # objective value (annual energy production or AEP)
    info[i] = String(Ipopt.ApplicationReturnStatus[status])    # optimization info

    # print optimization results
    println("Finished in : ", clk, " (s)")
    println("Info: ", info)
    println("End objective value: ", -fopt[i])

    # save the optimized coordinates to be used as the starting point for the next optimization run
    x = deepcopy(xopt[i])

end

t2t = time()        # stop optimization timer
clkt = t2t - t1t    # calculate total optimization time


# ==================================================================================
# ========================== SHOW OPTIMIZATION RESULTS =============================
# ==================================================================================

# print optimization results
println("\n\n============ FINAL OPTIMIZATION RESULTS =============\n")
println("Finished in : ", clkt, " (s)")
println("Info: ", info)
println("Final AEP: ", aep_wrapper(xopt[end])[1]*1e-6/obj_scale, " MWh")
println()

# extract final turbine locations
turbine_x = copy(xopt[end][1:nturbines])
turbine_y = copy(xopt[end][nturbines+1:end])

# save final turbine locations to a text file
open(final_layout_file, "w") do io
    write(io, "# final turbine coordinates (x,y)\n")
    writedlm(io, [turbine_x turbine_y])
end

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false,color="C1")) 
end

# add wind farm boundary to plot
plt.gcf().gca().plot([boundary_vertices[:,1];boundary_vertices[1,1]],[boundary_vertices[:,2];boundary_vertices[1,2]], color="C2")

# set up and save/show plot
axis("square")
xlim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
ylim(minimum(boundary_vertices) - (maximum(boundary_vertices)-minimum(boundary_vertices))/5, maximum(boundary_vertices) + (maximum(boundary_vertices)-minimum(boundary_vertices))/5)
savefig(final_layout_figure_file, dpi=600)
# plt.show()