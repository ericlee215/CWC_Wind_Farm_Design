#=
This script uses Ipopt (nonlinear solver) to optimize the turbine locations in a polygon-boundary wind farm.
=#

# cd("Optimization")

using Ipopt
using DelimitedFiles 
using PyPlot
import ForwardDiff
include("../cableCost/cbl_analysis2.jl")
include("plot_boundary.jl")

# layout_number = lpad(ARGS[1],3,"0")
layout_number = "000"
max_cable_length = 50000.0

# ==================================================================================
# ============================= USER INPUTS IN THIS BLOCK ==========================
# ==================================================================================

# SET THE BOUNDARY FILE:
boundary_file = "input/boundary_files/boundary_winner_colome_3_utm.txt"

# INITIAL LAYOUT FILE PATH:
# initial_layout_file = "input/initial_layouts/initial-layout-$layout_number.txt"
initial_layout_file = "input/initial_layouts/initial-layout-302-16turb-2.txt"

# TURBINE PARAMETERS FILE PATH:
turbine_params_file = "input/turbine_files/Vestas_V110_2MW/Vestas_V110_2MW_param.jl"

# WIND ROSE FILE PATH:
windrose_file = "input/wind_resource/winner_colome_windrose_24dirs_avg_speed.yaml"

# MODEL SET FILE PATH:
model_set_file = "input/model_sets/model_set_winner_colome_avg_speed.jl"

# INTERMEDIATE LAYOUT FILE PATH:
intermediate_layout_file = "output/intermediate_layouts/winner_colome_$(layout_number)_intermediate.txt"

# INTERMEDIATE LAYOUT FIGURE FILE PATH:
intermediate_layout_figure_file = "output/intermediate_layout_figures/winner_colome_$(layout_number)_intermediate"

# FINAL LAYOUT FILE PATH:
final_layout_file = "output/final_layouts/winner_colome_$layout_number.txt"

# FINAL LAYOUT FIGURE FILE PATH:
final_layout_figure_file = "output/final_layout_figures/winner_colome_$layout_number"
final_layout_figure_file_color = "output/final_layout_figures/winner_colome_$(layout_number)_color"


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
    return ff.ray_trace_boundary(boundary_vertices, boundary_normals, turbine_x, turbine_y, smooth_max_factor=1)
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
    return 4.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end


function cable_length_wrapper(x, params)

    # include relevant globals
    max_cable_length = params.max_cable_length
    cable_nodes = params.cable_nodes
    # substation_coordinates = params.substation_coordinates

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # # append the substation as a node in the cable network
    # append!(turbine_x, substation_coordinates[1])
    # append!(turbine_y, substation_coordinates[2])

    # cable_length = zeros(typeof(turbine_x[1]), 1)

    # for i = 1:length(cable_nodes)
    #     cable_length[1] += sqrt((turbine_x[cable_nodes[i][1]] - turbine_x[cable_nodes[i][2]])^2 + (turbine_y[cable_nodes[i][1]] - turbine_y[cable_nodes[i][2]])^2)
    # end

    ############################
    cable_length = zeros(typeof(turbine_x[1]), 1)

    for i = 1:length(cable_nodes)
        cable_length[1] += sqrt((turbine_x[cable_nodes[i][1]] - turbine_x[cable_nodes[i][2]])^2 + (turbine_y[cable_nodes[i][1]] - turbine_y[cable_nodes[i][2]])^2)
    end
    ############################

    return [cable_length[1] - max_cable_length]
end

# function turbine_exclusion_wrapper(x, params)

#     # include relevant globals
#     turbine_circle_exclusions = params.turbine_circle_exclusions
#     turbine_polygon_exclusions = params.turbine_polygon_exclusions

#     # get number of turbines
#     nturbines = Int(length(x)/2)

#     # extract x and y locations of turbines from design variables vector
#     turbine_x = x[1:nturbines]
#     turbine_y = x[nturbines+1:end]

#     # initialize constraint vector
#     c = zeros(typeof(turbine_x[1]), (length(turbine_circle_exclusions) + length(turbine_polygon_exclusions)) * nturbines)

#     # circle exclusions
#     k = 1    
#     for i in turbine_circle_exclusions
#         # get and return boundary distances
#         c[k:k+nturbines-1] = ff.circle_boundary(i.center, i.radius, turbine_x, turbine_y)*1e-6
#         k += nturbines
#     end

#     # polygon exclusions
#     for i in turbine_polygon_exclusions
#         # get and return boundary distances
#         c[k:k+nturbines-1] = ff.ray_trace_boundary(i.vertices, i.normals, turbine_x, turbine_y)
#         k += nturbines
#     end

#     return c
# end


# function cable_exclusion_wrapper(x, params)

#     # include relevant globals
#     cable_circle_exclusions = params.cable_circle_exclusions
#     cable_polygon_exclusions = params.cable_polygon_exclusions


# end


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

    # calculate cable length constraint
    cable_length_con = cable_length_wrapper(x)

    # combine constaint values and jacobians into overall constaint value and jacobian arrays
    g[:] = [spacing_con; boundary_con; cable_length_con]

end

# function con_with_exclusions(x, g)
#     # calculate spacing constraint value
#     spacing_con = spacing_wrapper(x)

#     # calculate boundary constraint
#     boundary_con = boundary_wrapper(x)

#     # calculate cable length constraint
#     cable_length_con = cable_length_wrapper(x)

#     # # calculate exclusion constraints
#     exclusion_con = turbine_exclusion_wrapper(x)

#     # combine constaint values and jacobians into overall constaint value and jacobian arrays
#     g[:] = [spacing_con; boundary_con; cable_length_con; exclusion_con]

# end

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

        # calculate cable length sonctraint jacobian
        dcl_dx = ForwardDiff.jacobian(cable_length_wrapper, x)

        # combine constaint jacobians into overall constaint jacobian arrays
        for i = 1:prob.m
            for j = 1:prob.n
                values[(i-1)*prob.n+j] = [ds_dx; db_dx; dcl_dx][i, j]
            end
        end
    end
end

function con_grad2(x)

    # calculate spacing constraint jacobian
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint jacobian
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # calculate cable length sonctraint jacobian
    dcl_dx = ForwardDiff.jacobian(cable_length_wrapper, x)

    [ds_dx; db_dx; dcl_dx]

end

# # constraint gradients function
# function con_grad_with_exclusions(x, mode, rows, cols, values)
#     if mode == :Structure
#         # report the sparcity structure of the jacobian
#         for i = 1:prob.m
#             rows[(i-1)*prob.n+1:i*prob.n] = i*ones(Int,prob.n)
#             cols[(i-1)*prob.n+1:i*prob.n] = 1:prob.n
#         end
#     else
#         # calculate spacing constraint jacobian
#         ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

#         # calculate boundary constraint jacobian
#         db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

#         # calculate cable length sonctraint jacobian
#         dcl_dx = ForwardDiff.jacobian(cable_length_wrapper, x)

#         # calculate exclusion constraint jacobian
#         de_dx = ForwardDiff.jacobian(turbine_exclusion_wrapper, x)

#         # combine constaint jacobians into overall constaint jacobian arrays
#         for i = 1:prob.m
#             for j = 1:prob.n
#                 values[(i-1)*prob.n+j] = [ds_dx; db_dx; dcl_dx; de_dx][i, j]
#             end
#         end

#     end

# end

# ==================================================================================
# ========================= SET UP FLOWFARM PARAMETERS =============================
# ==================================================================================

# import model set with wind farm and related details
include(model_set_file)

# scale objective to be approximately between 0 and 1
obj_scale = 1E-11

# set wind farm boundary parameters
boundary_vertices = readdlm(boundary_file, skipstart=1) 
include("input/boundary_files/boundary_normals_calculator.jl")
boundary_normals = boundary_normals_calculator(boundary_vertices)
lb = [minimum(boundary_vertices[:,1]); minimum(boundary_vertices[:,2])]
ub = [maximum(boundary_vertices[:,1]); maximum(boundary_vertices[:,2])]
bounding_box_vertices = [lb[1] lb[2]; ub[1] lb[2]; ub[1] ub[2]; lb[1] ub[2]]
bounding_box_normals = boundary_normals_calculator(bounding_box_vertices)
boundary_vertices = deepcopy(bounding_box_vertices)
boundary_normals = deepcopy(bounding_box_normals)

# get substation location
substation_coordinates = readdlm("input/exclusions/substation/substation_utm.txt", skipstart=1)

# create structs to hold exclusion parameters
struct CircleExclusion{}
    center
    radius
end
# houses
houses_file_name = "input/exclusions/houses/houses_utm.txt"
houses = readdlm(houses_file_name, skipstart=1)
n_houses = length(houses[:,1])
# ponds
ponds_file_name = "input/exclusions/ponds/ponds_utm.txt"
ponds = readdlm(ponds_file_name, skipstart=1)
n_ponds = length(ponds[:,1])
# misc
misc_file_names = ["input/exclusions/misc/misc_utm.txt", "input/exclusions/misc/substation_utm.txt"]
miscs = [readdlm(misc_file_names[i], skipstart=1) for i=1:length(misc_file_names)]
n_miscs = length(miscs)
# get number of circular exclusions
n_circle_exclusions = n_houses + n_ponds + n_miscs
# initialize circular exclusions vector
turbine_circle_exclusions = Vector{CircleExclusion}(undef, n_circle_exclusions)
# populate the vector
k = 1
for i = 1:n_houses
    house_x = houses[i,1]
    house_y = houses[i,2]
    turbine_circle_exclusions[k] = CircleExclusion([house_x, house_y], 370.0)
    k += 1
end
for i = 1:n_ponds
    ponds_x = ponds[i,1]
    ponds_y = ponds[i,2]
    ponds_r = ponds[i,3]
    turbine_circle_exclusions[k] = CircleExclusion([ponds_x, ponds_y], ponds_r)
    k += 1
end
for misc in miscs
    misc_x = misc[1]
    misc_y = misc[2]
    misc_r = misc[3]
    turbine_circle_exclusions[k] = CircleExclusion([misc_x, misc_y], misc_r)
    k += 1
end

struct PolygonExclusion{}
    vertices
    normals
end
PolygonExclusion(a) = PolygonExclusion(a, boundary_normals_calculator(a))
# rivers
rivers_file_names = ["input/exclusions/rivers/river_1_utm.txt"]
rivers = [readdlm(rivers_file_names[i], skipstart=1) for i=1:length(rivers_file_names)]
n_rivers = length(rivers)
# roads
roads_file_names = ["input/exclusions/roads/road_$(lpad(i,2,"0"))_utm.txt" for i=1:11]
roads = [readdlm(roads_file_names[i], skipstart=1) for i=1:length(roads_file_names)]
n_roads = length(roads)
# get number of polygon exclusions
n_polygon_exclusions = n_rivers + n_roads
# initialize polygon exclusions vector
turbine_polygon_exclusions = Vector{PolygonExclusion}(undef, n_polygon_exclusions)
# populate the vector
k = 1
for river in rivers
    turbine_polygon_exclusions[k] = PolygonExclusion(river)
    k += 1
end
for road in roads
    turbine_polygon_exclusions[k] = PolygonExclusion(road)
    k += 1
end

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
    max_cable_length
    cable_nodes
    substation_coordinates
    turbine_circle_exclusions
    turbine_polygon_exclusions
end

params = params_struct2(model_set, rotor_points_y, rotor_points_z, turbine_z, ambient_ti, 
    rotor_diameter, boundary_vertices, boundary_normals, obj_scale, hub_height, turbine_yaw, 
    ct_models, generator_efficiency, cut_in_speed, cut_out_speed, rated_speed, rated_power, 
    windresource, power_models, max_cable_length, [[1,1] for i=1:nturbines-1], substation_coordinates,
    turbine_circle_exclusions, turbine_polygon_exclusions)



# create legend objects
legend_elements_turbines = [
    matplotlib.lines.Line2D([0], [0], marker="o", color="white", markerfacecolor="black", markersize=5, label="Turbine"),
    matplotlib.lines.Line2D([0], [0], marker="s", color="white", markerfacecolor="tab:green", label="Substation"),
    matplotlib.lines.Line2D([0], [0], color="tab:green", alpha=0.5, linewidth=0.5, linestyle="--", label="Cable")
    ]

# plot intermediate farm layout
legend_elements_boundary = plot_winner_colome_boundary(
    gray_only = true,
    plot_houses = false,
    plot_ponds = false,
    plot_rivers = false,
    plot_roads = false,
    plot_misc = false,
    save_fig = false
    )
        
# add legend(s)
plt.gcf().gca().add_artist(legend(handles=legend_elements_turbines, bbox_to_anchor=(1.05, 1), loc="upper left", ncol=1))

# plot cable layout
nodes_x = [turbine_x; params.substation_coordinates[1]]
nodes_y = [turbine_y; params.substation_coordinates[2]]
for i = 1:length(params.cable_nodes)
    plot(nodes_x[params.cable_nodes[i]], nodes_y[params.cable_nodes[i]], color="tab:green", alpha=0.5, linewidth=0.5, linestyle="--")
end

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=true, color="black")) 
end

# add substation to plot
substation_x_length = 200.0
substation_y_length = 150.0
plt.gcf().gca().add_artist(plt.Rectangle(
    (params.substation_coordinates[1] - substation_x_length/2, params.substation_coordinates[2] - substation_y_length/2),
    substation_x_length, 
    substation_y_length,
    fill=true, 
    color="tab:green"
    ))

savefig(intermediate_layout_figure_file * "_0.png", dpi=400)
savefig(intermediate_layout_figure_file * "_0.pdf")

# ==================================================================================
# ============================= SET UP OPTIMIZATION ================================
# ==================================================================================

x = [deepcopy(turbine_x);deepcopy(turbine_y)]

# report initial objective value
println("Starting AEP: ", aep_wrapper(x, params)[1]*1e-6/obj_scale, " MWh")

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
n_constraints = n_spacingconstraints + n_boundaryconstraints + 1
n_constraints_with_exclusions = n_constraints + (n_polygon_exclusions + n_circle_exclusions)*nturbines

# set general lower and upper bounds for design variables
# lb = [ones(nturbines) * minimum(boundary_vertices[:,1]); ones(nturbines) * minimum(boundary_vertices[:,2])]
# ub = [ones(nturbines) * maximum(boundary_vertices[:,1]); ones(nturbines) * maximum(boundary_vertices[:,2])]
lb = ones(n_designvariables) * -Inf
ub = ones(n_designvariables) * Inf

# set lower and upper bounds for constraints
lb_g = ones(n_constraints) * -Inf
ub_g = zeros(n_constraints)
lb_g_with_exclusions = ones(n_constraints_with_exclusions) * -Inf
ub_g_with_exclusions = zeros(n_constraints_with_exclusions)

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)
cable_length_wrapper(x) = cable_length_wrapper(x, params)
# turbine_exclusion_wrapper(x) = turbine_exclusion_wrapper(x, params)

# cost, _, cable_nodes = Cable_Analysis.cbl_analysis([x[1:nturbines]; substation_coordinates[1]], [x[nturbines+1:end]; substation_coordinates[2]], 1.0, return_network=true)
# params.cable_nodes[:] = cable_nodes
# println("\n\nCable length = ", cost)
# initial_con_values = con_grad2(x)

# create the optimization problem
prob = createProblem(n_designvariables, lb, ub, n_constraints, lb_g, ub_g, n_designvariables*n_constraints, 0,
    obj, con, obj_grad, con_grad)
addOption(prob, "hessian_approximation", "limited-memory")
addOption(prob, "max_iter", 10)
addOption(prob, "tol", 1e-2)
addOption(prob, "output_file", "output/opt_history_files/opt_history")


# set up for WEC (wake expansion coefficient) optimization
wec_steps = 6
wec_max = 3.0
wec_end = 1.0
n_wec_repeats = 1
wec_values_ref = collect(LinRange(wec_max, wec_end, wec_steps))
wec_values = Float64[]
for wec_value_ref = wec_values_ref
    append!(wec_values, ones(n_wec_repeats)*wec_value_ref)
end
n_wec = length(wec_values)
println(wec_values)

# intialize 
xopt = fill(zeros(1), n_wec)
fopt = fill(0.0, n_wec)
info = fill("", n_wec)

# ==================================================================================
# ========================== RUN AND TIME OPTIMIZATION =============================
# ==================================================================================

# perform an optimization for each decreasing WEC value
for i = 1:length(wec_values)

    # get the cable network
    # cost, _, cable_nodes = Cable_Analysis.cbl_analysis([x[1:nturbines]; substation_coordinates[1]], [x[nturbines+1:end]; substation_coordinates[2]], 1.0, return_network=true)
    cost, _, cable_nodes = Cable_Analysis.cbl_analysis(x[1:nturbines], x[nturbines+1:end], 1.0, return_network=true)
    params.cable_nodes[:] = cable_nodes
    println("\n\nCable length = ", cost)

    # set the WEC value in FlowFarm
    println("Running with WEC = ", wec_values[i])
    params.model_set.wake_deficit_model.wec_factor[1] = wec_values[i]

    # warm start the optimization with the most recent turbine coordinates found
    global x
    global xopt

    # optimize
    t1 = time()
    status = solveProblem(prob)     # this runs the optimization
    addOption(prob, "tol", 1e-1)
    t2 = time()
    clk = t2-t1

    # get results from the optimizer
    xopt[i] = prob.x           # design variables (turbine coordinates)
    fopt[i] = prob.obj_val     # objective value (annual energy production or AEP)
    info[i] = String(Ipopt.ApplicationReturnStatus[status])    # optimization info

    # print optimization results
    println("\n\nFinished in : ", clk, " (s)")
    println("Info: ", info)
    println("End objective value: ", -fopt[i])

    # save the optimized coordinates to be used as the starting point for the next optimization run
    x = deepcopy(xopt[i])

    # save final turbine locations to a text file
    open(intermediate_layout_file, "w") do io
        write(io, "# final turbine coordinates (x,y)\n")
        writedlm(io, [turbine_x turbine_y])
    end

    # plot intermediate farm layout
    legend_elements_boundary = plot_winner_colome_boundary(
        gray_only = true,
        plot_houses = false,
        plot_ponds = false,
        plot_rivers = false,
        plot_roads = false,
        plot_misc = false,
        save_fig = false
        )

    # get turbine positions
    turbine_x = copy(x[1:nturbines])
    turbine_y = copy(x[nturbines+1:end])
            
    # add legend(s)
    plt.gcf().gca().add_artist(legend(handles=legend_elements_turbines, bbox_to_anchor=(1.05, 1), loc="upper left", ncol=1))

    # plot cable layout
    nodes_x = [turbine_x; params.substation_coordinates[1]]
    nodes_y = [turbine_y; params.substation_coordinates[2]]
    for i = 1:length(params.cable_nodes)
        plot(nodes_x[params.cable_nodes[i]], nodes_y[params.cable_nodes[i]], color="tab:green", alpha=0.5, linewidth=0.5, linestyle="--")
    end

    # add final turbine locations to plot
    for i = 1:length(turbine_x)
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=true, color="black")) 
    end

    # add substation to plot
    substation_x_length = 200.0
    substation_y_length = 150.0
    plt.gcf().gca().add_artist(plt.Rectangle(
        (params.substation_coordinates[1] - substation_x_length/2, params.substation_coordinates[2] - substation_y_length/2),
        substation_x_length, 
        substation_y_length,
        fill=true, 
        color="tab:green"
        ))

    savefig(intermediate_layout_figure_file * "_$i.png", dpi=400)
    savefig(intermediate_layout_figure_file * "_$i.pdf")

end


# create new optimization problem with exclusions
prob = createProblem(n_designvariables, lb, ub, n_constraints_with_exclusions, lb_g_with_exclusions, ub_g_with_exclusions, n_designvariables*n_constraints_with_exclusions, 0,
    obj, con_with_exclusions, obj_grad, con_grad_with_exclusions)
addOption(prob, "hessian_approximation", "limited-memory")

# perform an optimization for each decreasing WEC value
for i = 1:length(wec_values)

    # get the cable network
    cost, _, cable_nodes = Cable_Analysis.cbl_analysis([x[1:nturbines]; substation_coordinates[1]], [x[nturbines+1:end]; substation_coordinates[2]], 1.0, return_network=true)
    params.cable_nodes[:] = cable_nodes
    println("\n\nCable length = ", cost)

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
    println("\n\nFinished in : ", clk, " (s)")
    println("Info: ", info)
    println("End objective value: ", -fopt[i])

    # save the optimized coordinates to be used as the starting point for the next optimization run
    x = deepcopy(xopt[i])

    # save final turbine locations to a text file
    open(intermediate_layout_file, "w") do io
        write(io, "# final turbine coordinates (x,y)\n")
        writedlm(io, [turbine_x turbine_y])
    end

    # plot intermediate farm layout
    legend_elements_boundary = plot_winner_colome_boundary(
        gray_only = true,
        plot_houses = true,
        plot_ponds = true,
        plot_rivers = true,
        plot_roads = true,
        plot_misc = true,
        save_fig = false
        )

    # get turbine positions
    turbine_x = copy(x[1:nturbines])
    turbine_y = copy(x[nturbines+1:end])
            
    # add legend(s)
    plt.gcf().gca().add_artist(legend(handles=legend_elements_turbines, bbox_to_anchor=(1.05, 1), loc="upper left", ncol=1))

    # plot cable layout
    nodes_x = [turbine_x; params.substation_coordinates[1]]
    nodes_y = [turbine_y; params.substation_coordinates[2]]
    for i = 1:length(params.cable_nodes)
        plot(nodes_x[params.cable_nodes[i]], nodes_y[params.cable_nodes[i]], color="tab:green", alpha=0.5, linewidth=0.5, linestyle="--")
    end

    # add final turbine locations to plot
    for i = 1:length(turbine_x)
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=true, color="black")) 
    end

    # add substation to plot
    substation_x_length = 200.0
    substation_y_length = 150.0
    plt.gcf().gca().add_artist(plt.Rectangle(
        (params.substation_coordinates[1] - substation_x_length/2, params.substation_coordinates[2] - substation_y_length/2),
        substation_x_length, 
        substation_y_length,
        fill=true, 
        color="tab:green"
        ))

    savefig(intermediate_layout_figure_file * "_with_exclusions_$i.png", dpi=400)
    savefig(intermediate_layout_figure_file * "_with_exclusions_$i.pdf")

end


# ==================================================================================
# ========================== SHOW OPTIMIZATION RESULTS =============================
# ==================================================================================

# print optimization results
println("\n\n============ FINAL OPTIMIZATION RESULTS =============\n")
# println("Finished in : ", clkt, " (s)")
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

# plot gray farm boundary
legend_elements_boundary = plot_winner_colome_boundary(
    gray_only = true,
    plot_houses = true,
    plot_ponds = true,
    plot_rivers = true,
    plot_roads = true,
    plot_misc = true,
    save_fig = false
    )

# add legend(s)
plt.gcf().gca().add_artist(legend(handles=legend_elements_turbines, bbox_to_anchor=(1.05, 1), loc="upper left", ncol=1))

# plot cable layout
nodes_x = [turbine_x; params.substation_coordinates[1]]
nodes_y = [turbine_y; params.substation_coordinates[2]]
for i = 1:length(params.cable_nodes)
    plot(nodes_x[params.cable_nodes[i]], nodes_y[params.cable_nodes[i]], color="tab:green", alpha=0.5, linewidth=0.5, linestyle="--")
end

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=true, color="black")) 
end

# add substation to plot
substation_x_length = 200.0
substation_y_length = 150.0
plt.gcf().gca().add_artist(plt.Rectangle(
    (params.substation_coordinates[1] - substation_x_length/2, params.substation_coordinates[2] - substation_y_length/2),
    substation_x_length, 
    substation_y_length,
    fill=true, 
    color="tab:green"
    ))

savefig(final_layout_figure_file * ".png", dpi=600)
savefig(final_layout_figure_file * ".pdf")



# plot color farm boundary
legend_elements_boundary = plot_winner_colome_boundary(
    gray_only = false,
    plot_houses = true,
    plot_ponds = true,
    plot_rivers = true,
    plot_roads = true,
    plot_misc = true,
    save_fig = false
    )

# add legend(s)
plt.gcf().gca().add_artist(legend(handles=legend_elements_turbines, bbox_to_anchor=(1.05, 1), loc="upper left", ncol=1))
plt.gcf().gca().add_artist(legend(handles=legend_elements_boundary, title="Restricted area types", bbox_to_anchor=(1.05, 0.5), loc="upper left", ncol=1))

# plot cable layout
nodes_x = [turbine_x; params.substation_coordinates[1]]
nodes_y = [turbine_y; params.substation_coordinates[2]]
for i = 1:length(params.cable_nodes)
    plot(nodes_x[params.cable_nodes[i]], nodes_y[params.cable_nodes[i]], color="tab:green", alpha=0.5, linewidth=0.5, linestyle="--")
end

# add final turbine locations to plot
for i = 1:length(turbine_x)
    plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=true, color="black")) 
end

# add substation to plot
substation_x_length = 200.0
substation_y_length = 150.0
plt.gcf().gca().add_artist(plt.Rectangle(
    (params.substation_coordinates[1] - substation_x_length/2, params.substation_coordinates[2] - substation_y_length/2),
    substation_x_length, 
    substation_y_length,
    fill=true, 
    color="tab:green"
    ))

savefig(final_layout_figure_file_color * ".png", dpi=600)
savefig(final_layout_figure_file_color * ".pdf")
