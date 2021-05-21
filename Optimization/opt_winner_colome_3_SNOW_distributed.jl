#=
This script uses Ipopt (nonlinear solver) to optimize the turbine locations in a polygon-boundary wind farm.
=#

# cd("Optimization")

using Snopt
using SNOW
using DelimitedFiles 
using PyPlot
import ForwardDiff
using Distributed
using ClusterManagers

include("../cableCost/cbl_analysis2.jl")
include("plot_farm_winner_colome.jl")
include("exclusions_winner_colome.jl")
include("opt_functions_winner_colome.jl")
addprocs(SlurmManager(parse(Int, ENV["SLURM_NTASKS"])-1))


# ==================================================================================
# ============================= USER INPUTS IN THIS BLOCK ==========================
# ==================================================================================

layout_number = lpad(ARGS[1],3,"0")
max_cable_length = parse(Float64, ARGS[2])
plot_intermediate_layouts = false


# SET THE BOUNDARY FILE:
boundary_file = "input/boundary_files/boundary_winner_colome_3_utm.txt"

# INITIAL LAYOUT FILE PATH:
initial_layout_file = "input/initial_layouts/initial-layout-$layout_number.txt"
@everywhere initial_layout_file = $initial_layout_file

# TURBINE PARAMETERS FILE PATH:
turbine_params_file = "input/turbine_files/Vestas_V110_2MW/Vestas_V110_2MW_param.jl"
@everywhere turbine_params_file = $turbine_params_file

# WIND ROSE FILE PATH:
windrose_file = "input/wind_resource/winner_colome_windrose_24dirs_avg_speed.yaml"
@everywhere windrose_file = $windrose_file

# MODEL SET FILE PATH:
model_set_file = "input/model_sets/model_set_winner_colome_avg_speed.jl"
@everywhere model_set_file = $model_set_file

# INTERMEDIATE LAYOUT FILE PATH:
intermediate_layout_file = "output/pareto/$(Int(max_cable_length))_meter_cable/intermediate_layouts/winner_colome_$(layout_number)_intermediate.txt"

# INTERMEDIATE LAYOUT FIGURE FILE PATH:
intermediate_layout_figure_file = "output/pareto/$(Int(max_cable_length))_meter_cable/intermediate_layout_figures/winner_colome_$(layout_number)_intermediate"

# FINAL LAYOUT FILE PATH:
final_layout_file = "output/pareto/$(Int(max_cable_length))_meter_cable/final_layouts/winner_colome_$layout_number.txt"

# FINAL LAYOUT FIGURE FILE PATH:
final_layout_figure_file = "output/pareto/$(Int(max_cable_length))_meter_cable/final_layout_figures/winner_colome_$layout_number"
final_layout_figure_file_color = "output/pareto/$(Int(max_cable_length))_meter_cable/final_layout_figures/winner_colome_$(layout_number)_color"

# ==================================================================================



# ==================================================================================
# ========================= SET UP FLOWFARM PARAMETERS =============================
# ==================================================================================

# import model set with wind farm and related details
@everywhere include(model_set_file)

# scale objective to be approximately between 0 and 1
obj_scale = 1E-5

# set wind farm boundary parameters
boundary_vertices = readdlm(boundary_file, skipstart=1) 
include("input/boundary_files/boundary_normals_calculator.jl")
boundary_normals = boundary_normals_calculator(boundary_vertices)

# get substation location
substation_coordinates = readdlm("input/exclusions/substation/substation_utm.txt", skipstart=1)

# import exclusions
turbine_circle_exclusions, turbine_polygon_exclusions = get_winner_colome_exclusions()
n_circle_exclusions = length(turbine_circle_exclusions)
n_polygon_exclusions = length(turbine_polygon_exclusions)

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
    windresource, power_models, max_cable_length, [[1,1] for i=1:nturbines], substation_coordinates,
    turbine_circle_exclusions, turbine_polygon_exclusions)


# ==================================================================================
# ============================= SET UP OPTIMIZATION ================================
# ==================================================================================

# get starting point
x = [deepcopy(turbine_x);deepcopy(turbine_y)]

# plot initial layout
if plot_intermediate_layouts
    plot_farm_winner_colome(intermediate_layout_figure_file, x, params,
        figure_number = 1,
        ext = [".png"],#, ".pdf"],
        gray_only = true,
        plot_houses = false,
        plot_ponds = false,
        plot_rivers = false,
        plot_roads = false,
        plot_misc = false,
        save_fig = true
        )
end

# report initial objective value
initial_aep = -aep_wrapper(x, params)[1]
println("Starting AEP: ", initial_aep*1e-6/obj_scale, " MWh")

# get number of design variables
n_designvariables = length(x)

# get number of constraints
n_spacingconstraints = Int((nturbines-1)*nturbines/2)
n_boundaryconstraints = length(boundary_wrapper(x, params))
n_constraints = n_spacingconstraints + n_boundaryconstraints + 1
n_constraints_with_exclusions = n_constraints + (n_circle_exclusions + n_polygon_exclusions)*nturbines

# set general lower and upper bounds for design variables
lb = [ones(nturbines) * minimum(boundary_vertices[:,1]); ones(nturbines) * minimum(boundary_vertices[:,2])]
ub = [ones(nturbines) * maximum(boundary_vertices[:,1]); ones(nturbines) * maximum(boundary_vertices[:,2])]

# set lower and upper bounds for constraints
lb_g = ones(n_constraints) * -Inf
ub_g = zeros(n_constraints)
lb_g_with_exclusions = ones(n_constraints_with_exclusions) * -Inf
ub_g_with_exclusions = zeros(n_constraints_with_exclusions)

# generate wrapper function surrogates
spacing_wrapper(x) = spacing_wrapper(x, params)
@everywhere aep_wrapper(x) = aep_wrapper(x, params)
boundary_wrapper(x) = boundary_wrapper(x, params)
cable_length_wrapper(x) = cable_length_wrapper(x, params)
turbine_exclusion_wrapper(x) = turbine_exclusion_wrapper(x, params)

# set up WEC (wake expansion coefficient) optimization
wec_values = [3.0] #[3.0, 2.6, 2.2, 1.8, 1.4, 1.0]
# wec_values = Float64[]
# for wec_value_ref = wec_values_ref
#     append!(wec_values, ones(n_wec_repeats)*wec_value_ref)
# end
n_wec = length(wec_values)
println(wec_values)

# intialize 
xopt = fill(zeros(1), 2*n_wec+1)
xopt[1] = deepcopy(x)
fopt = fill(0.0, 2*n_wec+1)
fopt[1] = deepcopy(initial_aep)
info = fill("", 2*n_wec+1)
info[1] = "Initial layout"

# set optimizer options
snopt_opt_short = Dict(
    "Major iterations limit" => 2,
    "Major optimality tolerance" => 1e-5,
    "Minor feasibility tolerance" => 1e-6,
    "Summary file" => "output/pareto/$(Int(max_cable_length))_meter_cable/opt_history_files/snopt_summary.out",
    "Print file" => "output/pareto/$(Int(max_cable_length))_meter_cable/opt_history_files/snopt_print.out"
)

snopt_opt_full = Dict(
    "Major iterations limit" => 2,
    "Major optimality tolerance" => 1e-5,
    "Minor feasibility tolerance" => 1e-6,
    "Summary file" => "output/pareto/$(Int(max_cable_length))_meter_cable/opt_history_files/snopt_summary.out",
    "Print file" => "output/pareto/$(Int(max_cable_length))_meter_cable/opt_history_files/snopt_print.out"
)


# ==================================================================================
# ========================== RUN AND TIME OPTIMIZATION =============================
# ==================================================================================

# first loop => no exclusions
# second loop => with exclusions
for j = 1:2

    # perform an optimization for each decreasing WEC value
    for i = 1:n_wec

        # results index
        k = (j-1)*n_wec + i + 1

        # get the cable network
        cable_length, _, cable_nodes = Cable_Analysis.cbl_analysis([xopt[i][1:nturbines]; substation_coordinates[1]], [xopt[i][nturbines+1:end]; substation_coordinates[2]], 1.0, return_network=true)
        params.cable_nodes[:] = cable_nodes
        println("Cable length = ", cable_length)
        cable_length_wrapper(x) = cable_length_wrapper(x, params)

        # set the WEC value in FlowFarm
        println("Running with WEC = ", wec_values[i])
        params.model_set.wake_deficit_model.wec_factor[1] = wec_values[i]
        
        # cold start (round 1)
        println("\n\tRound 1")
        solver_cold_short = SNOPT(options=snopt_opt_short)
        options_cold_short = Options(derivatives=UserDeriv(), sparsity=DensePattern(), solver=solver_cold_short)
        if j == 1
            xopt[k], _, _, out = minimize(wind_farm_opt, xopt[k-1], n_constraints, lb, ub, lb_g, ub_g, options_cold_short)
        else
            xopt[k], _, _, out = minimize(wind_farm_opt_with_exclusions, xopt[k-1], n_constraints_with_exclusions, lb, ub, lb_g_with_exclusions, ub_g_with_exclusions, options_cold_short)
        end

        # get the cable network
        cable_length, _, cable_nodes = Cable_Analysis.cbl_analysis([xopt[i][1:nturbines]; substation_coordinates[1]], [xopt[i][nturbines+1:end]; substation_coordinates[2]], 1.0, return_network=true)
        params.cable_nodes[:] = cable_nodes
        println("Cable length = ", cable_length)
        cable_length_wrapper(x) = cable_length_wrapper(x, params)

        # warm start (round 2)
        println("\n\tRound 2")
        solver_warm_short = SNOPT(warmstart=out.warm, options=snopt_opt_short)
        options_warm_short = Options(derivatives=UserDeriv(), sparsity=DensePattern(), solver=solver_warm_short)
        if j == 1
            xopt[k], _, _, out = minimize(wind_farm_opt, xopt[k], n_constraints, lb, ub, lb_g, ub_g, options_warm_short)
        else
            xopt[k], _, _, out = minimize(wind_farm_opt_with_exclusions, xopt[k], n_constraints_with_exclusions, lb, ub, lb_g_with_exclusions, ub_g_with_exclusions, options_warm_short)
        end

        # get the cable network
        cable_length, _, cable_nodes = Cable_Analysis.cbl_analysis([xopt[i][1:nturbines]; substation_coordinates[1]], [xopt[i][nturbines+1:end]; substation_coordinates[2]], 1.0, return_network=true)
        params.cable_nodes[:] = cable_nodes
        println("Cable length = ", cable_length)
        cable_length_wrapper(x) = cable_length_wrapper(x, params)

        # warm start (round 3)
        println("\tRound 3")
        # solver_warm_full = SNOPT(warmstart=out.warm, options=snopt_opt_full)
        solver_warm_full = SNOPT(options=snopt_opt_full)
        options_warm_full = Options(derivatives=UserDeriv(), sparsity=DensePattern(), solver=solver_warm_full)
        if j == 1
            xopt[k], fopt[k], info[k], out = minimize(wind_farm_opt, xopt[k-1], n_constraints, lb, ub, lb_g, ub_g, options_warm_full)
        else
            xopt[k], fopt[k], info[k], out = minimize(wind_farm_opt_with_exclusions, xopt[k-1], n_constraints_with_exclusions, lb, ub, lb_g_with_exclusions, ub_g_with_exclusions, options_warm_full)
        end

        # print optimization results
        println()
        println("major iter = ", out.major_iter)
        println("iterations = ", out.iterations)
        println("solve time = ", out.run_time)
        println("Info: ", info[k])
        println("End objective value: ", -fopt[k])
        # println("End cable constraint value: ", out.gstar)

        # plot intermediate boundary
        if plot_intermediate_layouts
            if j == 1
                plot_farm_winner_colome(intermediate_layout_figure_file, xopt[k], params; 
                    figure_number = k,
                    ext = [".png"],#, ".pdf"],
                    gray_only = true,
                    plot_houses = false,
                    plot_ponds = false,
                    plot_rivers = false,
                    plot_roads = false,
                    plot_misc = false,
                    save_fig = true
                    )
            else
                plot_farm_winner_colome(intermediate_layout_figure_file, xopt[k], params; 
                    figure_number = k,
                    ext = [".png"],#, ".pdf"],
                    gray_only = true,
                    plot_houses = true,
                    plot_ponds = true,
                    plot_rivers = true,
                    plot_roads = true,
                    plot_misc = true,
                    save_fig = true
                    )
            end
        end
    end
end

# get the cable network
cable_length, _, cable_nodes = Cable_Analysis.cbl_analysis([xopt[end][1:nturbines]; substation_coordinates[1]], [xopt[end][nturbines+1:end]; substation_coordinates[2]], 1.0, return_network=true)
params.cable_nodes[:] = cable_nodes


# ==================================================================================
# ========================== SHOW OPTIMIZATION RESULTS =============================
# ==================================================================================

# print optimization results
println("\n\n============ FINAL OPTIMIZATION RESULTS =============\n")
# println("Finished in : ", clkt, " (s)")
println("Info: ", info)
final_aep = -aep_wrapper(xopt[end])[1]
println("Final AEP: ", final_aep*1e-6/obj_scale, " MWh")
println()

# extract final turbine locations
turbine_x = copy(xopt[end][1:nturbines])
turbine_y = copy(xopt[end][nturbines+1:end])

# save final turbine locations to a text file
open(final_layout_file, "w") do io
    write(io, "# final turbine coordinates (x,y)\n")
    writedlm(io, [turbine_x turbine_y])
end

plot_farm_winner_colome(final_layout_figure_file, xopt[end], params; 
    figure_number = Nothing,
    ext = [".png", ".pdf"],
    gray_only = true,
    plot_houses = true,
    plot_ponds = true,
    plot_rivers = true,
    plot_roads = true,
    plot_misc = true,
    save_fig = true
    )

plot_farm_winner_colome(final_layout_figure_file_color, xopt[end], params; 
    figure_number = Nothing,
    ext = [".png", ".pdf"],
    gray_only = false,
    plot_houses = true,
    plot_ponds = true,
    plot_rivers = true,
    plot_roads = true,
    plot_misc = true,
    save_fig = true
    )
