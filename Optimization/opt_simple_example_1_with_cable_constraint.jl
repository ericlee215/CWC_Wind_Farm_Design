using NLopt

include("problem_setup_simple_example_1.jl")
include("../cableCost/cbl_analysis.jl")
include("optimization_functions.jl")

# set up farm
farm_problem, flow_models, farm_constraints = problem_setup_simple_example_1(generate_random_layout=true)

# set optimizer parameters
optimizer_param = NLoptAlgorithmParam(:GN_ISRES, 1e-5, 1e-6, 100)

# set up optimization problem
opt, obj!, con!, obj_unscaled! = set_optimization_problem(optimizer_param, farm_problem, flow_models, farm_constraints, return_unscaled_obj=true)

# get initial turbine coordinates
initial_x = deepcopy([farm_problem.farm_description.turbine_x; farm_problem.farm_description.turbine_y])

# get number of turbines and number of constraints
nturbines = length(farm_problem.farm_description.turbine_x)
nconstraints = Int(nturbines + (nturbines-1)*nturbines/2 + 1)

# run optimization
clk1 = time()
(fopt, xopt, exit_code) = optimize(opt, initial_x)
clk2 = time()

numevals = opt.numevals # the number of function evaluations
final_aep = -obj_unscaled!(xopt, [])
println("\n#### OPTIMIZED FARM PERFORMANCE ####")
println("\nFinal AEP: ", final_aep*1e-6, " MWh")
println("Average farm power: ", final_aep*1e-6/(365*24), " MW")
println("Final cable length constraint value: ", con!(zeros(nconstraints), xopt, [])[end])

println("\n#### OPTIMIZATION INFO ####")
println("\nNumber of AEP function evaulations: ", numevals)
println("Exit code: ", exit_code)
println("Total run time: ", round(clk2-clk1, digits=5), " s")
println("Run time / number of function evaluations: ", (clk2-clk1)/numevals)

println(xopt)

Cable_Analysis.cbl_analysis(xopt[1:nturbines], xopt[nturbines+1:end], 1.0, plot_tree=true)
