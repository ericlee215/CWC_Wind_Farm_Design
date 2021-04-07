using PyPlot

include("input/pareto/pareto_functions.jl")
include("input/pareto/simple_test_functions.jl")

#######################################################################
# SET OPTIMIZATION ALGORITHM PARAMETERS:
# optimization algorithm, optimality tolerance, constraint tolerance, max function evaluations
#######################################################################
opt_algorithm = :GN_ISRES   # GN_ISRES ("G" = global, "N" = gradient-free, "ISRES" = Improved Stochastic Ranking Evolution Strategy)
opt_tolerance = 1e-6
con_tolerance = 1e-8
max_eval = Int(1e6)

opt_algorithm_param = NLoptAlgorithmParam(opt_algorithm, opt_tolerance, con_tolerance, max_eval)  


#######################################################################
# SET PARAMETERS FOR GENERATING THE PARETO FRONT
# objective 1, objective 2, epsilon
#######################################################################
objective_1 = rosenbrock_ndim
objective_2 = quadratic_ndim
epsilon = 0.1:0.1:3.0

pareto_param = ParetoFrontEpsilonConstraintParam(rosenbrock_ndim, quadratic_ndim, 0.1:0.1:3.0)


#######################################################################
# SET DESIGN VARIABLE PARAMETERS FOR OPTIMIZATION
# starting points, lower bounds, upper bounds
#######################################################################
starting_points = [[ 1.0,  2.0,  3.0],
                   [-1.0, -2.5, -4.0],
                   [ 6.0,  7.0,  8.0],
                   [-6.0, -7.5, -9.0]]  # four different starting points for each epsilon value
lower_bound = [-10.0, -10.0, -10.0]
upper_bound = [10.0, 10.0, 10.0]

design_variables_param = OptDesignVariablesParam(starting_points, lower_bound, upper_bound)


#######################################################################
# GET PARETO FRONT POINTS
#######################################################################

fopt, xopt, info = generate_2D_pareto_points(design_variables_param, pareto_param, opt_algorithm_param)


#######################################################################
# PLOT PARETO FRONT
#######################################################################
figure()
n_epsilon = length(epsilon)
n_samples = length(starting_points)
for i = 1:n_epsilon
    scatter(pareto_param.epsilon[i]*ones(n_samples), fopt[i,:], color="C0")
end
title("Example Pareto Front")
xlabel("Quadratic Function")
ylabel("Rosenbrock Function")
savefig("output/pareto/pareto1_example.png")