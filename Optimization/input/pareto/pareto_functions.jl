using NLopt, FillArrays

"""
    ParetoFrontEpsilonConstraintParam(objective_1, objective_2, epsilon)

Container to hold parameters and functions to create a pareto front using the epsilon-constraint method.

# Arguments
- `objective_1::Function`: First objective function. Takes in a vector of design variables, returns the objective value.
- `objective_2::Function`: Second objective function. Takes in a vector of design variables, returns the objective value. (This will be used to constrain the first objective.)
- `epsilon::Float64`
"""
struct ParetoFrontEpsilonConstraintParam
    objective_1
    objective_2
    epsilon
end

"""
    OptDesignVariablesParam(x, lower_bound, upper_bound)

Container to hold parameters for the design variable for an optimization.

# Arguments
- `x::Array{Array{Float64,1}}`: Array of design variable vectors. These are the initial starting points for the optimization(s).
- `lower_bound::Array{Float64,1}`: Vector with the lower bounds for all the design variables.
- `upper_bound::Array{Float64,1}`: Vector with the upper bounds for all the design variables.
"""
struct OptDesignVariablesParam
    x::Array{Array{Float64,1},1}
    lower_bound
    upper_bound
end

"""
    generate_2D_pareto_points(design_variables_param::OptDesignVariablesParam, pareto_param::ParetoFrontEpsilonConstraintParam, opt_algorithm_param::NLoptAlgorithmParam)

Generates the 2D Pareto front.

# Arguments
- `design_variables_param::OptDesignVariablesParam`: Container to hold parameters for the design variable for an optimization.
- `pareto_param::ParetoFrontEpsilonConstraintParam`: Container to hold parameters and functions to create a pareto front using the epsilon-constraint method.
- `opt_algorithm_param::NLoptAlgorithmParam`: Container to hold algorithm parameters for the NLopt optimization package.
"""
function generate_2D_pareto_points(design_variables_param::OptDesignVariablesParam, pareto_param::ParetoFrontEpsilonConstraintParam, opt_algorithm_param::NLoptAlgorithmParam; non_pareto_con! = [], non_pareto_con_tolerance=[])

    # get number of epsilon values, starting point samples, and design variables
    n_epsilon = length(pareto_param.epsilon)
    n_samples = length(design_variables_param.x)
    n_designvariables = length(design_variables_param.x[1])

    # define objective and constraint functions for NLopt optimizer
    obj!(x, gradient) = pareto_param.objective_1(x)
    con_epsilon!(x, gradient, epsilon) = pareto_param.objective_2(x) - epsilon   # to satisfy the constraint, the value of objective 2 must be below epsilon

    # initialize arrays to hold ouput values
    fopt = zeros(n_epsilon, n_samples)
    xopt = zeros(n_epsilon, n_samples, n_designvariables)
    info = fill(:EmptySymbol, n_epsilon, n_samples)

    # get the pareto front points
    for i = 1:n_epsilon
        epsilon = pareto_param.epsilon[i]
        # create NLopt optimizer object
        opt = Opt(opt_algorithm_param.opt_algorithm, n_designvariables)
        # pass objective and constraint functions into optimization object
        opt.min_objective = obj!
        inequality_constraint!(opt, (x,g) -> con_epsilon!(x,g,epsilon), opt_algorithm_param.con_tolerance)
        if non_pareto_con! != []
            inequality_constraint!(opt, non_pareto_con!, non_pareto_con_tolerance)
        end
        # set optimization tolerance, bounds for design variables, and maximum allowed function evaluations
        opt.xtol_rel = opt_algorithm_param.opt_tolerance
        opt.lower_bounds = design_variables_param.lower_bound
        opt.upper_bounds = design_variables_param.upper_bound
        opt.maxeval = opt_algorithm_param.max_eval
        # run optimizations
        for j = 1:n_samples
            x = design_variables_param.x[j]
            (fopt[i,j], xopt[i,j,:], info[i,j]) = optimize(opt, x)
        end
    end
    
    return fopt, xopt, info
end