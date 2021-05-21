
"""
NLoptAlgorithmParam(opt_algorithm, opt_tolerance, con_tolerance)

Container to hold algorithm parameters for the NLopt optimization package.

# Arguments
- `opt_algorithm::Symbol`: Symbol specifying the optimization algorithm, according to the NLopt.jl documentation (https://github.com/JuliaOpt/NLopt.jl#the-opt-type)
- `opt_tolerance::Float64`: Optimality tolerance
- `con_tolerance::Float64`: Constraint tolerance
"""
struct NLoptAlgorithmParam
    opt_algorithm
    opt_tolerance
    con_tolerance
    max_eval
end
NLoptAlgorithmParam() = NLoptAlgorithmParam(:GN_ISRES, 1e-6, 1e-8, Int(1e6))
NLoptAlgorithmParam(opt_algorithm) = NLoptAlgorithmParam(opt_algorithm, 1e-6, 1e-8, Int(1e6))

struct SnoptWECAlgorithm
    max_eval
    wec
    opt_tolerance
    checkgradients
    parallel_processing
end
SnoptWECAlgorithm() = SnoptWECAlgorithm(1e6, [3.0, 2.6, 2.2, 1.8, 1.4, 1.0, 1.0], [1e-3, 1e-3, 1e-4, 1e-4, 1e-5, 1e-6, 1e-6], true, false)


function set_optimization_problem(optimizer_param::NLoptAlgorithmParam, farm_problem::ff.AbstractWindFarmProblem, flow_models::ff.AbstractModelSet, farm_constraints::ff.AbstractWindFarmConstraints; return_unscaled_obj=true)
    
    function obj!(x, g, obj_scale, farm_problem::ff.AbstractWindFarmProblem, flow_models::ff.AbstractModelSet)

        # get number of turbines
        nturbines = Int(length(x)/2)

        # extract x and y locations of turbines from design variables vector
        farm_problem.farm_description.turbine_x[:] = x[1:nturbines]
        farm_problem.farm_description.turbine_y[:] = x[nturbines+1:end]

        # calculate AEP
        AEP = obj_scale*ff.calculate_aep!(farm_problem, flow_models)

        return -AEP
    end

    # calculate AEP for the initial turbine layout
    initial_x = deepcopy([farm_problem.farm_description.turbine_x; farm_problem.farm_description.turbine_y])
    initial_aep = -obj!(initial_x, [], 1.0, farm_problem, flow_models)

    # report farm energy and power calculation for the initial design
    println("\n#### INITIAL FARM PERFORMANCE ####")
    println("\nInitial AEP: ", initial_aep*1e-6, " MWh")
    println("Average farm power: ", initial_aep*1e-6/(365*24), " MW")

    # determine an appropriate scaling factor for the objective function 
    obj_scale = 10^round(log10(1/initial_aep))

    # create objective function wrapper for the optimizer
    obj!(x, g) = obj!(x, g, obj_scale, farm_problem, flow_models)

    function con!(c, x, g, con_scale, farm_problem::ff.AbstractWindFarmProblem, farm_constraints::ff.AbstractWindFarmConstraints)

        # get number of turbines
        nturbines = Int(length(x)/2)

        # extract x and y locations of turbines from design variables vector
        turbine_x = x[1:nturbines]
        turbine_y = x[nturbines+1:end]

        # extract other parameters
        boundaries = farm_constraints.boundaries
        max_rotor_diameter = maximum([farm_problem.farm_description.turbine_definitions[i].rotor_diameter[i] for i=1:length(farm_problem.farm_description.turbine_definitions)])
        min_turbine_spacing = farm_constraints.minimum_turbine_spacing
        max_cable_length = farm_constraints.maximum_cable_length

        # boundary constraint
        nboundaries = length(boundaries)
        #TODO initialize the boundary constraints vector with the correct length based on whether the boundary region is convex or concave
        boundary_con = zeros(2*nturbines*nboundaries)
        for i = 1:nboundaries
            boundary_con = boundary_wrapper(turbine_x, turbine_y, boundaries[i])
        end

        # turbine spacing constraint
        spacing_con = spacing_wrapper(turbine_x, turbine_y, max_rotor_diameter, min_turbine_spacing)

        # cable length constraint
        cable_length_con = cable_length_wrapper(turbine_x, turbine_y, max_cable_length)

        c[:] = [boundary_con; spacing_con; cable_length_con]
    end

    function boundary_wrapper(turbine_x, turbine_y, boundary::ff.PolygonBoundary)
        # if boundary.isconvex == true
        #     return convex_boundary(boundary.vertices, boundary.normals, turbine_x, turbine_y)
        # else
            return ff.ray_cast_boundary(boundary.vertices, boundary.normals, turbine_x, turbine_y)
        # end
    end

    function spacing_wrapper(turbine_x, turbine_y, rotor_diameter, min_turbine_spacing)
        rotor_diameter*min_turbine_spacing .- ff.turbine_spacing(turbine_x, turbine_y)
    end

    function cable_length_wrapper(turbine_x, turbine_y, max_cable_length)
        Cable_Analysis.cbl_analysis(turbine_x, turbine_y, 1.0) - max_cable_length
    end

    con!(c, x, g) = con!(c, x, g, [], farm_problem, farm_constraints)

    # get number of turbines
    nturbines = length(farm_problem.farm_description.turbine_x)

    # get bounds for the design variables
    turbine_x_lower_bound = minimum([minimum(farm_constraints.boundaries[i].vertices[:,1]) for i=1:length(farm_constraints.boundaries)])
    turbine_y_lower_bound = minimum([minimum(farm_constraints.boundaries[i].vertices[:,2]) for i=1:length(farm_constraints.boundaries)])
    turbine_x_upper_bound = maximum([maximum(farm_constraints.boundaries[i].vertices[:,1]) for i=1:length(farm_constraints.boundaries)])
    turbine_y_upper_bound = maximum([maximum(farm_constraints.boundaries[i].vertices[:,2]) for i=1:length(farm_constraints.boundaries)])
    lb = [zeros(nturbines) .+ turbine_x_lower_bound; zeros(nturbines) .+ turbine_y_lower_bound]
    ub = [zeros(nturbines) .+ turbine_x_upper_bound; zeros(nturbines) .+ turbine_y_upper_bound]

    # get number of constraints
    nconstraints = Int(nturbines + (nturbines-1)*nturbines/2 + 1)

    # initial cable length constraint value
    initial_con = zeros(nconstraints)
    con!(initial_con, initial_x, [], [], farm_problem, farm_constraints)

    println("Initial cable length constraint value: ", initial_con[end])

    # set optimizer parameters
    opt = Opt(optimizer_param.opt_algorithm, nturbines*2)
    opt.lower_bounds = lb
    opt.upper_bounds = ub
    opt.maxeval = optimizer_param.max_eval
    opt.xtol_rel = optimizer_param.opt_tolerance
    con_tolerance = zeros(nconstraints) .+ optimizer_param.con_tolerance

    # set optimization objective and constraints
    opt.min_objective = obj!
    inequality_constraint!(opt, con!, con_tolerance)

    if return_unscaled_obj
        # define unscaled objective function
        obj_unscaled!(x, g) = obj!(x, g)/obj_scale
        return opt, obj!, con!, obj_unscaled!
    else
        return opt, obj!, con!
    end
end