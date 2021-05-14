import FLOWFarm; const ff=FLOWFarm

function problem_setup_simple_example_1(; initial_layout_number=1, generate_random_layout=false)

    #### SET PROBLEM PARAMETERS FROM INPUT YAML FILES ####
    
    path_to_input_files_directory = "input/yaml_files/simple_example_1/"

    # constraints (spacing, boundary, cable length)
    farm_constraints = ff.WindFarmConstraints(path_to_input_files_directory * "constraints_simple_example_1.yaml", path_to_input_files_directory=path_to_input_files_directory)

    # farm design
    farm_definition = ff.WindFarm(path_to_input_files_directory * "farm_definition_simple_example_1.yaml", path_to_input_files_directory=path_to_input_files_directory)

    if generate_random_layout
        # generate and save a random starting layout
        initial_layout_number = lpad(initial_layout_number,3,"0")
        farm_definition = ff.generate_random_layout(farm_definition, farm_constraints.boundaries[1])
    end
    
    # wind resource
    wind_resource = ff.DiscretizedWindResource(path_to_input_files_directory * "wind_resource_winner_colome_24dirs_avgspeeds.yaml", average_speeds=true)
    
    # farm states
    nstates = length(wind_resource.wind_probabilities)
    nturbines = length(farm_definition.turbine_x)
    farm_states = Vector{ff.AbstractWindFarmModel}(undef, nstates)
    for i = 1:nstates
        farm_states[i] = ff.SingleWindFarmState(i, [zeros(nturbines) for i=1:9]..., zeros(Int64, nturbines))
    end

    # farm problem
    farm_problem = ff.WindFarmProblem(farm_definition, wind_resource, farm_states)

    # flow models
    flow_models = ff.WindFarmModelSet(path_to_input_files_directory * "flow_models.yaml")

    return farm_problem, flow_models, farm_constraints
end
