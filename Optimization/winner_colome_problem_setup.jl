import FLOWFarm; const ff=FLOWFarm
import YAML

function problem_setup_winner_colome(; initial_layout_number=1, generate_random_layout=false)

    #### SET PROBLEM PARAMETERS FROM INPUT YAML FILES ####
    
    path_to_input_files_directory = "input/yaml_files/winner_colome/"

    # constraints (spacing, boundary, cable length)
    farm_constraints = ff.WindFarmConstraints(path_to_input_files_directory * "constraints_winner_colome.yaml", path_to_input_files_directory=path_to_input_files_directory)

    # farm design
    farm_definition = ff.WindFarm(path_to_input_files_directory * "farm_initial_design_winner_colome_layout000.yaml", path_to_input_files_directory=path_to_input_files_directory)

    if generate_random_layout
        # generate and save a random starting layout
        initial_layout_number = lpad(initial_layout_number,3,"0")
        farm_definition = ff.generate_random_layout(farm_definition, farm_constraints.boundaries[1])
        data = YAML.load(open(path_to_input_files_directory * "farm_initial_design_winner_colome_layout000.yaml"))
        data["definitions"]["position"]["items"] = [[farm_definition.turbine_x[i], farm_definition.turbine_y[i]] for i=1:length(farm_definition.turbine_x)]
        YAML.write_file(path_to_input_files_directory * "farm_initial_design_winner_colome_layout$initial_layout_number.yaml", data)
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
