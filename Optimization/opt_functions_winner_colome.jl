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
    return 5.0*rotor_diameter[1] .- ff.turbine_spacing(turbine_x,turbine_y)
end


# cable length constraint wrapper
function cable_length_wrapper(x, params)

    # include relevant globals
    max_cable_length = params.max_cable_length
    cable_nodes = params.cable_nodes
    substation_coordinates = params.substation_coordinates

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # append the substation as a node in the cable network
    append!(turbine_x, substation_coordinates[1])
    append!(turbine_y, substation_coordinates[2])

    cable_length = zeros(typeof(turbine_x[1]), 1)

    for i = 1:length(cable_nodes)
        cable_length[1] += sqrt((turbine_x[cable_nodes[i][1]] - turbine_x[cable_nodes[i][2]])^2 + (turbine_y[cable_nodes[i][1]] - turbine_y[cable_nodes[i][2]])^2)
    end

    return [cable_length[1] - max_cable_length]
end


# restriced areas (exclusions) constraint wrapper
function turbine_exclusion_wrapper(x, params)

    # include relevant globals
    turbine_circle_exclusions = params.turbine_circle_exclusions
    turbine_polygon_exclusions = params.turbine_polygon_exclusions

    # get number of turbines
    nturbines = Int(length(x)/2)

    # extract x and y locations of turbines from design variables vector
    turbine_x = x[1:nturbines]
    turbine_y = x[nturbines+1:end]

    # initialize constraint vector
    c = zeros(typeof(turbine_x[1]), (length(turbine_circle_exclusions) + length(turbine_polygon_exclusions)) * nturbines)

    # circle exclusions
    k = 1    
    for i in turbine_circle_exclusions
        # get and return boundary distances
        c[k:k+nturbines-1] = -ff.circle_boundary(i.center, i.radius, turbine_x, turbine_y)*1e-6
        k += nturbines
    end

    # polygon exclusions
    for i in turbine_polygon_exclusions
        # get and return boundary distances
        c[k:k+nturbines-1] = -ff.ray_trace_boundary(i.vertices, i.normals, turbine_x, turbine_y)
        k += nturbines
    end

    return c
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
    return [-AEP]
end


# ==================================================================================
# ====================== FUNCTIONS FOR SNOPT (OPTIMIZER) ===========================
# ==================================================================================

function wind_farm_opt(g, df, dg, x)

    # objective
    f = aep_wrapper(x)[1]
    df[:] = ForwardDiff.jacobian(aep_wrapper,x)

    # calculate spacing constraint value
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # calculate cable length constraint
    cable_length_con = cable_length_wrapper(x)
    dcl_dx = ForwardDiff.jacobian(cable_length_wrapper, x)

    # combine constraint values and jacobians into overall constaint value and jacobian arrays
    g[:] = [spacing_con; boundary_con; cable_length_con]
    dg[:] = [ds_dx; db_dx; dcl_dx]

    return f

end


function wind_farm_opt_with_exclusions(g, df, dg, x)

    # objective
    f = aep_wrapper(x)[1]
    df[:] = ForwardDiff.jacobian(aep_wrapper,x)

    # calculate spacing constraint value
    spacing_con = spacing_wrapper(x)
    ds_dx = ForwardDiff.jacobian(spacing_wrapper, x)

    # calculate boundary constraint
    boundary_con = boundary_wrapper(x)
    db_dx = ForwardDiff.jacobian(boundary_wrapper, x)

    # calculate cable length constraint
    cable_length_con = cable_length_wrapper(x)
    dcl_dx = ForwardDiff.jacobian(cable_length_wrapper, x)

    # calculate turbine exclusion constraint
    exclusion_con = turbine_exclusion_wrapper(x)
    de_dx = ForwardDiff.jacobian(turbine_exclusion_wrapper, x)

    # combine constraint values and jacobians into overall constaint value and jacobian arrays
    g[:] = [spacing_con; boundary_con; cable_length_con; exclusion_con]
    dg[:] = [ds_dx; db_dx; dcl_dx; de_dx]

    return f

end
