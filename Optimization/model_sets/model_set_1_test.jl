import FlowFarm; const ff = FlowFarm
using DelimitedFiles 

# set initial turbine x and y locations
turbine_x = [-3.0, 0.0, 3.0, 0.0, 0.0, -1.5, 0.0, 1.5, 0.0].*150.0
turbine_y = [0.0, 3.0, 0.0, -3.0, 0.0, 0.0, 1.5, 0.0, -1.5].*150.0

# calculate the number of turbines
nturbines = length(turbine_x)

# set turbine base heights
turbine_z = zeros(nturbines)

# set turbine yaw values
turbine_yaw = zeros(nturbines)

# set turbine design parameters
rotor_diameter = zeros(nturbines) .+ 80.0 # m
hub_height = zeros(nturbines) .+ 70.0   # m
cut_in_speed = zeros(nturbines) .+4.  # m/s
cut_out_speed = zeros(nturbines) .+25.  # m/s
rated_speed = zeros(nturbines) .+16.  # m/s
rated_power = zeros(nturbines) .+2.0E6  # W
generator_efficiency = zeros(nturbines) .+0.944

# rotor swept area sample points (normalized by rotor radius)
rotor_points_y = [0.0]
rotor_points_z = [0.0]

# set flow parameters
wind_speed = 8.0
air_density = 1.1716  # kg/m^3
ambient_ti = 0.077
shearexponent = 0.15
winddirections = [275.0*pi/180.0, 0.0, pi]
windspeeds = [wind_speed, wind_speed, wind_speed]
windprobabilities = [1.0/3.0,1.0/3.0,1.0/3.0]
ambient_tis = [ambient_ti, ambient_ti, ambient_ti]
measurementheight = [hub_height[1], hub_height[1], hub_height[1]]

# load power curve
powerdata = readdlm("input_files/niayifar_vestas_v80_power_curve_observed.txt",  ',', skipstart=1)
pvelpoints = powerdata[:,1]
powerpoints = powerdata[:,2]*1E6

# initialize power model
power_model = ff.PowerModelPowerPoints(pvelpoints, powerpoints)
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# load thrust curve
ctdata = readdlm("input_files/predicted_ct_vestas_v80_niayifar2016.txt",  ',', skipstart=1)
ctvelpoints = ctdata[:,1]
ctpoints = ctdata[:,2]

# initialize thrust model
ct_model = ff.ThrustModelCtPoints(ctvelpoints, ctpoints)
ct_models = Vector{typeof(ct_model)}(undef, nturbines)
for i = 1:nturbines
    ct_models[i] = ct_model
end

# initialize wind shear model
wind_shear_model = ff.PowerLawWindShear(shearexponent)

# get sorted indecies 
sorted_turbine_index = sortperm(turbine_x)

# initialize the wind resource definition
windresource = ff.DiscretizedWindResource(winddirections, windspeeds, windprobabilities, measurementheight, air_density, ambient_tis, wind_shear_model)

# set up wake and related models
wakedeficitmodel = ff.GaussYaw()
wakedeflectionmodel = ff.GaussYawDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelNoLocalTI()

# initialize model set
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)