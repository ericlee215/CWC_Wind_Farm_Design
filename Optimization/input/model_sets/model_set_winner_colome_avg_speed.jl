import FLOWFarm; const ff = FLOWFarm
using DelimitedFiles 
import YAML

# get initial turbine x and y locations
turbine_coordinates = readdlm(initial_layout_file, skipstart=1)   # convert miles to meters
turbine_x = turbine_coordinates[:,1]
turbine_y = turbine_coordinates[:,2]

# calculate the number of turbines
nturbines = length(turbine_x)

# set turbine base heights
turbine_z = zeros(nturbines)

# set turbine yaw values
turbine_yaw = zeros(nturbines)

# set turbine design parameters
include("../../" * turbine_params_file)

# rotor swept area sample points (normalized by rotor radius)
rotor_points_y = [0.0]
rotor_points_z = [0.0]

# set flow parameters
wind_data = YAML.load_file(windrose_file)["definitions"]["wind_inflow"]["properties"]
winddirections = wind_data["direction"]["bins"]
windspeeds = wind_data["speed"]["bins"]
windprobabilities = wind_data["direction"]["frequency"]
ambient_ti = wind_data["turbulence_intensity"]["default"]
nstates = length(winddirections)
winddirections *= pi/180.0
air_density = 1.1716  # kg/m^3
shearexponent = 0.31
ambient_tis = zeros(nstates) .+ ambient_ti
measurementheight = zeros(nstates) .+ hub_height[1]

# load power curve
powerdata = readdlm(turbine_power_curve_file, skipstart=1)
pvelpoints = powerdata[:,1]
powerpoints = powerdata[:,2]*1E6

# initialize power model
if power_curve_with_cp_values
    power_model = ff.PowerModelCpPoints(pvelpoints, powerpoints)
else
    power_model = ff.PowerModelPowerPoints(pvelpoints, powerpoints)
end
power_models = Vector{typeof(power_model)}(undef, nturbines)
for i = 1:nturbines
    power_models[i] = power_model
end

# load thrust curve
ctdata = readdlm(turbine_thrust_curve_file, skipstart=1)
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
wakedeficitmodel = ff.GaussYawVariableSpread() # ff.GaussYaw()
wakedeflectionmodel = ff.GaussYawVariableSpreadDeflection()
wakecombinationmodel = ff.LinearLocalVelocitySuperposition()
localtimodel = ff.LocalTIModelNoLocalTI() # ff.LocalTIModelNoLocalTI() LocalTIModelMaxTI

# initialize model set
model_set = ff.WindFarmModelSet(wakedeficitmodel, wakedeflectionmodel, wakecombinationmodel, localtimodel)