# design parameters for the NREL 5MW offshore reference turbine
rotor_diameter = zeros(nturbines) .+ 126.0 # m
hub_height = zeros(nturbines) .+ 90.0   # m
cut_in_speed = zeros(nturbines) .+ 3.  # m/s
cut_out_speed = zeros(nturbines) .+ 25.  # m/s
rated_speed = zeros(nturbines) .+ 11.4  # m/s
rated_power = zeros(nturbines) .+ 5.0E6  # W
generator_efficiency = zeros(nturbines) .+ 0.944

# power curve file path
turbine_power_curve_file = "input/turbine_files/NREL_5MW/NREL_5MW_power_coefficient_curve.txt"
power_curve_with_cp_values = true

# thrust curve file path
turbine_thrust_curve_file = "input/turbine_files/NREL_5MW/NREL_5MW_thrust_coefficient_curve.txt"