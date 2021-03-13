# design parameters for the Vestas V80 2MW turbine
rotor_diameter = zeros(nturbines) .+ 80.0 # m
hub_height = zeros(nturbines) .+ 70.0   # m
cut_in_speed = zeros(nturbines) .+ 4.  # m/s
cut_out_speed = zeros(nturbines) .+ 25.  # m/s
rated_speed = zeros(nturbines) .+ 16.  # m/s
rated_power = zeros(nturbines) .+ 2.0E6  # W
generator_efficiency = zeros(nturbines) .+ 0.944

# power curve file path
turbine_power_curve_file = "input/turbine_files/vestas_V80_2MW/vestas_v80_2MW_power_curve.txt"
power_curve_with_cp_values = false

# thrust curve file path
turbine_thrust_curve_file = "input/turbine_files/vestas_V80_2MW/vestas_v80_2MW_thrust_coefficient_curve.txt"