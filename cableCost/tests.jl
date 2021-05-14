# a function to test the Cable_Analysis package 

include("cbl_analysis.jl")

x_cor = [330.314, 486.996, 966.240, 655.232, 490.497, 1163.653, 997.640, 301.532, 838.581]
y_cor = [391.361, 589.599, 1202.957, 443.494, 914.822, 1099.077, 329.300, 757.712, 1067.241]
linear_cost = 1

cost = Cable_Analysis.cbl_analysis(x_cor,y_cor,linear_cost)

print(cost)