using DelimitedFiles
using PyPlot

# layout_number_start = 201
# layout_number_end = 300

# final_aeps = zeros(layout_number_end-layout_number_start+1)
# for i = layout_number_start:layout_number_end

#     if isfile("output/final_layouts/winner_colome_$(lpad(i,3,"0"))_aep.txt")
#         final_aeps[i-layout_number_start+1] = readdlm("output/final_layouts/winner_colome_$(lpad(i,3,"0"))_aep.txt", skipstart=1)[1]
        
#         # save final turbine locations to a csv file
#         turbine_coordinates = readdlm("output/final_layouts/winner_colome_$(lpad(i,3,"0")).txt", skipstart=1)
#         open("output/final_layouts/winner_colome_$(lpad(i,3,"0")).csv", "w") do io
#             write(io, "# Winner Colome wind farm turbine coordinates (x,y)\n")
#             writedlm(io, turbine_coordinates, ',')
#         end

#     else
#         final_aeps[i-layout_number_start+1]
#     end

# end

# ranked_layouts = reverse(sortperm(final_aeps))
# display(final_aeps[ranked_layouts])


layout_numbers = [207, 217, 227, 282, 294, 295, 300]

for i in layout_numbers

    data = readdlm("output/final_layouts/winner_colome_$(lpad(i,3,"0")).txt", skipstart=1)
    turbine_x = data[:,1]
    turbine_y = data[:,2]

    cost, cables = ff.Cable_Analysis.cbl_analysis(turbine_x,turbine_y,1.0,return_network=true,plot_tree=false)
    
    boundary_vertices = readdlm("input/boundary_files/boundary_winner_colome_2.txt", skipstart=1)
    clf()
    axis("square")
    # plot([boundary_vertices[:,1]; boundary_vertices[1,1]], [boundary_vertices[:,2]; boundary_vertices[1,2]])
    xlim(minimum(boundary_vertices[:,1]) - (maximum(boundary_vertices[:,1])-minimum(boundary_vertices[:,1]))/5, maximum(boundary_vertices[:,1]) + (maximum(boundary_vertices[:,1])-minimum(boundary_vertices[:,1]))/5)
    ylim(minimum(boundary_vertices[:,2]) - (maximum(boundary_vertices[:,2])-minimum(boundary_vertices[:,2]))/5, maximum(boundary_vertices[:,2]) + (maximum(boundary_vertices[:,2])-minimum(boundary_vertices[:,2]))/5)

    for i = 1:length(cables)
        x = [cables[i][1][1], cables[i][2][1]]
        y = [cables[i][1][2], cables[i][2][2]]
        plot(x,y,"g-",alpha=0.5)
    end

    # add final turbine locations to plot
    for i = 1:length(turbine_x)
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), 110.0/2.0, fill=true, color="C1")) 
    end

    savefig("output/final_layout_figures/best_layouts_boundary_2/winner_colome_$(lpad(i,3,"0"))_cables.png", transparent=true, dpi=400)

end





include("../cableCost/cbl_analysis2.jl")

data = readdlm("output/final_layouts/winner_colome_$(lpad(207,3,"0")).txt", skipstart=1)
turbine_x = data[:,1]
turbine_y = data[:,2]

cost, cables, nodes = Cable_Analysis.cbl_analysis(turbine_x, turbine_y, 1.0, return_network=true)

function get_cable_length(x, nodes)

    n_turbines = Int(length(x)/2)

    turbine_x = x[1:n_turbines]
    turbine_y = x[n_turbines+1:end]

    cable_length = 0

    for i = 1:length(nodes)
        cable_length += sqrt((turbine_x[nodes[i][1]] - turbine_x[nodes[i][2]])^2 + (turbine_y[nodes[i][1]] - turbine_y[nodes[i][2]])^2)
    end

    return cable_length
end
