#=
Plot the boundary for the WInner Colome wind farm boundary (including exclusions).
=#

using DelimitedFiles
using PyPlot
using PyCall
@pyimport matplotlib.patches as patch

function plot_winner_colome_boundary(; 
    gray_only = true,
    plot_houses = true,
    plot_ponds = true,
    plot_rivers = true,
    plot_roads = true,
    plot_misc = true,
    save_fig = false
    )

    # if gray_only
    #     fig = figure(figsize=(4.5, 2.5))
    #     fig.add_axes([0.2, 0.2, 0.7, 0.7])
    # else
    #     fig = figure(figsize=(6, 2.5))
    #     fig.add_axes([0.13, 0.2, 0.55, 0.7])
    # end
    fig = figure(figsize=(6, 2.5))
    fig.add_axes([0.13, 0.2, 0.55, 0.7])
    axis("square")
    xlabel("Easting (m)")
    ylabel("Northing (m)")


    # ======= MAIN BOUNDARY =======
    # get coordinates
    boundary_file_name = "input/boundary_files/boundary_winner_colome_3_utm.txt"
    boundary_vertices = readdlm(boundary_file_name, skipstart=1)
    boundary_x_closed_loop = [boundary_vertices[:,1]; boundary_vertices[1,1]]
    boundary_y_closed_loop = [boundary_vertices[:,2]; boundary_vertices[1,2]]

    # plot
    _xlim = [minimum(boundary_vertices[:,1]) - (maximum(boundary_vertices[:,1])-minimum(boundary_vertices[:,1]))/8, maximum(boundary_vertices[:,1]) + (maximum(boundary_vertices[:,1])-minimum(boundary_vertices[:,1]))/8]
    _ylim = [minimum(boundary_vertices[:,2]) - (maximum(boundary_vertices[:,2])-minimum(boundary_vertices[:,2]))/8, maximum(boundary_vertices[:,2]) + (maximum(boundary_vertices[:,2])-minimum(boundary_vertices[:,2]))/8]
    xlim(_xlim)
    ylim(_ylim)
    fill([_xlim[1], _xlim[2], _xlim[2], _xlim[1]], [_ylim[1], _ylim[1], _ylim[2], _ylim[2]],  
        color="silver"
        )
    fill(boundary_x_closed_loop, boundary_y_closed_loop,
        color="white"
        )


    # ======= HOUSES =======
    if plot_houses
        # get coordinates
        houses_file_name = "input/exclusions/houses/houses_utm.txt"
        houses = readdlm(houses_file_name, skipstart=1)

        # plot
        for i = 1:length(houses[:,1])
            house_x = houses[i,1]
            house_y = houses[i,2]
            plt.gcf().gca().add_artist(plt.Circle([house_x, house_y], 370.0, 
                fill=true, 
                color= gray_only ? "silver" : "tab:olive", 
                alpha= gray_only ? 1.0 : 0.2,
                linewidth=0.5
                )) 
        end
    end


    # ======= PONDS =======
    if plot_ponds
        # get coordinates
        ponds_file_name = "input/exclusions/ponds/ponds_utm.txt"
        ponds = readdlm(ponds_file_name, skipstart=1)

        # plot
        for i = 1:length(ponds[:,1])
            pond_x = ponds[i,1]
            pond_y = ponds[i,2]
            pond_r = ponds[i,3]
            plt.gcf().gca().add_artist(plt.Circle([pond_x, pond_y], pond_r, 
                fill=true, 
                color= gray_only ? "silver" : "tab:blue", 
                alpha= gray_only ? 1.0 : 0.2,
                linewidth=0.5
                )) 
        end
    end


    # ======== RIVERS =======
    if plot_rivers
        # get coordinates
        rivers_file_names = ["input/exclusions/rivers/river_1_utm.txt"]
        rivers = [readdlm(rivers_file_names[i], skipstart=1) for i=1:length(rivers_file_names)]

        # plot
        for river in rivers
            river_x_closed_loop = [river[:,1]; river[1,1]]
            river_y_closed_loop = [river[:,2]; river[1,2]]
            fill(river_x_closed_loop, river_y_closed_loop,
                color= gray_only ? "silver" : "tab:blue",
                alpha= gray_only ? 1.0 : 0.2,
                linewidth=0.5
                )
        end
    end


    # ======= ROADS =======
    if plot_roads
        # get coordinates
        roads_file_names = ["input/exclusions/roads/road_$(lpad(i,2,"0"))_utm.txt" for i=1:11]
        roads = [readdlm(roads_file_names[i], skipstart=1) for i=1:length(roads_file_names)]

        # plot
        for road in roads
            road_x_closed_loop = [road[:,1]; road[1,1]]
            road_y_closed_loop = [road[:,2]; road[1,2]]
            fill(road_x_closed_loop, road_y_closed_loop,
                color= gray_only ? "silver" : "tab:brown",
                alpha= gray_only ? 1.0 : 0.2,
                linewidth=0.5
                )
        end
    end


    # ======= MISC =======
    if plot_misc
        # misc coordinates
        misc_file_names = ["input/exclusions/misc/misc_utm.txt", "input/exclusions/misc/substation_utm.txt"]
        miscs = [readdlm(misc_file_names[i], skipstart=1) for i=1:length(misc_file_names)]

        # plot
        for misc in miscs
            for i = 1:length(misc[:,1])
                misc_x = misc[i,1]
                misc_y = misc[i,2]
                misc_r = misc[i,3]
                plt.gcf().gca().add_artist(plt.Circle([misc_x, misc_y], misc_r, 
                    fill=true, 
                    color= gray_only ? "silver" : "tab:orange", 
                    alpha= gray_only ? 1.0 : 0.2,
                    linewidth=0.5
                    )) 
            end
        end
    end
    


    # create legend
    if !gray_only
        legend_elements = [
            patch.Patch(facecolor="tab:olive", alpha=0.2, label="Residence"),
            patch.Patch(facecolor="tab:blue", alpha=0.2, label="River/pond"),
            patch.Patch(facecolor="tab:brown", alpha=0.2, label="Road"),
            patch.Patch(facecolor="tab:orange", alpha=0.2, label="Misc"),
            ]
        legend(handles=legend_elements, bbox_to_anchor=(1.05, 1), loc="upper left", ncol=1)
    end

    # save figure
    if save_fig
        if gray_only
            savefig("boundary_winner_colome_3_grey.png", dpi=400)
            savefig("boundary_winner_colome_3_grey.pdf")
        else
            savefig("boundary_winner_colome_3_color.png", dpi=400)
            savefig("boundary_winner_colome_3_color.pdf")
        end
    end

    if gray_only
        return []
    else
        return legend_elements
    end

end
