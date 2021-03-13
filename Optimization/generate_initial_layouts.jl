#=
AUTHOR: Wesley Holt
PURPOSE: To generate intial turbine layouts for circular boundary farms.
=#

import FlowFarm; const ff = FlowFarm
using PyPlot
using DelimitedFiles


# =========================================================================================
# ============================= SET UP PARAMETERS =========================================
# =========================================================================================

# ================= This block contains all the user inputs. ==============================

# generate layout file, figure, or both
generate_layout_file = true
generate_layout_figure = true

# set up farm parameters
nturbines = 9   # number of turbines
boundary_center = [0.0 0.0] # x,y coordinates for the center of the farm
boundary_radius = 600.0
rotor_diameter = 80.0
layout_directory = "input/initial_layouts/"
figure_directory = "input/initial_layout_figures/"

# set how many layouts to generate/plot
start_layout_number = 1
end_layout_number = 10
layout_numbers = range(start_layout_number, stop=end_layout_number, step=1)

# =========================================================================================
# =========================================================================================
# =========================================================================================

function generate_random_layout_circle_boundary(layout_filename, nturbines, rotor_diameter, boundary_center, boundary_radius; layout_directory="", min_spacing=1.0, layout_number=0)

    # ----- SET UP -----
    # set bounding box
    xrange = [boundary_center[1] - boundary_radius; boundary_center[1] + boundary_radius]
    yrange = [boundary_center[2] - boundary_radius; boundary_center[2] + boundary_radius]

    # scale the minimum spacing
    min_spacing *= rotor_diameter

    # ----- GENERATE RANDOM WIND FARM LAYOUT -----
    # initialize array to hold turbine coordinates
    locations = zeros(Float64, nturbines, 2)
    # set counter to prevent while-loop from running forever
    count = 0
    for i = 1:nturbines
        good_point = false
        while !good_point && count < 1e8
            # set good_point initially to true
            good_point = true   # this will be set to false if the generated point is not feasible
            # print out the number of loops every 1e6 loops
            if mod(count,1e6)==0
                println("loop count: ", count)
            end
            # update the counter
            count += 1
            # generate random point in the bounding box
            locations[i,:] = (rand(1,2) .- 0.5).*[xrange[2]-xrange[1] yrange[2]-yrange[1]] .+ boundary_center
            # calculate the distance(s) from the point to the boundary(s)
            distances_to_boundaries = ff.circle_boundary(boundary_center, boundary_radius, locations[1:i,1], locations[1:i,2])
            # determine if the point is inside the wind farm boundary (+ distance is outside, - distance is inside)
            for j = 1:length(distances_to_boundaries)
                if distances_to_boundaries[j] > 0.0
                    good_point = false
                end
            end
            # determine if the point is far enough away from other points
            for turb = 1:i-1
                spacing = sqrt((locations[turb,1]-locations[i,1])^2 + (locations[turb,2]-locations[i,2])^2)
                if spacing < min_spacing
                    good_point = false
                end
            end
        end
    end

    # ----- WRITE TURBINE COORDINATES TO A FILE -----
    # store x and y corrdinates into separate vectors
    turbine_x = locations[:,1]
    turbine_y = locations[:,2]
    
    # # write turbine coordinates to a YAML file
    # ff.write_turb_loc_YAML(layout_directory * layout_filename, turbine_x, turbine_y; 
    # title="Randomly generated $nturbines-turbine layout for a circular boundary wind farm. Layout: $layout_number", 
    # titledescription="Contains randomly generated coordinates for $nturbines turbines arranged in a circular boundary wind farm with a boundary radius of $boundary_radius m. The farm is centered at the coordinate ($(boundary_center[1]), $(boundary_center[2])). Each turbine has a rotor diameter of $rotor_diameter m.", 
    # turbinefile="", locunits="m", wakemodelused="", windresourcefile="", aeptotal=[], 
    # aepdirs=[], aepunits="MWh", baseyaml="initial-layouts/default_turbine_layout.yaml")

    # write turbine coordinates to a text file
    open(layout_directory * layout_filename, "w") do io
        write(io, "# turbine_x, turbine_y\n")
        writedlm(io, [turbine_x turbine_y])
    end

end

function plot_initial_layout_circle(figure_filename, layout_filename, rotor_diameter, boundary_center, boundary_radius; figure_directory="", layout_directory="", save_fig=true, show_fig=false)

    # ----- IMPORT TURBINE LOCATIONS -----
    # layout_filedata = ff.get_turb_loc_YAML(layout_directory * layout_filename)
    turbine_locations = readdlm(layout_directory * layout_filename, skipstart=1)
    turbine_x = turbine_locations[:,1]
    turbine_y = turbine_locations[:,2]

    # ----- CREATE PLOT -----
    figure()
    # plot turbine locations
    for i = 1:length(turbine_x)
        plt.gcf().gca().add_artist(plt.Circle((turbine_x[i],turbine_y[i]), rotor_diameter[1]/2.0, fill=false, color="C0", linestyle="--")) 
    end
    # add wind farm boundary to plot
    plt.gcf().gca().add_artist(plt.Circle((boundary_center[1],boundary_center[2]), boundary_radius, fill=false,color="C2"))
    # set figure axes
    axis("square")
    xlim(-boundary_radius-100,boundary_radius+100)
    ylim(-boundary_radius-100,boundary_radius+100)

    # ----- SAVE/SHOW FIGURE -----
    # save figure (if specified)
    if save_fig
        savefig(figure_directory * figure_filename)
    end
    # show figure (if specified)
    if show_fig
        plt.show()
    end

end


# ========= Everything below this point normally does not need to be edited. ==============


# ----- GENERATE LAYOUT FILES ----
if generate_layout_file
    # iterate through each layout number
    for layout_number = layout_numbers
        layout_filename = "initial-layout-" * lpad(layout_number,3,"0") * ".txt"
        generate_random_layout_circle_boundary(layout_filename, nturbines, rotor_diameter, boundary_center, boundary_radius; layout_directory=layout_directory, min_spacing=1.8, layout_number=layout_number)
    end
end

# ----- PLOT LAYOUTS -----
if generate_layout_figure
    # iterate through each layout number
    for layout_number = layout_numbers
        layout_filename = "initial-layout-" * lpad(layout_number,3,"0") * ".txt"
        figure_filename = "initial-layout-" * lpad(layout_number,3,"0") * ".png"
        plot_initial_layout_circle(figure_filename, layout_filename, rotor_diameter, boundary_center, boundary_radius; figure_directory=figure_directory, layout_directory=layout_directory)
    end
end