# create structs to hold exclusion parameters
struct CircleExclusion{}
    center
    radius
end

struct PolygonExclusion{}
    vertices
    normals
end
PolygonExclusion(a) = PolygonExclusion(a, boundary_normals_calculator(a))



function get_winner_colome_exclusions()

    # houses
    houses_file_name = "input/exclusions/houses/houses_utm.txt"
    houses = readdlm(houses_file_name, skipstart=1)
    n_houses = length(houses[:,1])
    # ponds
    ponds_file_name = "input/exclusions/ponds/ponds_utm.txt"
    ponds = readdlm(ponds_file_name, skipstart=1)
    n_ponds = length(ponds[:,1])
    # misc
    misc_file_names = ["input/exclusions/misc/misc_utm.txt", "input/exclusions/misc/substation_utm.txt"]
    miscs = [readdlm(misc_file_names[i], skipstart=1) for i=1:length(misc_file_names)]
    n_miscs = length(miscs)
    # get number of circular exclusions
    n_circle_exclusions = n_houses + n_ponds + n_miscs
    # initialize circular exclusions vector
    turbine_circle_exclusions = Vector{CircleExclusion}(undef, n_circle_exclusions)
    # populate the vector
    k = 1
    for i = 1:n_houses
        house_x = houses[i,1]
        house_y = houses[i,2]
        turbine_circle_exclusions[k] = CircleExclusion([house_x, house_y], 370.0)
        k += 1
    end
    for i = 1:n_ponds
        ponds_x = ponds[i,1]
        ponds_y = ponds[i,2]
        ponds_r = ponds[i,3]
        turbine_circle_exclusions[k] = CircleExclusion([ponds_x, ponds_y], ponds_r)
        k += 1
    end
    for misc in miscs
        misc_x = misc[1]
        misc_y = misc[2]
        misc_r = misc[3]
        turbine_circle_exclusions[k] = CircleExclusion([misc_x, misc_y], misc_r)
        k += 1
    end

    # rivers
    rivers_file_names = ["input/exclusions/rivers/river_1_utm.txt"]
    rivers = [readdlm(rivers_file_names[i], skipstart=1) for i=1:length(rivers_file_names)]
    n_rivers = length(rivers)
    # roads
    roads_file_names = ["input/exclusions/roads/road_$(lpad(i,2,"0"))_utm.txt" for i=1:11]
    roads = [readdlm(roads_file_names[i], skipstart=1) for i=1:length(roads_file_names)]
    n_roads = length(roads)
    # get number of polygon exclusions
    n_polygon_exclusions = n_rivers + n_roads
    # initialize polygon exclusions vector
    turbine_polygon_exclusions = Vector{PolygonExclusion}(undef, n_polygon_exclusions)
    # populate the vector
    k = 1
    for river in rivers
        turbine_polygon_exclusions[k] = PolygonExclusion(river)
        k += 1
    end
    for road in roads
        turbine_polygon_exclusions[k] = PolygonExclusion(road)
        k += 1
    end

    return turbine_circle_exclusions, turbine_polygon_exclusions

end


