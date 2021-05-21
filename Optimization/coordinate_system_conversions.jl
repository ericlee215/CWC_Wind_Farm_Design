using Geodesy
using DelimitedFiles


function UTMfromLLAfile(lla_file_name, utm_file_name, utm_zone, northern_hemisphere=true)

    # import boundary coordinates 
    boundary_vertices = readdlm(lla_file_name, skipstart=1)

    # store coordinates into LLA objects
    n_vertices = length(boundary_vertices[:,1])
    lla_coordinates = Vector{LLA{Float64}}(undef, n_vertices)
    for i = 1:n_vertices
        lla_coordinates[i] = LLA(boundary_vertices[i,1], boundary_vertices[i,2], 0.0)
    end

    # convert to UTM
    utm_from_lla = UTMfromLLA(utm_zone, northern_hemisphere, wgs84)
    utm_coordinates = map(utm_from_lla, lla_coordinates)

    # store UTM coordinates into an n by 2 array
    boundary_vertices_utm = zeros(n_vertices, 2)
    for i = 1:n_vertices
        boundary_vertices_utm[i,:] = [utm_coordinates[i].x, utm_coordinates[i].y]
    end

    # write boundary coordinates to a file
    open(utm_file_name, "w") do io
        write(io, "# boundary coordinates (UTM zone $utm_zone, north=$northern_hemisphere) (x,y)\n")
        writedlm(io, boundary_vertices_utm)
    end

end

# for i = 1:11
#     # import boundary coordinates 
#     lla_file_name = "input/exclusions/roads/road_$(lpad(i,2,"0")).txt"
#     utm_file_name = "input/exclusions/roads/road_$(lpad(i,2,"0"))_utm.txt"

#     UTMfromLLAfile(lla_file_name, utm_file_name, 14, true)
# end

lla_file_name = "input/exclusions/substation/substation.txt"
utm_file_name = "input/exclusions/substation/substation_utm.txt"

UTMfromLLAfile(lla_file_name, utm_file_name, 14, true)
