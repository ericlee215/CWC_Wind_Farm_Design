"""

    boundary_normals_calculator(boundary_vertices)

Outputs the unit vectors perpendicular to each boundary of a shape, given the Cartesian coordinates for the shape's vertices.

# Arguments
- `boundary_vertices::Array{Float,1}` : n-by-2 array containing all the boundary vertices, counterclockwise

"""

function boundary_normals_calculator(boundary_vertices)

    # get number of vertices in shape
    nvertices = length(boundary_vertices[:, 1])

    # add the first vertex to the end of the array to form a closed loop
    boundary_vertices = [boundary_vertices; boundary_vertices[1,1] boundary_vertices[1,2]]

    # initialize array to hold boundary normals
    boundary_normals = zeros(nvertices, 2)

    # iterate over each boundary
    for i = 1:nvertices

        # create a vector normal to the boundary
        boundary_normals[i, :] = [-(boundary_vertices[i+1,2] - boundary_vertices[i,2]); boundary_vertices[i+1,1] - boundary_vertices[i,1]]
        
        # normalize the vector
        boundary_normals[i, :] = boundary_normals[i,:]/sqrt(sum(boundary_normals[i,:].^2))
    
    end

    return boundary_normals

end