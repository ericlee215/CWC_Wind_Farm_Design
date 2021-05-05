# BYU WIND ENERGY CLUB   ADAM WELKER      04/2021
#
# GraphComps.jl -- a series of structs to assist in the creation
# of undirected graphs on a 2 euclidiean plane.


module GraphComps

    struct Edge

        src
        dest
        length

    end

    struct Node

        node_id
        loc :: Array{Float64}
        edges :: Vector

    end


    function node_addEdge(node :: Node, edge :: Edge)

        append!(node.edges, edge)

    end

    function getDistance(a :: Node, b :: Node)

        return sqrt((a.loc[1] - b.loc[1]) ^ 2 + (a.loc[2] - b.loc[2])^2)

    end



    function edge_getLength(edge :: Edge)

        edge.length = getDistance(edge.src, edge.dest)

    end


end

