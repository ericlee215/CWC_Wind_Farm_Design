# BYU WIND ENERGY CLUB        ADAM WELKER      04/2021
#
# clb_analysis.jl -- a series of functions that analyzes the
# most efficient orientation of cabling for a wind farm
# 
#
# ASSUMPTIONS -- cabling has a linear cost per foot
#             -- cabling is done in parrellel, meaning that back-tracking does not cost money
#             -- Wind farm exists on a flat, homogenous space (no rocks/ obstacles to get around)
#             -- Wind turbines can be hooked up from any angle
#
# APPROACH -- The approach will be to make a minimum spanning tree between the turbines using Prim's Algorithm
# and then use the tree cost multiplied by the linear cost per m or ft of the wiring

module Cable_Analysis

    using DataStructures
    using PyPlot
    include("GraphComps.jl")


    # Given a newtwork of nodes in Graph(V,E)
    # this function will preform prim's algorithm
    # on the graph and will return the set of edges E
    # that make up the minimum spanning tree
    function primAnalysis(network; return_network=false, plot_tree=false)

        # make the cost vector
        cost = []

        for i in 1:(length(network))

            push!(cost, Inf)

        end 

        cost[1] = 0

        # make the prev vector
        prev = []

        for i in 1:length(network)

            push!(prev, -1)

        end

        #make the priority queue and insert all values
        queue = PriorityQueue()

        for i in 1:length(network)

            queue[network[i]] = cost[i]

        end

        # begin the main while loop of prim's Algorithm
        while length(queue) != 0

            turbine = dequeue!(queue)

            # for each edge
            for edge in turbine.edges

                #if the edge cost is less from this node
                if cost[edge.dest.node_id] > edge.length

                    # and the edge is not being used
                    if prev[edge.src.node_id] != edge.dest.node_id

                        # check for loops
                        previous = deepcopy(prev)
                        previous[edge.dest.node_id] = edge.src.node_id

                        if !loopCheck(previous, edge.dest.node_id)

                            # change the arrangement
                            prev[edge.dest.node_id] = edge.src.node_id
                            cost[edge.dest.node_id] = edge.length
                            queue[edge.dest] = edge.length

                        end
  
                    end

                end

            end
            
        end

        #plot the tree

        if plot_tree

            for i in network

                PyPlot.plot(i.loc[1],i.loc[2],"ro")

            end

            #link each node on the plotted tree
            for i in 2:length(prev)
                
                x = [network[i].loc[1], network[prev[i]].loc[1]]
                y = [network[i].loc[2], network[prev[i]].loc[2]]
                
                PyPlot.plot(x,y,"b-")

            end

            PyPlot.show()

        end


        if return_network

            # get cable layout
            n_cables= length(network) - 1
            cables = [[zeros(2),zeros(2)] for i=1:n_cables]
            nodes = [zeros(Int, 2) for i=1:n_cables]

            for i in 2:length(prev)

                cables[i-1][1][:] = network[i].loc
                cables[i-1][2][:] = network[prev[i]].loc

                nodes[i-1][1] = network[i].node_id
                nodes[i-1][2] = network[prev[i]].node_id
                
            end

        end


        # sum up the cost array
        total_cost = 0

        for i in cost

            total_cost += i

        end

        if return_network

            return total_cost, cables, nodes

        else

            return total_cost

        end

    end


    # Given a vectors of x and y coordinates, this
    # function will turn the sets in to a graph and 
    # return the cost of the minimum spanning tree multiplied
    # by a constant cost variable
    function cbl_analysis(x_coordinates, y_coordinates, linear_cost; return_network=false, plot_tree=false)

        cost = 0

        windFarm = []

        #make each turbine and connect it to the others
        for i in 1:length(x_coordinates)
          
            node_id = i
            loc = [x_coordinates[i],y_coordinates[i]]
            turbine = GraphComps.Node(node_id,loc,[])

            push!(windFarm, turbine)

        end

        # Connect turbines to each other
        for i in (1:length(windFarm))

            for j in (1:length(windFarm))

                if i != j
                    
                    edge = GraphComps.Edge(windFarm[i],windFarm[j],GraphComps.getDistance(windFarm[i],windFarm[j]))

                    push!(windFarm[i].edges,edge)

                end

            end

        end

        # find the distance cost through prim's algorithm and then
        # multiply it by the linear cost    
        if return_network
            cost, cables, nodes = primAnalysis(windFarm, return_network=true, plot_tree=plot_tree)
        else
            cost = primAnalysis(windFarm, plot_tree=plot_tree)
        end


        # #############################
        # #plot the tree

        # if plot_tree

        #     for i in windFarm

        #         PyPlot.plot(i.loc[1],i.loc[2],"ro")

        #     end

        #     #link each node on the plotted tree
        #     for i in 2:length(prev)
                
        #         x = [network[i].loc[1], network[prev[i]].loc[1]]
        #         y = [network[i].loc[2], network[prev[i]].loc[2]]
                
        #         PyPlot.plot(x,y,"b-")

        #     end


        #     PyPlot.show()

        # end
        # ###############################

        if return_network

            return cost * linear_cost, cables, nodes

        else

            return cost * linear_cost

        end

    end


    #A function to check for loops in a previous index array
    function loopCheck(prev, index)

        found = falses(length(prev))

        indexExists = prev[index] != -1

        while  indexExists && !found[index]

            found[index] = true
            index = prev[index]
            indexExists = prev[index] != -1

        end

        if found[index]

            return true

        else

            return false

        end

    end
end