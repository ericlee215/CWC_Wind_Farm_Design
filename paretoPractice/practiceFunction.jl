#  Practice for making a epsilon method pareto front

using Ipopt


function obj1(x)

    x^2 - 2*x - 25
end


function obj2(x)

    x^4 - 4*x^3 + 50
end


epsilon = [200, 150, 100, 50, 25 , 10 , 5] # an epsilon value for obj2
results = []


for e in epsilon

    for x in domain

        

    end

end