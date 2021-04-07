#=
Functions to test simple optimization algorithms.
=#

"""
    quadratic_ndim(x; shift=0.0)

Multidimensional quadratic function. Global minimum occurs when inputs are all zeros.

# Arguments
- `x::Array{Float64,1}`: Input to the quadratic function.
- `shift::Float64`: Amount by which to shift the global minimum in each direction of x. Default is no shift.
"""
function quadratic_ndim(x; shift=0.0)

    n = length(x)
    shift_vec = zeros(length(x)) .+ shift

    return sum((x .- shift_vec).^2)
end


"""
    rosenbrock_ndim(x)

Multidimensional Rosenbrock function. Commonly used to test optimization algorithms. (https://en.wikipedia.org/wiki/Rosenbrock_function#Multidimensional_generalisations)
- The 2-dimensional Rosenbrock function has a global minimum at [1,1]
- The 3-dimensional Rosenbrock function has a global minimum at [1,1,1]
- The 4-dimensional and higher Rosenbrock function has a global minimum at [1,1,1,1,...] and a local minimum near [-1,1,1,1,...]

# Arguments
- `x::Array{Float64,1}`: Input to the Rosenbrock function. Should have at least two elements.
"""
function rosenbrock_ndim(x)

    n = length(x)
    f = 0.0
    for i = 1:n-1
        f += 100*(x[i+1] - x[i]^2)^2 + (1 - x[i])^2
    end

    return f
end
