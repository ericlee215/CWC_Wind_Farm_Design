#################################################################
# Simple example of pseudo object-oriented programming in Julia
#################################################################

struct ExampleObject
    value_1
    value_2
    string_1
    string_2
end

function sum_values(example_object)
    return example_object.value_1 + example_object.value_2
end

function combine_strings(example_object)
    return example_object.string_1 * example_object.string_2
end

example_object_1 = ExampleObject(1., 2., "hello, ", "world")

println(sum_values(example_object_1))
println(combine_strings(example_object_1))


#################################################################
# pop and append
#################################################################

array_1 = [1,2,3,4,5]
println("\nOriginal array: ", array_1)
pop!(array_1)   # exclamation mark is Julia's naming convention for a function that mutates its inputs
append!(array_1, 10)
println("Modified array: ", array_1)
println()
