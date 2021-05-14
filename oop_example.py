#################################################################
# Simple example of object-oriented programming in Python
#################################################################

class ExampleObject:

    def __init__(self, value_1, value_2, string_1, string_2):
        self.value_1 = value_1
        self.value_2 = value_2
        self.string_1 = string_1
        self.string_2 = string_2

    def sum_values(self):
        return self.value_1 + self.value_2
        
    def combine_strings(self):
        return self.string_1 + self.string_2

example_object_1 = ExampleObject(1., 2., "hello, ", "world")

print(example_object_1.sum_values())
print(example_object_1.combine_strings())


#################################################################
# pop and append
#################################################################

list_1 = [1,2,3,4,5]
print("\nOriginal list: ", list_1)
list_1.pop()
list_1.append(10)
print("Modified list: ", list_1)
print()