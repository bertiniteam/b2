import pybertini as pb
import numpy as np

x = pb.Variable('x')
y = pb.Variable('y')
z = pb.Variable('z')

vg = pb.VariableGroup([x,y,z])

sys = pb.System()

sys.add_function(x-y)
sys.add_function(x**2 + y**2 - 1)
sys.add_function(5*x**3 + 16*x*y**4 - 17*x*y*z - z**3)

sys.add_variable_group(vg)

C = pb.multiprec.Complex # unpack a type for convenience

gen = lambda : pb.random.complex_in_minus_one_to_one()
variable_values = np.array([gen(), gen(), gen()])

print(variable_values)

result = sys.eval(variable_values)


