import pybertini as pb
import numpy as np

x = pb.Variable('x')
y = pb.Variable('y')
z = pb.Variable('z')

vg = pb.VariableGroup([x,y,z])

sys = pb.System()

sys.add_function(x-y)
sys.add_function(x**2 + y**2 - 1)
sys.add_function(5*x**3 + 16*x*y**4 - 17x*y*z - z**3)

sys.add_variable_group(vg)

C = pb.multiprec.Complex # unpack a type for convenience

variable_values = np.array([C(0), C(0), C(0)])

result = sys.eval(variable_values)


