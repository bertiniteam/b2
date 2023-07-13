import pybertini as pb
import numpy as np

pb.default_precision(500)

x = pb.Variable('x')
y = pb.Variable('y')
z = pb.Variable('z')

vg = pb.VariableGroup([x,y,z])

sys = pb.System()

sys.add_function(x-y)
sys.add_function(x**2 + y**2 - 1)
sys.add_function(5*x**3 + 16*x*y**4 - 17*x*y*z - z**3)

sys.add_variable_group(vg)

random_complex = lambda : pb.random.complex_in_minus_one_to_one()

n_iterations = 10000
print(f'evaluating system at random complex point {n_iterations} times at precision {pb.default_precision()}')
for ii in range(n_iterations):
	variable_values = np.array([random_complex(), random_complex(), random_complex()])
	result = sys.eval(variable_values)


