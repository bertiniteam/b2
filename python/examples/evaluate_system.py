import pybertini as pb
import numpy as np

x = pb.Variable('x')
y = pb.Variable('y')

vg = pb.VariableGroup([x,y])

sys = pb.System()

sys.add_function(x-y)
sys.add_function(x**2 + y**2 - 1)

sys.add_variable_group(vg)

C = pb.multiprec.Complex

variable_values = np.array([C(0), C(0)])

result = sys.eval(variable_values)

sys.homogenize()
sys.auto_patch()

variable_values = np.array([C(1), C(1), C(1)])

sys.rescale_point_to_fit_patch_in_place(variable_values) # failing due to non-const reference


solver = pb.nag_algorithm.ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)

solver.solve()

for soln in solver.solutions():

	print(sys.dehomogenize_point(soln))