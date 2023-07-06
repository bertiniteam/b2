import pybertini as pb

x, y, z = pb.Variable('x'), pb.Variable('y'), pb.Variable('z')

f1 = x**2 + y**2 - 1
f2 = x+y

sys = pb.System()
sys.add_function(f1)
sys.add_function(f2)

sys.add_variable_group(pb.VariableGroup([x,y]))

solver = pb.nag_algorithms.ZeroDimCauchyAdaptivePrecisionTotalDegree(sys)


import timeit

result = timeit.timeit('solver.solve()', number=10000, globals=globals())

print(result)
	


