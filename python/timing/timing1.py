import pybertini as pb




import timeit

result = timeit.timeit('solver.solve()', number=10000, globals=globals())

print(result)
	


