import bertini as pb

x = pb.Variable('x')
y = pb.Variable('y')

sys = pb.System()

sys.add_function(x-y)
sys.add_function(x**2 + y**2 - 1)

sys.add_variable_group(pb.VariableGroup([x,y]))

# print(sys)

sys.homogenize()
sys.auto_patch()

# solver = pb.nag_algorithm.ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='rootsofunity')


import cProfile
import pstats


# Or profile a block of code
with cProfile.Profile() as pr:

    for n in range(10000):
        print(n)
        print(sys)
        solver = pb.nag_algorithm.ZeroDimSolver(sys, endgame='cauchy', mptype='adaptive', startsystem='rootsofunity')
        solver.solve()
    pr.dump_stats('profile_stats')
    # pr.sort_stats('cumulative')
    pr.print_stats()




print('\nsolutions in homogeneous space:')
for soln in solver.all_solutions():
    print(soln)

print('\nsolutions in original affine space:')
for soln in solver.all_solutions():
    print(sys.dehomogenize_point(soln))
