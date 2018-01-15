give a polynomial system and solve it

::

x = pybertini.function_tree.symbol.Variable("x") #yes, you can make a variable not match its name...
y = pybertini.function_tree.symbol.Variable("y")
f = x**2 + y**2 -1
g = x+y

sys = pybertini.System()
sys.add_function(f)
sys.add_function(g)

grp = _pybertini.VariableGroup()
grp.append(x)
grp.append(y)
sys.add_variable_group(grp)

td = pybertini.start_system.TotalDegree(sys)

zd = pybertini.algorithms.ZeroDim(sys, td)



zd.num_successful_paths()
zd.solutions[5]
zd.num_nonsingular_solutions()
zd.nonsingular_solutions()

zd.real_solutions[4]
zd.set_real_tolerance()
# must use new tolerance to reclassify solutions
zd.real_solutions[4]



#######
regen extension

solve a system.

solve it

