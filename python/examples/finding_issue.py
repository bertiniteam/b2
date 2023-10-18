# finding issue

import pybertini as pb

x = pb.Variable('x')
y = pb.Variable('y')

f = pb.function_tree.root.Function(x)
f.name = 'f'
x.set_current_value(1)

print(f.eval_d())

vg = pb.VariableGroup([x,y])

sys = pb.System()
sys.add_function(f)

print(sys)
# sys.add_function(x**2 + y**2 - 1)

sys.add_variable_group(vg)
