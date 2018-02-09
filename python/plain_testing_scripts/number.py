import pybertini


def exercise_ring(a, b):
	c1 = a+b
	c2 = b+a

	d1 = a*b
	d2 = b*a

	e1 = a-b
	e2 = b-a

	print('{} {} passed ring checks'.format(type(a),type(b)))

def exercise_field(a,b):
	exercise_ring(a,b)
	
	f1 = a/b
	f2 = b/a

	print('{} {} passed field checks'.format(type(a),type(b)))

def compare(a, b):
	a < b
	a <= b
	b < a
	b <= a

	a > b
	a >= b
	b > a
	b >= a

	print('{} {} passed compare checks'.format(type(a),type(b)))

def eq(a,b):
	a == b
	a != b
	b == a
	b != a
	a == a
	a != a
	b == b
	b != b

	print('{} {} passed eq checks'.format(type(a),type(b)))

rings = [pybertini.numbers.int(2)]

fields = [pybertini.numbers.float(3),pybertini.numbers.rational(3,4),pybertini.numbers.complex(5,6)]

ordereds = [pybertini.numbers.int(2), pybertini.numbers.float(3),pybertini.numbers.rational(3,4)]

all_types = [pybertini.numbers.int(2), pybertini.numbers.float(3),pybertini.numbers.rational(3,4), pybertini.numbers.complex(5,6)]


for ii in rings:
	for jj in rings:
		exercise_ring(ii,jj)

for ii in fields:
	for jj in fields:
		exercise_field(ii,jj)

# then mixed products
for ii in rings:
	for jj in fields:
		exercise_field(ii,jj)


for ii in ordereds:
	compare(ii, 1)

for ii in all_types:
	eq(ii,ii)

