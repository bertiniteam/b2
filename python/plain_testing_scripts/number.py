import pybertini


def exercise_ring(a, b):
	c0 = a+a
	c1 = a+b
	c2 = b+a

	d0 = a*a
	d1 = a*b
	d2 = b*a

	e0 = a-a
	e1 = a-b
	e2 = b-a


def exercize_field(a,b):
	exercize_ring(a,b)

	f0 = a/a
	f1 = a/b
	f2 = b/a




rings = [pybertini.numbers.int(2), pybertini.numbers.float(3),pybertini.numbers.rational(3,4),pybertini.numbers.complex(3,4)]

fields = [pybertini.numbers.float(3),pybertini.numbers.rational(3,4),pybertini.numbers.complex(3,4)]


for ii in rings:
	for jj in rings:
		exercise_ring(ii,jj)

for ii in rings:
	for jj in fields:
		exercise_field(ii,jj)

for ii in fields:
	for jj in fields:
		exercise_field(ii,jj)


