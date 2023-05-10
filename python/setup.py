from setuptools import find_packages, setup

EXCLUDE_FROM_PACKAGES = []

setup(name='pybertini',
      version='1.0.alpha3',
      description='Software for numerical algebraic geometry',
      url='http://github.com/bertiniteam/b2',
      author='Bertini Team',
      author_email='amethyst@uwec.edu',
      license='GPL3 with permitted additional clauses',
      packages=find_packages(exclude=EXCLUDE_FROM_PACKAGES),
      package_dir = {'pybertini': 'pybertini'},
      zip_safe=False)

# dependencies to add
# sphinxcontrib-bibtex
