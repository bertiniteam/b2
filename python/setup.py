from setuptools import find_packages, setup

EXCLUDE_FROM_PACKAGES = []

setup(name='pybertini',
      version='1.0.alpha4',
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


# from setuptools.command.egg_info import egg_info

# class EggInfoCommand(egg_info):

#     def run(self):
#         if "build" in self.distribution.command_obj:
#             build_command = self.distribution.command_obj["build"]

#             self.egg_base = build_command.build_base

#             self.egg_info = os.path.join(self.egg_base, os.path.basename(self.egg_info))

#         egg_info.run(self)

# setup(
#     # ...
#     cmdclass={
#         "egg_info": EggInfoCommand,
#     },
#     #...
# )
