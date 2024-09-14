from setuptools import find_packages, setup

with open('README.md', 'r') as f:
    readme = f.read()

# ext_modules = [
#     Extension(
#         'cpp_bindings',
#         sources=['src/num_solvers/cpp_bindings/binder/bindings.cpp'],
#         include_dirs=[pybind11.get_include(), "src/num_solvers/cpp_bindings/include"],
#         language='c++',
#         extra_compile_args=['-std=c++11'],
#     ),
# ]

setup(name='num_solvers',
      version='0.1.0',
      packages=find_packages(where=["src", "src/num_solvers/cpp_bindings/binder"]),
      url='https://github.com/syedalimohsinbukhari/Solvers',
      license='MIT',
      author='Syed Ali Mohsin Bukhari, Astrophysics and Python',
      author_email='syedali.b@outlook.com, astrophysicsandpython@gmail.com',
      description='A library for some basic numerical num_solvers.',
      long_description=readme,
      long_description_content_type="text/markdown",
      python_requires=">=3.9",
      install_requires=["numpy==1.26.4", "setuptools", "custom-inherit", "umatrix"],
      include_package_data=True,
      # ext_modules=ext_modules,
      package_data={
          'num_solvers': ['cpp_bindings/binder/*.cpp', 'cpp_bindings/include/*']
      },
      zip_safe=False,
      classifiers=[
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3.9",
          "Programming Language :: Python :: 3.10",
          "Programming Language :: Python :: 3.11"],
      )
