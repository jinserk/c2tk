scalapack = True
fftw = True

libraries = ['openblas', 'xc', 'fftw3', 'scalapack-openmpi']
#include_dirs = ['/usr/include', '/usr/include/x86_64-linux-gnu']
#library_dirs = ['/usr/lib', '/usr/lib/x86_64-linux-gnu']

extra_compile_args += ['-fopenmp']
extra_link_args += ['-fopenmp']
