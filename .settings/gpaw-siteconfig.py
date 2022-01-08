scalapack = True
fftw = True
elpa = True
libvdwxc = True

libraries = ['openblas', 'xc', 'fftw3', 'scalapack', 'vdwxc', 'elpa_openmp']
include_dirs = ['/usr/include', '/usr/local/include', '/usr/include/elpa_openmp-2021.05.002']
#library_dirs = ['/usr/lib', '/usr/local/lib']

extra_compile_args += ['-fopenmp']
extra_link_args += ['-fopenmp']

