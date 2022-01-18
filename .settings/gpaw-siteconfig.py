scalapack = True
fftw = True
elpa = True
libvdwxc = True

libraries = ['openblas', 'xc']
include_dirs = ['/usr/include', '/usr/local/include']
#library_dirs = ['/lib', '/lib/x86_64-linux-gnu', '/usr/lib', '/usr/lib/x86_64-linux-gnu', '/usr/local/lib']

if scalapack:
    libraries.append('scalapack-openmpi')

if fftw:
    libraries.append('fftw3')

if libvdwxc:
    libraries.append('vdwxc')

if elpa:
    import glob
    libraries.append('elpa')
    incs = glob.glob('/usr/include/elpa*')
    include_dirs.extend(incs)

extra_compile_args += ['-fopenmp']
extra_link_args += ['-fopenmp']

