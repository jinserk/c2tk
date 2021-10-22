import os
import subprocess as sp

import matplotlib


try:
    # load environment-modules env
    exec(open('/usr/share/Modules/init/python.py').read())
    module('load', 'mpi/openmpi-x86_64')
except:
    pass


def is_mpi_enabled():
    from shutil import which
    return which('mpirun') is not None

matplotlib.use('TKAgg')
