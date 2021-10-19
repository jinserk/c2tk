import os
import subprocess as sp
import shlex
from collections.abc import Iterable

from ase.calculators.calculator import (
    Calculator, FileIOCalculator,
    all_changes, CalculationFailed
)

from .. import settings, is_mpi_enabled


class C2TKFileIOCalculator(FileIOCalculator):

    def __init__(self, mpi_embed_cmd=False, *args, **kwargs):
        if 'directory' not in kwargs:
            kwargs['directory'] = settings.SCRATCH_PATH
        super().__init__(*args, **kwargs)
        self.mpi_embed_cmd = mpi_embed_cmd

    def calculate(self, atoms=None, properties=['energy'],
                  system_changes=all_changes):
        Calculator.calculate(self, atoms, properties, system_changes)
        self.write_input(self.atoms, properties, system_changes)
        self.execute()
        self.read_results()

    def execute(self):
        if self.command is None or not isinstance(self.command, str):
            raise CalculatorSetupError(
                'Please set ${} environment variable '
                .format('ASE_' + self.name.upper() + '_COMMAND') +
                'or supply the command keyword')
        command = self.command
        command = command.replace('PREFIX', self.prefix)

        if is_mpi_enabled() and not self.mpi_embed_cmd:
            command = f'mpiexec -np {settings.NPROC} {command}'

        #real_command = shlex.split(command)
        #print(real_command)

        try:
            proc = sp.Popen(command, cwd=self.directory, shell=True, env=settings.env)
        except OSError as err:
            # Actually this may never happen with shell=True, since
            # probably the shell launches successfully.  But we soon want
            # to allow calling the subprocess directly, and then this
            # distinction (failed to launch vs failed to run) is useful.
            msg = 'Failed to execute "{}"'.format(command)
            raise EnvironmentError(msg) from err

        errorcode = proc.wait()

        if errorcode:
            path = os.path.abspath(self.directory)
            msg = ('Calculator "{}" failed with command "{}" failed in '
                   '{} with error code {}'.format(self.name, command,
                                                  path, errorcode))
            raise CalculationFailed(msg)
