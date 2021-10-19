import ase.calculators.orca as aco

from .calculator import C2TKFileIOCalculator
from .. import settings


class ORCA(C2TKFileIOCalculator, aco.ORCA):

    def __init__(self, *args, **kwargs):
        super().__init__(mpi_embed_cmd=True, *args, **kwargs)

        self.command = f'{settings.ORCA_PATH}/orca PREFIX.inp >> PREFIX.out 2> PREFIX.err'

