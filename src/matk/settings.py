import os
from pathlib import Path

from dotenv import load_dotenv, dotenv_values


# code default
# all values defined here should be str, bytes or os.PathLike objs
default_env = dict(
    SCRATCH_PATH="/matk/scratch",
    NPROC=f"{os.cpu_count()}",
)


# import env
env = {
    # default value
    **default_env,
    # load .env if exists
    **dotenv_values('.env'),
    # or use env vars
    **os.environ,
}


# actual configs
SCRATCH_PATH = env["SCRATCH_PATH"]
Path(SCRATCH_PATH).mkdir(mode=0o755, parents=True, exist_ok=True)

NPROC = env["NPROC"]
