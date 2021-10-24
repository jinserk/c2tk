import os
from pathlib import Path

from dotenv import load_dotenv, dotenv_values


c2tk_path = Path("~/c2tk").expanduser()
scratch_path = c2tk_path.joinpath("scratch")

# code default
# all values defined here should be str, bytes or os.PathLike objs
default_env = dict(
    C2TK_PATH=str(c2tk_path),
    SCRATCH_PATH=str(scratch_path),
    NPROC=f"{os.cpu_count()}",
    ORCA_PATH="/opt/orca501",
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
C2TK_PATH = env["C2TK_PATH"]
SCRATCH_PATH = env["SCRATCH_PATH"]
Path(SCRATCH_PATH).mkdir(mode=0o755, parents=True, exist_ok=True)

NPROC = env["NPROC"]

ORCA_PATH = env["ORCA_PATH"]
