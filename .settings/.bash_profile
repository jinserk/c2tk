. /etc/profile
module load mpi/openmpi-x86_64

export PYENV_ROOT="$HOME/.pyenv"
export PATH="$PYENV_ROOT/bin:$PATH"
eval "$( pyenv init --path )"
