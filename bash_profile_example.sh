## THIS IS AN EXAMPLE. YOU NEED TO SAVE THIS AS .bash_profile IN YOUR HOME DIR##
# .bash_profile

# Get the aliases and functions
if [ -f ~/.bashrc ]; then
        . ~/.bashrc
fi

# User specific environment and startup programs

PATH=$PATH:$HOME/.local/bin:$HOME/bin

export PATH

conda activate uroconda381
export http_proxy=http://172.28.7.1:3128
export https_proxy=http://172.28.7.1:3128
export all_proxy=http://172.28.7.1:3128
export no_proxy=localhost,*.hpc.mssm.edu,*.chimera.hpc.mssm.edu,172.28.0.0/16

module load python
module load gcc
module load singularity