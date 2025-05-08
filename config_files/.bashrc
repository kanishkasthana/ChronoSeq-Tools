export http_proxy="mizar.ucsd.edu:4254"
export https_proxy=$http_proxy
export ftp_proxy=$http_proxy
export PATH=/opt/wangcluster/git/2.12.1/bin/:$PATH
export PATH=/home/kanishk/bin/:$PATH
export PATH=/home/kanishk/Drop-seq_tools-2.4.0/:$PATH

function startJupyter(){
    jupyter lab --no-browser --port=$1
}

free -h

function clusterInfo(){
    sinfo -N -l
}

function jobList(){
    squeue
}

function startBash(){
    srun -w $1 --pty bash
}

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/kanishk/anaconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/kanishk/anaconda3/etc/profile.d/conda.sh" ]; then
        . "/home/kanishk/anaconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/kanishk/anaconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

