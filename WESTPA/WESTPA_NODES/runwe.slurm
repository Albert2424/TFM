#!/bin/bash
#SBATCH -J WESTPA
#SBATCH -e WEST.err
#SBATCH -o WEST.out
#SBATCH --gres=gpu:2
#SBATCH --cpus-per-task=32
#SBATCH -N 2
#SBATCH --ntasks-per-node=2
#SBATCH --time=0-01:00:00
#SBATCH --mem=48G


#set -x
cd $SLURM_SUBMIT_DIR
source env.sh || exit 1
#env | sort

rm -f *.log
cd common_files
python analyse.py #generate pickle of the proteins.
cd ..


date

if [ ! -d "traj_segs" ]; then
    ./init.sh
fi

cd $WEST_SIM_ROOT
SERVER_INFO=$WEST_SIM_ROOT/west_zmq_info-$SLURM_JOBID.json

# start server
w_run --work-manager=zmq --n-workers=0 --zmq-mode=master --zmq-write-host-info=$SERVER_INFO --zmq-comm-mode=tcp &> west.log &

# wait on host info file up to one minute
for ((n=0; n<60; n++)); do
    if [ -e $SERVER_INFO ] ; then
        echo "== server info file $SERVER_INFO =="
        cat $SERVER_INFO
        break
    fi
    sleep 1
done

echo 'Number of GPUs:' $SLURM_GPUS_ON_NODE

# exit if host info file doesn't appear in one minute
if ! [ -e $SERVER_INFO ] ; then
    echo 'server failed to start'
    exit 1
fi

# start clients, with the proper number of cores on each

scontrol show hostname $SLURM_NODELIST >& SLURM_NODELIST.log
node_id=-1
for node in $(scontrol show hostname $SLURM_NODELIST); do
    echo 'lalalala'
    ((node_id++))
    srun -lN 2 -r $node_id -n $SLURM_GPUS_ON_NODE $PWD/westpa_scripts/node.sh $SLURM_SUBMIT_DIR $SLURM_JOBID $node $CUDA_VISIBLE_DEVICES --work-manager=zmq --n-workers=$SLURM_GPUS_ON_NODE --zmq-mode=client --zmq-read-host-info=$SERVER_INFO --zmq-comm-mode=tcp & #MODIFY --n-workers to the same number of gpus you have!
    echo 'after lalalala'
done


wait

date
