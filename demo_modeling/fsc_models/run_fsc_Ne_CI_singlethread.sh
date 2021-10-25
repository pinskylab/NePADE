#!/bin/bash 
#SBATCH --job-name=fastsimcoal26
#SBATCH --partition=p_eoas_1,p_deenr_1,main
#SBATCH -N1
#SBATCH -n1                        ## 28 cores on amarel
#SBATCH --cpus-per-task=1         ## 28 cores on amarel 
############################
#SBATCH --mem=5000               ## 128GB /node of amarel
#SBATCH --time=24:00:00            ## max time is 3 days:  3-00:00:00  or  72:00:00
#SBATCH --export=ALL
#SBATCH --array=1-100%50
#SBATCH --requeue       ## optional: this prevent the job from restarting if the node fails or the jobs is preempted.
#SBATCH -o /scratch/jah414/logs/slurm-%A_%a.out   ## %A is the jobid and %a is the array index
#SBATCH -e /scratch/jah414/logs/slurm-%A_%a.err   ## /dev/null
#SBATCH --mail-type=END                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jahoey13@gmail.com # Email to which notifications will be sent

#############################################################
# standard output is saved in a file:  slurm-$SLURM_JOBID.out 
############################################################# 

echo "Starting task $SLURM_ARRAY_TASK_ID"
DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" dir_names.txt)
cd $DIR

mkdir -p /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID
cp /home/$USER/fsc26_linux64/fsc26 /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID
cp /home/$USER/Ne279_1068loci/model6/model6.tpl /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID/Ne279_1068loci.tpl
cp /home/$USER/Ne279_1068loci/model6/model6.est /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID/Ne279_1068loci.est
cp /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/Ne279_1068loci_MSFS.obs /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID
cp /home/$USER/Ne279_1068loci/model6/model6.pv /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID/Ne279_1068loci.pv
cd /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID

#Create a personal dirctory on the node scratch disk if needed.
#mkdir -p /mnt/scratch/$USER/$SLURM_JOBID


# this will tell you when the job started and the host. 
date=`date "+%Y.%m.%d-%H.%M.%S"`
hostname=`hostname`

# to print the variable -- echo
echo $date
echo $hostname

# obtain the environment variables 
# this is useful for reference and troubleshooting issues.

env >               /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-all.out
env | grep SLURM >  /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-SLURM.out


# start time
date

# start the simulation

srun ./fsc26 -t Ne279_1068loci.tpl -e Ne279_1068loci.est -n 100000 -m -u -M -L 40 -0 --numBatches 1 --cores 1 --initValues Ne279_1068loci.pv  > /projects/f_mlp195/hoey/bootSFS_279fish_1068loci_model6/$DIR/slurm-out/$SLURM_JOBID/$SLURM_JOBID-fastsimcoal26.$SLURM_NNODES.$SLURM_NODELIST.$date.a.txt

#cp -R /mnt/scratch/$USER/$SLURM_JOBID/* /home/$USER/slurm-out/$SLURM_JOBID/

#end time
date
