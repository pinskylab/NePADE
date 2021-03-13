#!/bin/bash 
#SBATCH --job-name=fastsimcoal26
#SBATCH --partition=main
#SBATCH -N1
#SBATCH -n1                        ## 28 cores on amarel
#SBATCH --cpus-per-task=1         ## 28 cores on amarel 
############################
#SBATCH --mem=2000               ## 128GB /node of amarel
#SBATCH --time=30:00            ## max time is 3 days:  3-00:00:00  or  72:00:00
#SBATCH --export=ALL
#SBATCH --requeue       ## optional: this prevent the job from restarting if the node fails or the jobs is preempted.
#SBATCH -o /projects/f_mlp195/hoey/logs/slurm-%j.out   ## the default file name is "slurm-%j.out"
#SBATCH -e /projects/f_mlp195/hoey/logs/slurm-%j.err
#SBATCH --mail-type=END                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jahoey13@gmail.com # Email to which notifications will be sent

#############################################################
# standard output is saved in a file:  slurm-$SLURM_JOBID.out 
############################################################# 

mkdir -p /home/$USER/boot_model6/slurm-out/$SLURM_JOBID
cp /home/$USER/fsc26_linux64/fsc26 /home/$USER/boot_model6/model6_maxL.par /home/$USER/boot_model6/slurm-out/$SLURM_JOBID
cd       /home/$USER/boot_model6/slurm-out/$SLURM_JOBID

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

env >               /home/$USER/boot_model6/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-all.out
env | grep SLURM >  /home/$USER/boot_model6/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-SLURM.out


# start time
date

# start the simulation

srun ./fsc26 -i model6_maxL.par -n 100 -m -j -s1196 -x -I -0 > /home/$USER/boot_model6/slurm-out/$SLURM_JOBID/$SLURM_JOBID-fastsimcoal26.$SLURM_NNODES.$SLURM_NODELIST.$date.a.txt

#cp -R /mnt/scratch/$USER/$SLURM_JOBID/* /home/$USER/slurm-out/$SLURM_JOBID/

#end time
date
