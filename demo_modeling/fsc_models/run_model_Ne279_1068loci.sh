#!/bin/bash 
#SBATCH --job-name=fastsimcoal26
#SBATCH --partition=p_deenr_1,main,p_eoas_1
#SBATCH -N1
#SBATCH -n1                        ## 28 cores on amarel
#SBATCH --cpus-per-task=1         ## 28 cores on amarel 
############################
#SBATCH --mem=5000               ## 128GB /node of amarel
#SBATCH --time=24:00:00            ## max time is 3 days:  3-00:00:00  or  72:00:00
#SBATCH --export=ALL
#SBATCH --requeue       ## optional: this prevent the job from restarting if the node fails or the jobs is preempted.
#SBATCH -o /projects/f_mlp195/hoey/logs/slurm-%j.out   ## the default file name is "slurm-%j.out"
#SBATCH -e /projects/f_mlp195/hoey/logs/slurm-%j.err
#SBATCH --mail-type=NONE                # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=jahoey13@gmail.com # Email to which notifications will be sent

#############################################################
# standard output is saved in a file:  slurm-$SLURM_JOBID.out 
############################################################# 

mkdir -p /home/$USER/Ne279_1068loci/model6/slurm-out/$SLURM_JOBID
cp /home/$USER/fsc26_linux64/fsc26 /home/$USER/Ne279_1068loci/model6/model6.tpl /home/$USER/Ne279_1068loci/model6/model6.est /home/$USER/Ne279_1068loci/model6/model6_MSFS.obs /home/$USER/Ne279_1068loci/model6/slurm-out/$SLURM_JOBID
cd       /home/$USER/Ne279_1068loci/model6/slurm-out/$SLURM_JOBID

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

env >               /home/$USER/Ne279_1068loci/model6/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-all.out
env | grep SLURM >  /home/$USER/Ne279_1068loci/model6/slurm-out/$SLURM_JOBID/slurm-$SLURM_JOBID-env-SLURM.out


# start time
date

# start the simulation

srun ./fsc26 -t model6.tpl -e model6.est -n 100000 -m -u -M -L 40 -0 --numBatches 1 --cores 1 > /home/$USER/Ne279_1068loci/model6/slurm-out/$SLURM_JOBID/$SLURM_JOBID-fastsimcoal26.$SLURM_NNODES.$SLURM_NODELIST.$date.a.txt

#cp -R /mnt/scratch/$USER/$SLURM_JOBID/* /home/$USER/slurm-out/$SLURM_JOBID/

#end time
date
