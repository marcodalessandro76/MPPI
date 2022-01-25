"""
This module manages the parameters used to define the mpi and omp parellelization strategy.
"""
import os

def build_slurm_header(pars):
    """
    Define the header of the slurm script. Note that the name variable is not present in the
    RunRules parameter but it is added by the Calculator in the run_options dictionary.

    Args:
        pars (:py:class:`dict`) : dictionary with the structure of an instance of
            the :class:`RunRules`

    Return:
        :py:class:`list` : a list with the lines of the header of the script

    """
    name = pars.get('name','default')
    job = 'job_'+name

    lines = []
    lines.append('#!/bin/bash')
    if pars['account'] is not None:
        lines.append('#SBATCH --account %s'%pars['account'])
    lines.append('#SBATCH --nodes=%s              ### Number of nodes'%pars['nodes'])
    lines.append('#SBATCH --ntasks-per-node=%s    ### Number of MPI tasks per node'%pars['ntasks_per_node'])
    lines.append('#SBATCH --cpus-per-task=%s      ### Number of HT per task'%pars['cpus_per_task'])
    lines.append('#SBATCH --mem %s             ### Memory per node'%pars['memory'])
    if pars['partition'] is not None:
        lines.append('SBATCH --partition %s'%pars['partition'])
    if pars['time'] is not None:
        lines.append('#SBATCH --time %s         ### Walltime, format: HH:MM:SS'%pars['time'])
    lines.append('#SBATCH --job-name=%s'%job)
    lines.append('#SBATCH --output=%s.out'%job)
    lines.append('')
    lines.append('export OMP_NUM_THREADS=%s'%pars['omp_num_threads'])
    lines.append('')
    lines.append('echo "Cluster name $SLURM_CLUSTER_NAME"')
    lines.append('echo "Job name $SLURM_JOB_NAME "')
    lines.append('echo "Job id $SLURM_JOB_ID"')
    lines.append('echo "Job nodelist $SLURM_JOB_NODELIST"')
    lines.append('echo "Number of nodes $SLURM_JOB_NUM_NODES"')
    lines.append('echo "Number of tasks $SLURM_NTASKS"')
    lines.append('echo "Number of tasks per node $SLURM_TASKS_PER_NODE"')
    lines.append('echo "Number of threads per task $SLURM_CPUS_PER_TASK"')
    lines.append('echo "OMP_NUM_THREADS=$OMP_NUM_THREADS"')
    lines.append('')
    lines.append('echo "###############End of the header section###############"')
    lines.append('')

    return lines

def mpi_command(pars):
    """
    Define the mpi run command.

    Args:
        pars (:py:class:`dict`) : dictionary with the structure of an instance of
            the :class:`RunRules`

    Return:
        :py:class:`string` : command that defines the mpi run

    """
    mpi_run = 'mpirun'
    if pars['scheduler'] == 'direct':
        mpi_run += ' ' + '-np %s'%pars['rules']['mpi']
    if pars['scheduler'] == 'slurm':
        ntasks = pars['ntasks_per_node']*pars['nodes']
        mpi_run += ' ' + '-np %s'%ntasks
        if pars['map_by'] is not None:
            mpi_run += ' ' + '--map-by %s:PE=%s'%(pars['map_by'],pars['pe'])
        if pars['rank_by'] is not None:
            mpi_run += ' ' + '--rank-by %s'%pars['rank_by']

    return mpi_run

class RunRules(dict):
    """
    Defines and manage a dictionary with all the parameters needed to set up the rules
    for parallel computing.

    Parameters:
        scheduler (:py:class:`string`) : choose the scheduler used to submit the job
        omp_num_threads (:py:class:`int`) : the value of the environment variable OMP_NUM_THREADS
        mpi (:py:class:`int`) : number of mpi processes. This parameter is used only if `scheduler` is
            ``direct``, otherwise the `nodes` and `ntasks_per_node` parameters as used.
        nodes (:py:class:`int`) : slurm nodes variable
        ntasks_per_node (:py:class:`int`) : slurm ntasks-per-node variable
        cpus_per_task (:py:class:`int`) : slurm cpus-per-task variable
        partition (:py:class:`string`) : slurm parition variable
        memory (:py:class:`string`) : slurm mem variable
        account (:py:class:`string`) : slurm output variable

    """

    def __init__(self,scheduler='direct',omp_num_threads=os.environ.get('OMP_NUM_THREADS', 1),mpi=4,
                nodes=1,ntasks_per_node=1,cpus_per_task=1,partition=None,memory='124GB',
                account=None,time=None,map_by=None,pe=1,rank_by=None):
        if scheduler == 'direct':
            rules = dict(mpi=mpi,omp_num_threads=omp_num_threads)
            dict.__init__(self,scheduler=scheduler,**rules)
        if scheduler == 'slurm':
            rules=dict(nodes=nodes,ntasks_per_node=ntasks_per_node,cpus_per_task=cpus_per_task,
            omp_num_threads=omp_num_threads,partition=partition,memory=memory,
            account=account,time=time,map_by=map_by,pe=pe,rank_by=rank_by)
            dict.__init__(self,scheduler=scheduler,**rules)
