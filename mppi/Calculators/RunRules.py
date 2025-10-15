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
    lines.append('#SBATCH --nodes=%s              ### Number of nodes'%pars['nodes'])
    lines.append('#SBATCH --ntasks-per-node=%s    ### Number of MPI tasks per node'%pars['ntasks_per_node'])
    lines.append('#SBATCH --cpus-per-task=%s      ### Number of HT per task'%pars['cpus_per_task'])
    if pars['gpus_per_node'] is not None:
        lines.append('#SBATCH --gpus-per-node=%s      ### Number of GPUS per node'%pars['gpus_per_node'])
    if pars['gres_gpu'] is not None:
        lines.append('#SBATCH --gres=gpu:%s      ### Value of the --gres=gpu: variable'%pars['gres_gpu'])
    if pars['memory'] is not None:
        lines.append('#SBATCH --mem %s             ### Memory per node'%pars['memory'])
    if pars['time'] is not None:
        lines.append('#SBATCH --time %s         ### Walltime, format: HH:MM:SS'%pars['time'])
    if pars['partition'] is not None:
        lines.append('#SBATCH --partition %s'%pars['partition'])
    if pars['account'] is not None:
        lines.append('#SBATCH --account %s'%pars['account'])
    if pars['qos'] is not None:
        lines.append('#SBATCH --qos %s'%pars['qos'])
    lines.append('#SBATCH --job-name=%s'%job)
    lines.append('#SBATCH --output=%s.out'%job)
    lines.append('')
    lines.append('export OMP_NUM_THREADS=%s'%pars['omp_num_threads'])
    if pars['omp_places'] is not None:
        lines.append('export OMP_PLACES=%s'%pars['omp_places'])
    if pars['omp_proc_bind'] is not None:
        lines.append('export OMP_PROC_BIND=%s'%pars['omp_proc_bind'])
    lines.append('')
    lines.append('echo "Cluster name $SLURM_CLUSTER_NAME"')
    lines.append('echo "Job name $SLURM_JOB_NAME "')
    lines.append('echo "Job id $SLURM_JOB_ID"')
    lines.append('echo "Job nodelist $SLURM_JOB_NODELIST"')
    lines.append('echo "Number of nodes $SLURM_JOB_NUM_NODES"')
    lines.append('echo "Number of tasks $SLURM_NTASKS"')
    lines.append('echo "Number of tasks per node $SLURM_TASKS_PER_NODE"')
    lines.append('echo "Number of threads per task $SLURM_CPUS_PER_TASK"')
    lines.append('echo "OMP_NUM_THREADS : $OMP_NUM_THREADS"')
    lines.append('echo "# Info GPU"')
    lines.append("export NUM_GPUS=$(echo $SLURM_JOB_GPUS | tr ',' '\n' | wc -l)")
    lines.append('echo "GPUs allocated (total): $NUM_GPUS"')
    lines.append('')
    lines.append('echo " "')
    if pars['pre_processing'] is not None:
        with open(pars['pre_processing']) as f:
            for l in f:
                lines.append(l)
    lines.append('')
    lines.append('echo " "')
    lines.append('echo "############### End of the header section ###############"')
    lines.append('echo " "')
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
        mpi_run += ' ' + '-np %s'%pars['mpi']
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
        gpus_per_node (:py:class:`int`) : slurm gpus-per-node variable
        gres_gpu (:py:class:`int`) : value of the --gres=gpu slurm variable
        memory (:py:class:`string`) : slurm mem variable
        time (:py:class:`string`) : slurm time variable, format 'HH:MM:SS'
        partition (:py:class:`string`) : slurm parition variable
        account (:py:class:`string`) : slurm account variable
        qos (:py:class:`string`) : slurm qos variable
        omp_places (:py:class:`string`) : the OMP_PLACES option, can be `cores` or `socket`
        omp_proc_bind (:py:class:`string`) : the OMP_PROC_BIND option, can be `close` or `spread`
        map_by (:py:class:`string`) : the mpi unit for the --map-by option of mpirun
        pe (:py:class:`int`) : number of `processing elements` in the --map-by:unit:PE=n option of mpirun
        rank_by (:py:class:`string`) : the unit for the --rank-by option of mpirun
        pre_processing (:py:class:`string`) : name of the file with pre-processing actions peformed by
            the script before running the computation. For instance, it can be used to load the module
            needed by the running applications

    """

    def __init__(self,scheduler='direct',omp_num_threads=os.environ.get('OMP_NUM_THREADS', 1),mpi=1,
                nodes=1,ntasks_per_node=1,cpus_per_task=1,gpus_per_node=None,gres_gpu=None,memory=None,
                time=None,partition=None,account=None,qos=None,omp_places=None,omp_proc_bind=None,
                map_by=None,pe=1,rank_by=None,pre_processing=None):
        if scheduler == 'direct':
            rules = dict(mpi=mpi,omp_num_threads=omp_num_threads)
            dict.__init__(self,scheduler=scheduler,**rules)
        if scheduler == 'slurm':
            rules=dict(nodes=nodes,ntasks_per_node=ntasks_per_node,cpus_per_task=cpus_per_task,
            omp_num_threads=omp_num_threads,gpus_per_node=gpus_per_node,gres_gpu=gres_gpu,memory=memory,
            time=time,partition=partition,account=account,qos=qos,omp_places=omp_places,
            omp_proc_bind=omp_proc_bind,map_by=map_by,pe=pe,rank_by=rank_by,pre_processing=pre_processing)
            dict.__init__(self,scheduler=scheduler,**rules)
