import os
from mppi import Calculators as C
import pytest
from mppi.Calculators.RunRules import direct_command, build_slurm_header, mpi_command, environment_info

def test_runrules_direct():
    rr = C.RunRules(mpi=2,omp_num_threads=1)
    assert rr == {'scheduler':'direct','mpi':2,'omp_num_threads':1,'pre_processing':None}
    assert mpi_command(rr) == 'mpirun -np 2'
    assert direct_command(rr,'run','pw.x') == 'cd run ; pw.x'

def test_direct_command_pre_processing(tmp_path):
    env = tmp_path/'env.sh'
    env.write_text('module load qe\n')
    rr = C.RunRules(mpi=2,pre_processing=str(env))
    # the output of the pre_processing file (e.g. module list) is discarded
    assert direct_command(rr,'run','pw.x') == 'source %s > /dev/null 2>&1 ; cd run ; pw.x'%os.path.abspath(str(env))

def test_slurm_header_pre_processing(tmp_path):
    env = tmp_path/'env.sh'
    env.write_text('module purge\nmodule load qe\n')
    rr = C.RunRules(scheduler='slurm',ntasks_per_node=4,partition='debug',pre_processing=str(env))
    lines = build_slurm_header(dict(rr,name='test'))
    assert '#SBATCH --partition debug' in lines
    assert '#SBATCH --job-name=job_test' in lines
    assert 'module purge' in lines and 'module load qe' in lines
    assert mpi_command(rr) == 'mpirun -np 4'

@pytest.mark.skipif(not os.path.isfile('/bin/bash'),reason='needs /bin/bash')
def test_environment_info(tmp_path):
    env = tmp_path/'env.sh'
    env.write_text('echo this output is discarded\nexport PATH=%s:$PATH\n'%tmp_path)
    exe = tmp_path/'my_code.x'
    exe.write_text('#!/bin/bash\n')
    exe.chmod(0o755)
    info = environment_info(C.RunRules(mpi=2,pre_processing=str(env)),'my_code.x')
    assert 'pre_processing : %s'%env in info
    assert 'executable : %s'%exe in info
    assert 'discarded' not in info
