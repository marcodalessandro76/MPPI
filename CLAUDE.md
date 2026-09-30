# MPPI — notes for Claude Code

The user (Marco D'Alessandro, the author) writes in Italian: reply in Italian. Code, docstrings and commit messages stay in English.

## What the package does
Python interface to run and post-process QuantumESPRESSO and Yambo computations.
Flow: `InputFiles` (PwInput/PhInput/YamboInput) → `Calculators` (QeCalculator/YamboCalculator, both subclasses of
`Runner`: pre_processing → process_run → post_processing; `RunRules` sets MPI/OMP and the scheduler `direct` or `slurm`)
→ `Datasets` (many runs in parallel, `fetch_results`, `seek_convergence`) → `Parsers` (PwParser on
data-file-schema.xml, Yambo parsers on o-* files and netCDF ndb.*; `YamboParser` aggregates them) → analysis in
`Optics` (linear response, χ⁽ⁿ⁾ single frequency / frequency mixing), `Models` (Gaussian pulses, two-level system),
`Utilities` (Dos, BandStructure, FT, lattice utils, constants).
Tutorials (the de-facto documentation and integration tests) are the notebooks in `sphinx_source/tutorials`.

## Machines and workflow
- **Windows laptop** (this repo at `D:\Projects\Research\MPPI`): no QE/Yambo. Python with deps:
  `C:/Users/Marco/miniconda3/python.exe` (the Git Bash `python` has no numpy). Work here on pure-python code
  (inputs, parsers, analysis) and test it against `sphinx_source/tutorials/Reference_data/`.
- **Linux clusters** (reached by ssh, slurm): QE and Yambo installed. Needed to validate Calculators, Datasets and
  to re-run the notebooks. Never run heavy computations on login nodes: use `RunRules(scheduler='slurm', ...)`.
- Git (GitHub `origin`, branches `master` and `devel`) is the only sync channel between machines. Commit and push
  from one machine, pull on the other. Claude's local memory is per-machine: anything that must survive goes in this
  file.
- **Cluster `ismhpc`** (host in `~/.ssh/config`, jump through `narro`): CentOS 7 / glibc 2.17, so VS Code Remote-SSH
  does NOT work there. Claude runs on the laptop and drives the cluster with
  `ssh -o BatchMode=yes -o ClearAllForwardings=yes ismhpc '...'` (ClearAllForwardings avoids the LocalForward 4444
  clash). There: repo `~/Applications/MPPI` (pip editable install, it is the copy the user's notebooks use), miniconda
  python 3.13, `pw.x` (qe-7.0) and `yambo` (Lumen fork 2.1.0) in PATH, slurm. Edits are made on the laptop, pushed,
  then `git pull` on the cluster (the cluster never commits).
  Environments: the user's `pw.x` needs Intel MPI + MKL, yambo needs OpenMPI, so they cannot share one shell env.
  The user keeps module files in `~/module_script/` (`qe_module`, `yambo_module`) and passes them as
  `RunRules(pre_processing=...)`: included in slurm scripts and (since fix/bugs) sourced before `direct` runs.
  Without them pw.x fails with `libmkl_gf_lp64.so: cannot open shared object file`. Small direct runs (mpi 2-4,
  omp 1) on the login node `frontend` are OK for the tutorials; slurm partition `debug` (2h) for short test jobs.
  The user's production RunRules (use them for the slurm runs of the tutorials; ntasks_per_node*cpus_per_task = 32
  cores per node), always with `partition='all12h', memory='125000'` and `activate_BeeOND=True`:
  - QE: `time='11:59:00', ntasks_per_node=16, cpus_per_task=2, omp_num_threads=2,
    pre_processing='/home/dalessandro/module_script/qe_module'`, `QeCalculator(rr, activate_BeeOND=True)`
  - Yambo: `ntasks_per_node=32, cpus_per_task=1, omp_num_threads=1,
    pre_processing='/home/dalessandro/module_script/yambo_module'`,
    `YamboCalculator(rr, executable='yambo_nl', activate_BeeOND=True)`
- Tests: `python -m pytest tests` (fixtures in `tests/conftest.py`). Tests that need the executables are marked
  `@pytest.mark.requires_qe` / `@pytest.mark.requires_yambo` and are skipped when `pw.x`/`yambo` are not in PATH, so
  the same suite runs on both machines. `Reference_data/nl_results` (~1.3 GB) is not in git: the tests that use it
  are skipped when it is missing. Every bug fix gets a test.

## Re-running the tutorial notebooks on the cluster
The notebooks are executed IN PLACE in `~/Applications/MPPI/sphinx_source/tutorials` on ismhpc: their results
(`<Class>_tutorial/` folders, ignored by git via `*_tutorial/`) stay there because later tutorials (e.g. Yambo)
reuse them. Results never go to git or to the laptop.
1. Write/modify the notebook on the laptop (in the scratchpad, NOT in the laptop repo) and `scp` it into the
   cluster `tutorials/` folder. Before running, `rm -rf` that tutorial's `<Class>_tutorial/` folder (clean run).
2. `cd ~/Applications/MPPI/sphinx_source/tutorials && jupyter nbconvert --to notebook --execute --inplace
   --ExecutePreprocessor.kernel_name=python3 --NotebookClient.record_timing=False <nb>.ipynb`.
   Light notebooks run on the login node (direct runs with mpi 2-4, omp 1 are OK); slurm runs use the production
   RunRules above.
3. The executed notebooks stay modified and uncommitted on the cluster until the end of the review. MPPI code
   changes still go laptop → commit/push → `git pull` on the cluster (the pull works as long as the laptop commits
   do not touch those notebooks, so do NOT commit notebooks from the laptop meanwhile).
4. At the end of the review: commit the notebooks on the cluster, push, `git pull` on the laptop.

The user wants clean runs from zero, never reusing the outputs of old runs. Each tutorial starts with a cell
printing the execution date, is kept short (show the main features of the class, no exhaustive tour), writes its
files in a `<Class>_tutorial/` folder, and avoids `obj.method?` cells (nbconvert does not capture the pager: use a
markdown pointer instead). The old run folders `QeCalculator_test/` and `Si_gs_convergence/` were deleted from the
cluster copy on 2026-09-30.

Tutorial status: Tutorial_PwInput, Tutorial_QeCalculator and Tutorial_PwParser (parses the QeCalculator_tutorial results) rewritten (shorter) and executed on ismhpc
(QeCalculator: direct runs with mpi=2 on the login node + one slurm job on all12h with BeeOND). Still to review and run:
all the others. Analysis_BandStructure still uses the old `build_kpath` (now `mppi.Calculators.Tools.build_pw_kpath`).

## Conventions
- Match the existing style: classes that inherit from `dict`, Sphinx-style docstrings with `:py:class:` types,
  `verbose` flags.
- Keep the public API backward compatible (the user's notebooks and scripts on the clusters depend on it). Flag it
  explicitly when a change breaks it.

## Status: bug-fix checklist
Found in the analysis of 2026-09-30. Tick an item once it is fixed AND tested. "verified" = reproduced by running
the code; the other items come from reading it.

Runnable on the laptop (done on branch `fix/bugs`, covered by `tests/`):
- [x] `PhInput.store` calls `self._slicefile`, but the method is `slicefile` → parsing any file fails
- [x] `Dataset.fetch_results` matched ids by substring (`{'ecut':4}` selected `ecut_40`). Now uses `id_matches`:
      dict ids are compared as (key,value) subsets, other ids by '-'-separated tokens of the name. Post-processing is
      called once
- [x] `Utils.dict_merge` used `collections.Mapping` (removed in py3.10)
- [x] `YamboOutputParser`: added the `orbt_magnetization` key (the extension Yambo actually writes)
- [x] `Dos`: default `broad_kind` was the function `lorentzian` → no broadening. Now the string (a callable is also
      accepted)
- [x] `Dos`: `get_evals(set_gap,set_direct_gap)` was positional → the gap was used as a scissor
- [x] `Xn_*.from_file`: `cls(data,verbose)` passed verbose as X_order
- [x] `Xn_frequency_mixing.check_harmonic_reliability` used the stale `ifreq`; debug prints removed
- [x] `Utils.file_parser` ignored its `skip` argument
- [x] `YamboInput.read_file`: undefined `filename` + `exit()` → raises IOError
- [x] `PwInput`/`PhInput` namelists: quoted values were truncated (`prefix = 'ecut:100,k:9'` → `'ecut`). Now
      `parse_namelist_variables` reads quoted strings as a whole, skips `!` comments; namelist names case insensitive
- [x] `PwInput`: raises ValueError instead of `exit()`; parses if_pos flags (stored as a 3rd element of the atom and
      written back), `K_POINTS gamma` and `K_POINTS type` without braces; files closed
- [x] `PwParser`: also catches FileNotFoundError/ParseError; lsda reads nbnd_up+nbnd_dw (with a warning: the gap and
      transitions methods do not support lsda yet, see the Dos/spin item of the roadmap)
- [x] `YamboDipolesParser.get_info` printed dip_v for dip_p; netCDF files closed in all the parsers
- [x] `BandStructure.get_high_sym_positions` passed allclose atol/rtol swapped
- [x] `fit_sum_frequencies`: IndexError when lstsq returned an empty residual
- [x] numpy 2: `np.array(ncVariable)` → `np.array(ncVariable[:])` (DeprecationWarning); invalid escape sequences
      in regexes and LaTeX docstrings (SyntaxWarning in py3.12+) fixed with r-prefixes
- [ ] `MergeQPndb.merge_qp` does not close its input Datasets (do it together with the YamboQPParser work)
- [ ] `Constants.high_sym_fcc['K'] = [0,1,1]` is equivalent to X (the K point is (3/4,3/4,0)); the path X→(0,1,1)→Γ
      used in Analysis_BandStructure passes through K anyway. Ask the user before changing it, notebooks depend on it
- [ ] Physics to check with the user: difference-frequency field not conjugated in `Xn_frequency_mixing.eval_Ew`;
      dephasing 12/damp vs 6/damp in the two Xn classes; missing `dt` and t0 phase in `LRoptics`
- [ ] `NLanalysisYamboPy.py`: many latent NameErrors; decide whether to fix it or remove it

Need the cluster (QE/Yambo/slurm):
- [ ] `Dataset.run_the_calculations` and `Utilities.Parallel.loop` use nested functions with multiprocessing →
      they fail under `spawn` (Windows/macOS, and the Linux default from Python 3.14)
- [ ] slurm `run_ended`: waits forever if the job dies before `JOB_DONE`; IndexError on an empty .out; no
      squeue/sacct check
- [ ] `run_job`: `job` is undefined for an unknown scheduler; `os.system('rm ...')` breaks on paths with spaces
- [ ] `RunRules`: default `omp_num_threads` is read from the environment at import time

## Roadmap after the bug fixes
Refactoring: a common base class for the two calculators, `logging` instead of print, exceptions instead of
print+return, a common base for the two Xn classes. Features (see also `Todo-list.md`): complete YamboQPParser and
MergeQPndb, spin in Dos and `Dos.from_Yambo`, k-point expansion with weights, `update_from_remote` (rsync), Hubbard
support in PwInput, ph.x output parser, BSE exciton parser, `pyproject.toml`.
