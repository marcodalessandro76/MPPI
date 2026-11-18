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

## Related project
The research work on LiF (pump-probe, yambo_nl non-linear chi vs RT transient absorption) moved on 2026-10-05 to its
own repository and sessions: `D:\RICERCA\DFT AND MANY BODY\SIMULATIONS\LiF` (laptop) and `~/work/LiF` (ismhpc),
GitHub `marcodalessandro76/LiF`, with its own CLAUDE.md. Keep this file for the development of MPPI.

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
  Python env on ismhpc: ~/miniconda3 base, python 3.13; the Jupyter stack, numpy, scipy, matplotlib, netCDF4 are
  pip-installed there (update them with pip, not conda, to avoid duplicated copies). conda is 26.3.2: newer conda
  (>=26.5.2) needs glibc 2.28 (CentOS 7 has 2.17). anaconda-anon-usage 0.8.1, Anaconda ToS accepted (2026-10-02).
  Laptop: conda 26.9.0, JupyterLab 4.6.4 from conda (defaults, ToS accepted). After a conda update check `conda
  info`: an old anaconda-anon-usage plugin breaks conda 26 (fix: update only that package, see git log/notes).
  `export OMPI_MCA_btl=^openib` is set in `~/.bashrc` and in `yambo_module`: it only hides the Open MPI
  "error initializing an OpenFabrics device" warning (openib is never used for IB: btl_openib_allow_ib=false,
  the PML is UCX). Yambo `reformat=True` (header in the inputs) must stay the default.
  Environments: the user's `pw.x` needs Intel MPI + MKL, yambo needs OpenMPI, so they cannot share one shell env.
  The user keeps module files in `~/module_script/` (`qe_module`, `yambo_module`) and passes them as
  `RunRules(pre_processing=...)`: included in slurm scripts and (since v1.3) sourced before `direct` runs,
  with its output discarded (the module files end with `module list`: the user does not want it printed at every
  run). `code.show_environment()` prints the loaded modules and the executable on request.
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
  the same suite runs on both machines. `Reference_data/nl_results` contains only LiF-delta_pulse (tracked in git): the heavy LiF
  runs (~1.3 GB, never in git) were deleted on 2026-10-02, the Optics module is tested on the anharmonic oscillator. Every bug fix gets a test.

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

Tutorial status (on `master`, version 1.3; the work started on the branch `fix/bugs`, merged in master on
2026-10-01; the previous master is saved in the branch `v1.2`):
- DONE, rewritten (shorter) and executed on ismhpc: Tutorial_PwInput, Tutorial_QeCalculator (direct runs with
  mpi=2 on the login node + one slurm job on all12h with BeeOND), Tutorial_PwParser (parses the
  `QeCalculator_tutorial/` results, which stay on the cluster) — committed on 2026-09-30. Tutorial_Datasets
  (2026-10-01, QE only: ecut dataset, post-processing, fetch_results, seek_convergence on k points, a slurm dataset)
  — executed in place on the cluster, NOT committed yet. The old Yambo HF dataset part of Tutorial_Datasets was
  dropped: show a Yambo dataset in the Yambo tutorials.
- Tutorial_YamboInput rewritten and executed (2026-10-01): run_dir `Yambo_tutorial` (SHARED by all the Yambo
  tutorials, they use the same SAVE from `QeCalculator_tutorial/out_nscf`) built with
  `Tools.init_yambo_dir(yambo_dir, input_dir)` (the current API; the MoS2 notebooks use the older
  `make_p2y(source_dir)` + `init_yambo_run_dir`, removed in March 2023). YamboInput.py was rewritten (line based
  parser, readable writer) and checked against 14 yambo/ypp/yambo_nl/yambo_rt inputs (tests/data/yambo_inputs).
  With `reformat=True` yambo keeps only the variables of the `-V` verbosity of args.
- IMPORTANT for every Yambo tutorial: the nscf used for p2y must be run with `force_symmorphic=True` (since
  2026-10-01 it is the default of PwInput and of set_scf/set_nscf/set_bands, it was False before). With
  non-symmorphic symmetries in the SAVE this Lumen yambo silently activates no runlevel (generated inputs contain
  only setup variables, runs end with a report named `r-..._ypp`). Tutorial_QeCalculator now does so.
  The new p2y always exits with MPI_ABORT after writing the wavefunctions, but the SAVE works.
- Tutorial_YamboCalculator rewritten and executed (2026-10-01) in `Yambo_tutorial`: HF direct runs (name/jobname,
  skip, clean_restart replicas), GW ppa, ypp bands, a Yambo Dataset (HF gap vs EXXRLvcs with PP.yambo_get_gap)
  and a slurm run with the production RunRules and BeeOND. HF direct gap of Si at Gamma: 7.9501 eV.
- Tutorial_YamboParser rewritten and executed (2026-10-01): parses the `Yambo_tutorial` results (GW ppa with
  ExtendOut, HF; no band structures: they belong to Analysis_BandStructure) and `Reference_data/rt_results` (RT o- files, ndb.RT_carriers). GW direct gap of Si
  at Gamma 3.32 eV (DFT 2.57). Notes: the yambo vector alat of a fcc cell is alat/2, so the `rescale=True`
  lattice of YamboDftParser differs by 2 from the PwParser one (equal in a.u.); in ypp bands do NOT set the
  `BANDS_path` labels (ypp then uses its own high-symmetry points, wrong for this cell): use `BANDS_kpts`. Use
  `INTERP_mode='BOLTZ'` for smooth bands (the default NN gives step-like bands; the `GfnQP_*INTERP*` variables
  only act on the interpolation of QP corrections from a database, they do not change the DFT bands).
- Analysis_BandStructure rewritten and executed (2026-10-01) in `BandStructure_tutorial`: GaAs, GaAs with
  spin-orbit (rel pseudos, 16 bands), graphene (smearing, G-M-K-G), GaAs-SO with ypp BOLTZ (works, the old
  'interpolation error' issue is gone). Path L-G-X-W-K-G from the corrected Constants.high_sym_fcc.
- Analysis_Dos rewritten and executed (2026-10-01) in `Dos_tutorial`: Si nscf 12x12x12, DOS normalization
  (2 states per band, electrons up to the gap), lorentzian vs gaussian, set_gap, JDOS from get_transitions,
  generic levels and rescale.
- Tutorial_YamboNLDBParser rewritten and executed (2026-10-02): parses `Reference_data/nl_results/LiF-delta_pulse`
  (tracked in git, 756 KB).
- Analysis_Optics rewritten (2026-10-02) on the analytical anharmonic oscillator (no LiF data): linear response,
  single frequency, frequency mixing 2nd and 3rd order vs Boyd's formulas, third harmonic and vanishing even orders of the
  centrosymmetric oscillator (|P_even/P(w)| ~1e-14 single frequency, ~1e-9 mixing: the noise floor of the harmonic fit;
  the same check applies to centrosymmetric crystals such as LiF, where chi(1,+-1) must vanish). The old LiF version (in git history) is meant
  to move to the LiF repository (`~/work/LiF` on ismhpc, github marcodalessandro76/LiF), where the pump-probe
  analysis chi vs RT transient absorption continues.
- Model_AnharmonicOscillator (2026-10-02, new, pure python: written and executed on the laptop): the intro explains the
  two independent parts of the class (numerical P(t) by RK4 of the full equation; analytical chi from perturbation
  theory, Boyd's chi without D). Sections: chi1, transient e^{-gamma t}, harmonics and inversion symmetry (FFT),
  numerical harmonics vs analytical chi (centrosymmetric: chi2=0, even components ~1e-9 of the linear one), validity of
  the perturbative regime (E0^n scaling, fields > ~1e-2 au escape the cubic well). The sections on
  the resonances of chi2 (Miller's rule) and on the signed frequencies were removed at the user's request. When a file is added, add its rst page / notebooks.rst entry to sphinx_source, and the link to a new notebook in the "Tutorials and examples" list of README.md.
- **NEXT**: Analysis_Electron-phonon, Analysis_FourierTransform and Model_TLS_optical_absorption. They can reuse the QE results in `QeCalculator_tutorial/` on the cluster (the nscf in
  `out_nscf/si_scf.save` for p2y). Yambo runs need `pre_processing='/home/dalessandro/module_script/yambo_module'`.
  Lumen was rebuilt on 2026-10-01 (`~/Applications/Lumen`: sources in `src`, build in `gpl-gcc_10.2` from its
  `config_file`, libraries in `lumen-libs`): core, nl-project and rt-project compiled. PETSc 3.24 needs
  `module load cmake-3.26.3` (system cmake is 2.8.12); when a library build fails yambo still writes its
  `*.stamp` files in `gpl-gcc_10.2/lib/<lib>/`, so remove them (and the extracted source dir) before rebuilding.
  The home quota is 19.5 GB (it filled up once during the build). The user wants to give instructions before the
  first yambo tests: ask before running p2y/yambo. The old Tutorial_YamboInput uses the removed `U.build_SAVE` (now
  `mppi.Calculators.Tools.init_yambo_dir`) and the nonexistent `set_GbndRange`/`set_BndsRnXp`.
- Still to review and run after Yambo: the Analysis_* notebooks,
  Model_TLS_optical_absorption.
- Work on `master` (version 1.3). The branches v1.0, v1.1, v1.2 keep the old versions. Do NOT delete the branch
  `fix/bugs`: the user keeps it for future rounds of fixes like this one (to be merged again into master).

## Conventions
- Match the existing style: classes that inherit from `dict`, Sphinx-style docstrings with `:py:class:` types,
  `verbose` flags.
- Spell check: `cspell.json` (VS Code Code Spell Checker) checks only comments and docstrings of the python code
  (also in notebook code cells), never the code. Add new technical terms (Yambo/QE variables...) to its `words`
  list, and write comments/docstrings without typos.
- Keep the public API backward compatible (the user's notebooks and scripts on the clusters depend on it). Flag it
  explicitly when a change breaks it.

## Status: bug-fix checklist
Found in the analysis of 2026-09-30. Tick an item once it is fixed AND tested. "verified" = reproduced by running
the code; the other items come from reading it.

Runnable on the laptop (done in v1.3, covered by `tests/`):
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
- [x] `Constants.high_sym_fcc` was inconsistent (K=(0,1,1) equivalent to X, W not on the face of X). Fixed on
      2026-10-01 (user approved): X=(0,1,0), W=(0,1,1/2), K=(0,3/4,3/4), U=(1/4,1,1/4), L=(1/2,1/2,1/2), the same
      points of high_sym_fcc_crystal (which was already correct; U added)
- [x] Optics physics (2026-10-02, user approved, checked against Boyd *Nonlinear Optics* ch. 1 and the analytical
      anharmonic oscillator of `mppi.Models.AnharmonicOscillator`, tests in `tests/test_optics.py`):
      `eval_Ew` now uses the signed field amplitudes E(-w)=E(w)^* (`field_amplitude`), so the keys with negative orders
      (e.g. (1,-1), (1,-2)) and the zero-th order of `Xn_single_frequency` (now / |E|^2) changed w.r.t. v1.3 and YamboPy
      (YamboPy does not conjugate: its odd negative orders have the opposite sign). The susceptibilities include Boyd's
      degeneracy factor D (documented in the docstrings, user wants only the docs): chi(1,+-1)=2chi2, chi(1,+-2)=3chi3,
      chi0=2chi2(0;w,-w), third order part of (1,0) = 6chi3(w;w,wP,-wP)|E_P|^2 with |E(wP)|^2=E_P^2/4 (in the LiF
      notebook x11m1 must be divided by EP**2/4, not by Ew[(0,2)]). Dephasing 12/damp also in `Xn_frequency_mixing`
      (6/damp gave ~3% errors on the 3rd order keys). `Linear_Response`: P(w) = dt*sum P(t)e^{iwt} (yambo writes the
      DELTA field as one step of value E0/dt); the old factor 2 was ~dt only for the 0.05 fs IO step (eps-1 changes by
      dt/2). The comparison with YamboPy is not important for the user: deviate from it when it is not standard

Need the cluster (QE/Yambo/slurm):
- [x] `Dataset.run_the_calculations` and `Utilities.Parallel.loop` used nested functions with multiprocessing →
      failed under `spawn`. Now module-level `calculator_run`/`func_loop`; a failing run returns None; results
      collected before join. Tested with spawn (laptop) and fork (cluster). Note: the `verbose` (and the other
      global options) of a Dataset are passed to the runs, so they override the verbose of the calculator
- [x] `PostProcessing.pw_get_energy`/`pw_get_gap` crashed on a failed run → None

### OPEN BUGS (to fix in a future round)
Checked again in the code on 2026-11-15 (commit a7f4e25): all still present. The user postponed them; fix them on
`fix/bugs` (then merge into master), each with a test, ticking them here.

Runnable on the laptop:
- [ ] `RunRules.__init__` (`mppi/Calculators/RunRules.py`, signature): the default
      `omp_num_threads=os.environ.get('OMP_NUM_THREADS', 1)` is evaluated once at import time (an `export` done after
      the import is ignored) and it is a string when it comes from the environment. Fix: default None, resolved in the
      body.
- [ ] `run_job` and `run_ended` of `QeCalculator` and `YamboCalculator`: with an unknown scheduler `run_job` prints
      'scheduler unknown' and then `return job` raises UnboundLocalError (`is_ended` is undefined in `run_ended` too).
      Fix: raise a clear ValueError, possibly already in `RunRules`.
- [ ] Shell commands built with `os.system` break on paths with spaces: `rm`/`rm -r` in `clean_slurm_out` and
      `clean_run` of both calculators, and the `rm`/`cd ...; p2y`/`cd ...; yambo` strings of `Calculators/Tools.py`
      (`init_yambo_dir` and the fixsyms helpers). Fix: `os.remove`/`shutil.rmtree`, `subprocess.run(..., cwd=...)` or
      `shlex.quote`.
- [ ] `MergeQPndb.merge_qp` (`mppi/Utilities/MergeQPndb.py`): opens the input netCDF files
      (`datasets = [Dataset(f) for f in files]`) and closes only the output `fout`. Do it together with the
      YamboQPParser/MergeQPndb work of the roadmap.
- [ ] `mppi/Optics/NLanalysisYamboPy.py` (copy of the yambopy nl_analysis kept for comparison): imported nowhere (not
      in `Optics/__init__`, not in tests or notebooks); undefined names `tol` and `residuals_func` in
      `solve_lin_system`, `ifrq` in `update_time_range`, `values` in `reconstruct_signal`; the module docstring imports
      it from `mppi.Utilities` (it is in `mppi.Optics`); it imports `tqdm`, which is not in the requirements of
      `setup.py`. Decide whether to fix it or remove it (the comparison with YamboPy is no longer important for the
      user, see the Optics physics item).

Need the cluster (slurm):
- [ ] slurm `run_ended` (both calculators) and the `wait` loop (1 s polling, no timeout): `wait` loops forever if the
      job dies before writing `JOB_DONE` (time limit, OOM, node failure, scancel) or if `sbatch` itself fails (the
      `Popen` of the submission is never checked); IndexError on an empty `.out` (`f.read().splitlines()[-1]`).
      Fix: read the job id (`sbatch --parsable`), check `squeue`/`sacct` when `JOB_DONE` is missing, optional timeout.

## Roadmap after the bug fixes
Refactoring: a common base class for the two calculators, `logging` instead of print, exceptions instead of
print+return, a common base for the two Xn classes. Features (see also `Todo-list.md`): complete YamboQPParser and
MergeQPndb, spin in Dos and `Dos.from_Yambo`, k-point expansion with weights, `update_from_remote` (rsync), Hubbard
support in PwInput, ph.x output parser, BSE exciton parser, `pyproject.toml`.
The k-point expansion (IBZ -> BZ) of YamboDftParser is done (2026-10-06): `expand_IBZ_kpoints` (attributes
`kpoints_bz`, `kpoints_bz_crystal`, `ibz_index`, `sym_index`, `weights` normalized to one; time reversal from
`DIMENSIONS[9]`, new attribute `time_reversal`) and `get_minus_k_indexes` (index of the point at -k, used by the LiF
project to check the inversion symmetry of a fixsym/NoTr SAVE). Yambo stores the k points in single precision:
equivalent points are found with a tolerance (default 1e-4) on the crystal coordinates, not by rounding. Tests in
`tests/test_parsers.py` (WSe2 12x12x3 with the pw.x weights, rt_results without inversion). Still to do: the same
for PwParser.
Fix (2026-11-14): `YamboDftParser.syms` is now the TRANSPOSE of the netCDF variable SYMMETRY, i.e. the cartesian matrices
dl_sop that Yambo applies to the k points (k' = S k, `k_expand.F` with rl_sop = dl_sop*alat_i/alat_j; netCDF reverses the
dimensions of the Fortran array dl_sop(3,3,nsym)), as in YamboPy (`np.transpose(SYMMETRY,(0,2,1))`). Before, `syms[i]`
was the inverse of the i-th operation of Yambo: the set of the BZ points and the weights were right (the group contains
the inverses: the LiF inversion check is not affected), `sym_index` and the order of `kpoints_bz` change. Checked on the Si
SAVE of the tutorials: in crystal coordinates the Yambo matrices equal the first nsym `PwParser.syms` one by one (p2y keeps
the QE order). New attributes `natoms` and `atomic_positions` (N_ATOMS, ATOM_POS, atomic_numbers; same atoms as PwParser).
`PwParser.syms` contains all the nrot lattice symmetries of the xml (48 for Si), of which only the first nsym (12) are
crystal symmetries (`<info>crystal_symmetry</info>`): new attributes `nsym` and `crystal_syms` (user's request, 2026-11-14;
`syms` unchanged). Test: on WSe2 the crystal_syms equal the YamboDftParser syms in crystal coordinates, in the same order.

## Next work: review of the Optics chi classes (from the LiF project, 2026-10-07)
The LiF project (`NL_Chi/NL-Chi_Analysis.ipynb` in the LiF repository) uses `Linear_Response`, `Xn_single_frequency`
and `Xn_frequency_mixing` on yambo_nl runs with up to 201 probe frequencies (now 155 frequencies, 10-25 eV, 100 fs,
damping 0.3 eV). Points found there, to be reviewed here (with tests on `mppi.Models.AnharmonicOscillator`, whose
susceptibilities are analytic and accept complex frequencies):

1. DONE (2026-10-07). **A posteriori Lorentzian broadening**: `Optics.Utils.lorentzian_broadening(freqs, chi, delta,
   renormalize=True)` (trapezoidal weights, also for non uniform grids; increasing frequencies required) and the options
   `broadening=None` (eV) and `renormalize=True` of `compute_Xn` of both classes (added at the end of the arguments, backward
   compatible), applied only to the keys linear in the probe (key 1 of `Xn_single_frequency`, keys (1,m) of
   `Xn_frequency_mixing`); the other keys are returned unchanged (no error, the reason is in the docstrings: n >= 2 would get
   n*Delta on the n-photon resonance, the zero-th order and the pump-only keys are not analytic/do not depend on the probe).
   The convolution gives chi(w + i Delta): Delta is added to the denominators that contain the probe frequency, the pump-only
   ones are unchanged. The error of the missing tails, ~(Delta/pi)(1/d1+1/d2), affects the whole window (~3% for Delta 0.01
   Ha on 0.30-0.75 Ha), larger if a resonance is close to an edge (10% for (1,2) in Analysis_Optics). Tests on the oscillator
   (exact chi at the complex frequency); a short section "A posteriori broadening" at the end of Analysis_Optics.
2.-4. DONE (2026-10-07, tests in `tests/test_optics.py`). **Outliers**: the default time sampling of
   `Xn_frequency_mixing` started at `Tend - Tw` (Tw = 8*2pi/min distance between the fitted frequencies) even before the
   dephasing time 12/damp, so the free oscillations (not in the fit) spoiled it: error 7.0 on (1,2) at 0.32 Ha, where
   w-3wP = 0.17 is close to 3wP = 0.15. Now it starts at max(Tend - Tw, deph) (back to ~1e-2). **Warnings**: the old ones on
   LiF were spurious, from floating point rounding ((t[-1]-T)+T > t[-1]); times are now compared with a tolerance dt/2,
   and `compute_Xn`/`eval_Pw`/`check_harmonic_reliability` print one summary line per warning (number and indexes of the
   frequencies) via `Optics.Utils.print_sampling_warnings`. **Efficiency**: the fits are stored per frequency
   (`_harmonic_analysis`, key with X_order, Trange, Trange_units, tol, inactive_harmonics): one fit per frequency and
   direction instead of 18; LiF 201 frequencies, X_order (1,3): 13.1 s -> 2.3 s, results identical bit by bit (no LiF
   frequency started before the dephasing time). Still possible: fit the three directions with one lstsq, a
   reliability check for isolated outliers, the common base class of the two Xn classes. Note: `generate_frequencies`
   drops a key whose frequency coincides with another one (e.g. (1,-3) when w = 6wP): the component goes in the kept key.
5. DONE (2026-10-07): note on the INVINT integrator of yambo_nl in the docstring of `Linear_Response` (kick effectively
   at t0 + dt/2 with dt = NLstep, phase w*dt/2 removable by adding dt/2 to efield['initial_time']; red shift
   E_eff = (2/dt) arctan(E dt/2); with NLstep 0.01 fs at 17 eV: shift ~0.1 eV, phase 0.13 rad, not negligible for the
   mixing of Re and Im of eps).

## Exclude nodes in the slurm scripts (from the LiF project, DONE 2026-10-08)
On ismhpc the node wnode07 silently corrupts some yambo_nl runs (see the LiF CLAUDE.md, open point 5). `RunRules` has
the slurm-only parameters `exclude=None` and `nodelist=None` (a slurm node list such as `'wnode07'` or `'wnode[07-08]'`,
or a list of names joined with commas); `build_slurm_header` writes `#SBATCH --exclude=...` / `#SBATCH --nodelist=...`
reading them with `pars.get`, so old RunRules dictionaries keep working and they can also be passed as options of a
single run (`code.run(..., exclude='wnode07')`, the run options update the RunRules). Test in `tests/test_calculators.py`,
notes in the docstrings of RunRules and of the two calculators.

## Transient response chi^neq (from the LiF project; DONE 2026-10-20, commits f0f8add and 8e5df12)
Port of the construction of the pump-probe transient response dchi from the equilibrium chi (theory: LiF repository,
`Notes/main.tex` Secs. 1.4-1.5 and 3.2; reference implementation `Notes/make_figures.py`). Decisions of the user: general
and simple implementation; the corrections of the yambo_nl INVINT integrator stay OUT of MPPI (applied in LiF to the results
of the class, with its `invint_correction`); the derivative term of the QSFIELD pump is neglected; second order terms
included for non-centrosymmetric systems, dropped with `centrosymmetric=True`; user units eV and fs, fields in au.
- `mppi/Optics/TransientResponse.py`: `Transient_Response(freqs (eV), chi1, keys, E0P (au, of the mixing runs), wP (eV),
  centrosymmetric=False, verbose=True)` and `from_Xn(mix, single, broadening=None, renormalize=True, direction=0,
  centrosymmetric=False)` (same probe grid required); attributes `orders` ([2], [3] or [2,3], from the keys available),
  `w` (output grid, at 2wP from the edges, wP for second order only); `components()` (chi_0, chi_+-, chi2_+- on w);
  `dchi(taus (fs), E0P, pump='monochromatic'|'probe time'|'accumulated', sigma (fs), wP=None (eV, phases only), orders=None)`
  -> (w, dchi (n_tau, n_w), au; Delta eps = 4 pi dchi). Sine pump E0P f(t-t_P) sin(wP(t-t_P)), tau = t_p - t_P, signs of the
  notes. The step of the grid must divide m*wP for the keys (1,m) used (ValueError otherwise). `envelope_weight(tau, t,
  sigma, power, model)`: f^p(tau) or the mean of f^p between the probe and t (p = 2 third order, 1 second order).
  `Optics.Utils.time_weighted(freqs, chi, weight)`: int dt e^{iwt} W(t) R(t) with DFTs on the uniform grid (identity for W=1).
- Validation on the oscillator vs the exact dynamics (delta probe + pump, `AnharmonicOscillator.delta_probe_response`,
  new; `solve` has the new optional `kick_index`/`kick_velocity`): monochromatic third and second order within ~1e-2
  (sign check: the opposite signs of the +-2wP terms give 0.3); Gaussian pump with sigma = 2/gamma, third order: accumulated
  2-6% (16% at the tail, tau = +300 au), probe time 26-180%; second order: probe time 7-53%. The accumulated envelope is
  NOT available for the second order terms (ValueError, user's decision 2026-10-20): its natural extension (mean of the
  field envelope f) gave 6-28% on the oscillator, not derived; to be studied if needed. Tests in `tests/test_optics.py`; section "Transient response" at the end of Analysis_Optics
  (before the closing note on the units), sphinx page `transientResponse.rst`.

## Next work: wavefunction parser and 2D visualization (from the 2d-layer-substrate-effective-screening project, 2026-11-14)
Origin: project `~/work/2d-layer-substrate-effective-screening` on ismhpc (GitHub
marcodalessandro76/2d-layer-substrate-effective-screening, own CLAUDE.md, open point "LAYER PROJECTION OF THE FULL
SYSTEM"). There the states of a bilayer (hBN) are projected on the states of one layer, with a prototype reader of the QE
wavefunctions in `hBN_bilayer/Projection_Analysis.ipynb` (cell 2: `read_xml`, `read_wfc`, `gkey`, `overlap`,
`degenerate_groups`; cell 6: `top_state`, rotation/translation/time reversal of the states; validated, see below). The
user wants general tools to read and visualize wavefunctions in MPPI, to understand where a state lives (which layer,
which atoms); the projection itself may become an MPPI module later, built on this class.

Decisions of the user (2026-11-14):
- ONE class with TWO constructors, from the QE output (`<prefix>.save`) and from the Yambo `SAVE`; it holds the data
  and the post-processing methods. It goes in `mppi/Parsers`. Name (agreed): `WfcParser` (`mppi/Parsers/WfcParser.py`,
  `WfcParser.from_qe(save_dir)`, `WfcParser.from_yambo(save_dir)`).
- Band and k indices in the public methods start from 1, as in QE/Yambo and in the user's notebooks ("band 8").
- The class stays GENERIC: one instance = one run. Operations between two runs are generic primitives and must work
  between different instances (e.g. `overlap` of the states of a bilayer and of a monolayer, same k and same cell;
  `transform` of the states of one run). The projection logic (weights P_mn = |<phi_m|psi_n>|^2, averages over
  degenerate groups, layer energies, attribution) stays OUT of MPPI, in the notebooks of the project; a separate
  class for projections may come later, once the procedure is stable.
- Scope of this round (work order points 1 and 2): (1) reading, scalar products, transformations of the states;
  (2) tools for 2D systems: |psi|^2 along z and in the plane, weight of a state in a z window, real space, plots.
  NO k-space tools for now (the user has no use for them yet): implement what the projects need.
- Tests and tutorial use ONLY data produced by the MPPI tutorials (no new heavy runs, no external data); advanced
  tests later with the data of the project.

Formats (checked on ismhpc 2026-11-14 on the hBN monolayer of the project, QE 7.0 + Lumen p2y):
- QE `<prefix>.save/wfc<ik>.dat` (`wfcup/wfcdw` with lsda; `.hdf5` if QE is built with HDF5: not our case, raise a
  clear error), Fortran unformatted, `scipy.io.FortranFile`: record 1 `ik, xk(3), ispin, gamma_only, scalef`
  (`'i4,3f8,i4,i4,f8'`); record 2 `ngw, igwx, npol, nbnd`; record 3 `b1, b2, b3` (reciprocal vectors); record 4 Miller
  indices `(igwx, 3)`; then one record per band, `npol*igwx` complex (double) coefficients (spinor components stacked).
  `xk` and `b` are cartesian in bohr^-1 (checked: k of the IBZ point K equals the yambo one; phases e^{-i(k+G).tau}
  with tau in bohr reproduce the states of a translated run). gamma_only (half of the G stored): raise in the first
  version. Metadata (k points, energies in Ha, lattice, atoms in bohr, nbnd, electrons) from `data-file-schema.xml`
  (reuse PwParser where possible).
- Yambo `SAVE` (netCDF): `ns.db1`: `LATTICE_VECTORS` (3, 3) with the lattice vectors as COLUMNS (bohr);
  `LATTICE_PARAMETER` = the yambo vector alat (alat_x, alat_y, alat_z); `G-VECTORS` (3, nG) and `K-POINTS` (3, nk, IBZ)
  cartesian in units of 2pi/alat_i component by component (so Miller = G_cart . a_j / 2pi, integers within 2e-5);
  `EIGENVALUES` (spin, nk, nb) Ha; `WFC_NG` (nk) and `WFC_GRID` (nk, max ng) = 1-based indices into `G-VECTORS` of the
  components of each k. `ns.wf`: header (`SPIN_VARS`, `FRAGMENTED`, `BAND_GROUPS`); coefficients in
  `ns.wf_fragments_<ik>_<band group>`, variable `WF_COMPONENTS_@_SP_POL<is>_K<ik>_BAND_GRP_<ig>`, shape
  (nb, nspinor, ng, 2) float32 (real, imaginary): single precision. Check: QE vs Yambo states of the same run, G matched
  by Miller indices, |<qe_n|yambo_n>| = 1.000000 at Gamma and K (same states, same phases). Reference implementation:
  yambopy `YamboWFDB` (`~/Applications/YamboPy/yambopy/dbs/wfdb.py`, also rotations to the full BZ and spinor SU(2)
  matrices): use it as a cross-check, not as a dependency.

Data model and API (proposal, MPPI style: dict-like class, Sphinx docstrings, `verbose`): global data (lattice and
reciprocal lattice in bohr / bohr^-1, atoms, k points cartesian and crystal, energies in eV, nk, nbands, nspinor,
electrons, source 'qe'/'yambo'); per k, read LAZILY and only for the requested bands (the save of the hBN bilayer is
several GB): k, Miller indices (ng, 3) int, coefficients (nbands, nspinor, ng) complex. Band and k indices in the
public methods 1-based (see above).
1. `get_states(ik, bands)` -> object with k, Miller indices, coefficients.
2. `overlap(states_1, states_2)`: <phi_m|psi_n> = sum_G c1_m*(G) c2_n(G), G matched by the Miller indices (same k and
   same reciprocal lattice required: check and raise). Norm-conserving pseudopotentials only (no USPP/PAW
   augmentation: warn).
3. `transform(states, rotation, translation, time_reversal)`: psi'(r) = psi(R^-1(r - tau)) gives k' = Rk and
   c'(RG) = e^{-i(Rk+RG).tau} c(G); time reversal: k -> -k, c(G) -> c*(-G). R must map the G set onto itself
   (lattice operation: check Miller integrality). Spinors need SU(2) matrices: NotImplementedError in the first
   version. The project uses C2z + time reversal (k -> k): c'(Gx,Gy,Gz) = e^{-i(k+G).tau} c*(Gx,Gy,-Gz).
4. 2D tools (layer in the xy plane, third cell vector along z):
   - z profile p_n(z) = int_A |psi_n|^2 dx dy = (1/c) sum_{G_par} |sum_{G_z} c_n(G_par, G_z) e^{i G_z z}|^2
     (sum over spinor components; int_0^c p_n dz = 1): no 3D grid needed, one 1D transform per G_par;
   - weight of a state in a z window, int_{z1}^{z2} p_n(z) dz (e.g. the fraction of the state on one layer of a
     bilayer; the project computed it with pp.x + average.x);
   - in-plane map of |psi_n|^2 integrated over a z window (or at one z), on a few cells, atoms and cell drawn;
   - `to_real_space(ik, band, grid)`: psi(r) on an FFT grid (1/sqrt(Omega) normalization); cube writer later;
   - plots (matplotlib): z profile with the atomic planes marked, z shown in [-c/2, c/2) so a layer at z = 0 is not
     split at the cell edge; in-plane map.
   - density of a set of states with k weights (planar average of the total density): later, to compare with the
     density of QE (check with the data of the project).

Tests and tutorial (data of the MPPI tutorials, on ismhpc only, `sphinx_source/tutorials/`):
- `BandStructure_tutorial/out_graphene/graphene.save` (scf, 12x12x1, 19 IBZ k, 8 bands, ecut 60 Ry, graphene at z = 0,
  c = 15 A) and `out_graphene_bands/graphene.save` (bands, 91 k on G-M-K-G): z profiles of sigma vs pi vs pi* states
  (pi: node on the plane), weight within +-1 A of the plane, in-plane maps; transformation check without new runs: a
  C3 rotation about an atom maps K onto K, so the degenerate pi pair at K must be mapped onto itself (sum of
  |<.|.>|^2 over the pair = 2); inversion + time reversal maps K onto K as well.
- QE vs Yambo: Si, `QeCalculator_tutorial/out_nscf/si_scf.save` and the SAVE of `Yambo_tutorial` (made from it):
  overlaps = identity within single precision.
- Spinors: `BandStructure_tutorial/out_gaas_so/gaas_so.save` (npol = 2): reading and normalization.
- Unit tests in `tests/` must not depend on the tutorial outputs (not in git): synthetic coefficients (e.g. one G
  component -> flat z profile; two G_z components -> cosine profile, analytic weights) and a tiny wfc file written in
  the test with FortranFile in the QE format.
- New short tutorial notebook (name after the class), MPPI tutorial rules (execution-date cell, clean run, executed in
  place on ismhpc, rst entry in sphinx_source).
- Later, in the project session: reproduce Projection_Analysis (top-layer states built from ML_bot vs the ML_top runs,
  max deviation 4e-6) and the layer character of DFT_Analysis (AB at K: bands 8 and 10 on the top layer, 98-99%).

Where to work (user, 2026-11-14): the session is opened on the LAPTOP, but this work is done DIRECTLY ON THE CLUSTER
through ssh (`ssh -o BatchMode=yes -o ClearAllForwardings=yes ismhpc '...'`), in `~/Applications/MPPI`, because the
test data (tutorial outputs) are only there: write the code there, run the tests and the tutorial there, commit and
push from the cluster (exception, for this work, to the rule "edits on the laptop, the cluster never commits"). At the
start check that laptop, cluster and origin are at the same commit (`git pull` on the laptop first); do not edit the
laptop copy meanwhile; at the end `git pull` on the laptop. Old git on ismhpc (1.8.3: no `git -C`, use `cd`). Write
files on the cluster with heredocs or scp from the laptop scratchpad.

Status (2026-11-14): point 1 DONE on the cluster (no spin in this round, user's decision: lsda, spinors, gamma_only and
hdf5 files raise NotImplementedError). `mppi/Parsers/WfcParser.py`: module function `read_qe_wfc(filename, bands)`;
`WfcParser.from_qe(save_dir)` (metadata from data-file-schema.xml through PwParser, warning if the type in the
PP_HEADER of a UPF is not NC/SL) and `WfcParser.from_yambo(save_dir)` (metadata through YamboDftParser; from ns.db1
only the G vectors, converted to Miller indices, WFC_NG and WFC_GRID; coefficients from ns.wf or its fragments, band
groups, the first WFC_NG[ik] of the max-ng arrays); attributes lattice and reciprocal_lattice (rows), volume, atoms
(bohr), kpoints (bohr^-1), kpoints_crystal, evals (eV), nkpoints, nbands, num_electrons; `get_states(ik, bands)` (1-based)
-> dict {ik, k, b, mill, bands, c (nbands, ng)}; static methods `overlap(s1, s2)` (G matched by Miller indices, k equal
modulo a reciprocal lattice vector, check of the reciprocal lattice) and `transform(states, rotation, translation,
time_reversal)` (psi(R^-1(r - tau)), conjugate with time reversal; checks orthogonality and Miller integrality). Tests
`tests/test_wfcparser.py` (synthetic only). `Tutorial_WfcParser.ipynb` (executed on ismhpc, in git): Si QE orthonormal
2e-15, QE vs Yambo SAVE = identity within 4e-8 at the 32 k (same phases), energies 1e-6 eV; graphene scf, K = ik 19
(crystal (1/3,-1/3)): C3 about atom A or B and inversion about the bond center + time reversal map each degenerate group
(bands 1-6) onto itself (sums 2.000000/1.000000), C3 about the bond center does not (1.05, 0.004, 0.84, 0.10). The last
computed group (bands 7-8, 1.99997) is excluded.
Point 2 DONE (2026-11-14, on the cluster), static methods on the states (work across runs) unless noted; 2D cell required
(a1, a2 in the xy plane, a3 along +z: `_layer_cell`, ValueError otherwise); lengths in bohr: `z_profile(states, nz)` ->
(z in [-c/2,c/2), p (bands, nz)), one 1D FFT along z per G_par, int p dz = 1; `z_weight(states, z1, z2)` analytic (exact)
weight in a window, also across the cell boundary; `to_real_space(states, grid)` -> psi (bands, n1, n2, n3) on the grid
r = (i/n1)a1 + (j/n2)a2 + (l/n3)a3, sum |psi|^2 dV = 1, default grid 2 max|m_i| + 1, None in the tuple = minimum along
that vector (e.g. (36, 36, None)); `inplane_map(states, z_window, grid)` = sum over the z planes in the window of |psi|^2
times c/n3; instance methods (need the atoms) `plot_z_profile(states, bands, nz, ax)` (planes of the atoms marked) and
`plot_inplane_map(states, band, z_window, grid, ncells, ax)` (pcolormesh on ncells, one cell drawn, atoms folded in the
cell, filled if inside the window). Tests on synthetic states (flat and cosine profiles, analytic weights, real space vs
direct sum, normalizations, plots with Agg). Tutorial, section "Tools for 2D systems": graphene at Gamma, weight within
1 A: sigma bands 1, 3, 4 0.96-0.99, pi band 2 0.80 (node on the plane, maxima at +-0.82 bohr), bands 5-8 vacuum states
0.02-0.05 (their z profiles shown too, user's request 2026-11-15: nearly free electron states over the whole vacuum,
band 7 86% at |z| > c/4, even bands 5 and 7 / odd bands 6 and 8 under z -> -z); maps of sigma (bonds) and pi (atoms, upper half); norms 1 on the (36, 36, 69) grid. Not done yet: density of
a set of states with k weights, cube writer, spin.

## Next work: band interpolation as in ypp BOLTZ (from the 2d-layer-substrate-effective-screening project, 2026-11-17)
NOT started; to do in a future MPPI session (user, 2026-11-17). Goal: a python interpolator of any per-band quantity
f_b(k) known on the IBZ of a regular grid (QP energies, QP corrections dE = E - Eo, energies of layer-projected states)
onto arbitrary k points (a band path), the same as ypp `-s b` with `INTERP_mode= "BOLTZ"`, so that the projects do not
have to write a modified ndb.QP and rerun ypp to get smooth bands.

Why (project, open point LAYER PROJECTION, point 16 of its CLAUDE.md): the GW bands of the top layer extracted from the
bilayer hBN are known only at the 37 IBZ points of the k = 18 grid (18 of them on G-M-K-G) and are plotted joined by
straight lines; the ypp GW bands are smooth. Needs there: interpolate the corrections of bands 2-5 (valence 2-4 known at
all 37 k; band 5 only at 22 of them: support a MASK of missing points, at least to show where the result is not
constrained by data) and add them to the exact KS bands of a pw.x bands run on the path.

Algorithm, as read in the Lumen source on ismhpc (`~/Applications/Lumen/src`, branch bug-fixes, 2026-11-17;
`src/interpolate/INTERPOLATION_BZ_setup.F`, `INTERPOLATION_BZ_coeff.F`, `INTERPOLATION_fourier.F`,
`src/modules/mod_interpolate.F`, `mod_interpolate_tools.F`) = Shankland-Koelling-Wood with the Pickett-Krakauer-Allen
roughness (PRB 38, 2721 (1988); the BoltzTraP scheme):
- Lattice vectors R (integer coordinates on the direct lattice) with |R| <= r, r^3 = INTERP_Shell_Fac * N_k(BZ) *
  3 V / (4 pi) (V = cell volume, bohr^3; default Shell_Fac 20), sorted by |R| and grouped in stars under the symmetry
  operations written on the direct-lattice crystal basis (int_sop, integer matrices from dl_sop). Yambo's nsym
  INCLUDES the time-reversal operations -S (12 = 6 + 6 for the hBN monolayer); MPPI `YamboDftParser.syms` has the
  spatial ones and `time_reversal` separately (add -S), `PwParser.crystal_syms` the crystal ones.
- Star function S_R(k) = (1/n_star) sum_{R' in star(R)} exp(2 pi i k.R'), k in crystal (reduced) coordinates of the
  reciprocal lattice; f(k) = Re sum_R c_R S_R(k).
- Fit, exact at the N data points (the IBZ points), last point N as reference: delta_k = f(k) - f(k_N),
  S~_R(k) = S_R(k) - S_R(k_N) (R != 0), rho(R) = (1 - c1 x)^2 + c2 x^3 with x = |R|^2/|R_min|^2, c1 = c2 = 0.75,
  R_min = first non-zero shell; H_kk' = sum_R S~_R(k') conj(S~_R(k)) / rho(R) (N-1 x N-1), solve H lambda = delta,
  c_R = sum_k conj(S~_R(k)) lambda_k / rho(R) (R != 0), c_0 = f(k_N) - sum_{R!=0} c_R S_R(k_N). After the fit yambo
  logs the max/average error at the data points (`k_compare`), which should be ~0.
- ypp bands with a QP database (`electrons_bands_interpolate.F`, ypp/electrons): if the energies come from a db with
  Eo (GfnQPdb) and a circuit SAVE is given (`CIRCUIT_E_DB_path`, a pw.x bands run through `p2y -nowf`), ypp
  interpolates only dE = E - Eo (`FILL_SEED_E`, l_INTERP_dE) and adds it to the exact KS energies of the circuit;
  the circuit energies are first shifted by delta_align = (interpolated E of the first band at the first circuit
  point) - (circuit KS energy there), which gives the rigid offset found in the project (592.5 meV ML hBN, removed
  there by matching at K). MPPI should not reproduce this alignment: just KS(path) + dE_interp(path).
- 2D cells: yambo takes R inside a SPHERE, so for hBN (c = 40 bohr, r ~ 106 bohr at k = 18) it also includes R with
  out-of-plane components (i3 = +-1, +-2); with k_z = 0 they repeat the in-plane star functions with a different
  roughness (583 stars in the log of the ML run). To reproduce ypp exactly keep them; a cleaner 2D option is in-plane
  R only (circle of radius r in the plane): an option, to compare with ypp.
- Known limit (project, GW_Analysis): the interpolation is per band index; where bands of different character cross
  (pi*/NFE in hBN) dE jumps and the interpolation is off by 0.2-0.3 eV between grid points. Not a bug of the method.

Validation data on ismhpc (project folder `~/work/2d-layer-substrate-effective-screening/hBN_bilayer/MB_analysis/`):
`kpoints_18-nbnds_250-ML_bot` (ML, bands 2-8), `kpoints_18-nbnds_400-AA1`, `kpoints_18-nbnds_400-AB` (bands 3-14):
full-IBZ GW run `pp-...-full_ibz-b<b0>_<b1>/ndb.QP` (37 k; read with YamboQPParser), ypp runs `ypp_bands_ks/` and
`ypp_bands_qp/` (inputs `ypp_bands_*.in`, outputs `o-*.bands_interpolated`, logs `r-*` with the number of shells),
circuit SAVEs `bands_path-<system>` (181 k on G-M-K-G), pw.x bands runs in `hBN_bilayer/DFT_analysis/bands-<system>.save`.
Test: MPPI KS(path) + dE_interp(path) = ypp `bands_interpolated` minus a constant (the alignment), band by band; exact
at the grid points that lie on the path. MPPI tools to reuse: YamboQPParser (ndb.QP), YamboDftParser (IBZ k points,
syms, time_reversal, lattice), PwParser (bands run), BandStructure (plots; `from_Ypp` already reads the ypp output).
API to agree with the user (e.g. an interpolator built from the IBZ k points, lattice and symmetries, with fit(values,
mask=None) and evaluate(kpoints)), with tests on an analytic periodic function and on the hBN data.

## Next work: screening database parser (from the 2d-layer-substrate-effective-screening project, 2026-11-18)
Origin: project `~/work/2d-layer-substrate-effective-screening` on ismhpc (own CLAUDE.md, entry "YAMBO SCREENING
DATABASES"). There a prototype in `hBN_bilayer/Embedding_Analysis.ipynb` (sections 1, 2, 4) reads the Yambo screening
databases, moves a layer by a rotation/translation, builds the W of a layer with a substrate and writes a new database
that yambo then uses for GW. The user wants the GENERAL tools in MPPI, then the project notebook rewritten on top of them.

Decisions of the user (2026-11-18):
- ONE class for a screening database, in `mppi/Parsers`; one instance = one database. Proposed name
  `YamboScreeningParser` (`mppi/Parsers/YamboScreeningParser.py`; to confirm with the user).
- In MPPI: reading, quantities derived from X, Coulomb interaction from `ndb.cutoff`, symmetry transformation, writer,
  2D tools (list below). NOT in MPPI in this phase: the construction of the new W of a layer with a substrate (the
  Dyson equation W = vt + vt chi0^L W of the project notes, Eq. 19, symmetrized as
  X_new = (1 - (1 + X^S) chi0~^L)^-1 (1 + X^S) - 1), the choice of what is layer and what is substrate, the R and tau
  of a given system, the analysis (g(Q), real space, comparisons): it stays in the project notebook and may move to
  MPPI later, once the procedure is stable. The writer takes any X, whatever built it.
- Indices of q in the public methods start from 1, as in yambo (`ndb.pp_fragment_<iq>`) and in the notebooks.

Formats (checked 2026-11-17/18 on the Lumen source and on the hBN runs of the project, Lumen 5.4 r27515):
- Same layout for `ndb.em1s` (static, 1 frequency), `ndb.em1d` (n frequencies) and `ndb.pp` (PPA, 2 frequencies):
  `src/io/io_X.F` (db_desc from X%whoami = 2, 3, 4). Header `ndb.<kind>`: `X_PARS_1` = [ng, first band, last band,
  ...]; `X_PARS_3[3:6]` = "Field direction" q0 (the direction of q -> 0, LongDrXp; |q0| = 1e-5 in the hBN runs:
  (7.07e-6, 7.07e-6, 0) for LongDrXp [1,1,0]); `X_RL_vecs` (3, ng) = the first ng G vectors of the global list of
  `SAVE/ns.db1` (`G-VECTORS`, same order), in yambo units (cartesian, 2pi/alat_i component by component: cartesian
  bohr^-1 = iku * 2pi/alat with the alat VECTOR `LATTICE_PARAMETER` of ns.db1; the header has no alat, so the SAVE is
  needed); `HEAD_QPT` (3, nq) the IBZ q points, same units; `X_TDDFT_KERNEL` (char, 'HARTREE' = RPA); `CUTOFF` (char,
  e.g. 'slab z', 'none'). Read char variables as b''.join(np.ma.filled(var[:], b' ').ravel()).decode().strip().
- Fragments `ndb.<kind>_fragment_<iq>`: `X_Q_<iq>` with netCDF shape (n_freq, ng, ng, 2) float32 (re, im) written
  from the Fortran array X(G, G', w): the two G axes must be TRANSPOSED when read with netCDF4
  (X[w, G, G'] = x[w, G', G]); `FREQ_sec_iq<iq>` (n_freq, 2) in Ha (hBN PPA: 0 and i 27.211 eV).
- Content (`pol_function/X_redux.F`): Dyson X = (1 - chi0 v)^-1 chi0 in the ng block, then multiplied by
  4pi/(bare_qpg(G) bare_qpg(G')): X = v^1/2 chi v^1/2 = eps^-1 - 1 (symmetrized), v = the (cut) Coulomb interaction,
  v(q+G) = 4pi/bare_qpg^2. chi0 is not stored, but it is recovered exactly in the same block (RPA kernel only):
  chi0~ = v^1/2 chi0 v^1/2 = 1 - (1 + X)^-1.
- `ndb.cutoff` (written by the run that computes the Coulomb cutoff; in the project it is in the shared folder
  `ndb_rim_cutoff`, NOT in the screening folder: pass its path explicitly): `GEOMETRY`, `DESCRIPTION` (char, 'slab z'),
  `PARS` ([5] = ng_cut, 12311 for hBN), `HEAD_QPT` (same q), `CUT_BARE_QPG` (ng_cut, nq, 2) float32, real part only
  (imaginary 0 for the slab); its first ng G are those of X (same global order). Slab cutoff (`coulomb/cutoff_slab.F`):
  bare_qpg = |q+G| / sqrt(1 - e^{-Q_par lcut} cos(q_z lcut)), lcut = c/2, Q_par = |q+G|_par regularized to
  >= q0_def_norm; at iq = 1 the |q+G| use q0; without RIM, bare_qpg(1,1) is replaced by an average over a circle
  (RIM_qpg not allocated). Checked on hBN: file = analytic formula within 2e-7 (q != 0), 2e-5 at iq = 1 (with q0).
  => read the file (exact values of yambo, any geometry: slab, box, cylinder), the analytic slab formula only as a test.
  No cutoff ('none'): v = 4pi/|q+G|^2 computed from G and q, with q0 for G = 0 at iq = 1.
- `ndb.RIM_W` (in the same folder as the screening when RIM_W is used) is built by `qp/QP_interpolate_W.F` FROM the X
  of the screening database, and a file with a matching header is reused without checking X: a database with a
  modified X must NEVER carry an old ndb.RIM_W (yambo recomputes it if absent).
- Other yambo screening databases may live in the same folder (e.g. ndb.dipoles, ndb.RIM_W): copy only what is needed.
- Reference implementation (cross-check only, not a dependency): yambopy `YamboStaticScreeningDB`
  (`~/Applications/YamboPy/yambopy/dbs/em1sdb.py`): reads only the static frequency, reads ndb.cutoff, its writer
  `saveDBS` deletes the destination folder (rmtree): do not copy that behaviour.

API (proposal, MPPI style: dict-like class, Sphinx docstrings, `verbose`):
1. Constructor `YamboScreeningParser(folder, save, cutoff=None, verbose=True)`: `folder` with ndb.em1s / ndb.em1d /
   ndb.pp (kind detected, error if more than one or none), `save` = the SAVE (alat, lattice through YamboDftParser),
   `cutoff` = path of ndb.cutoff (warn if CUTOFF != 'none' and no path given). Attributes: kind, ng, bands, kernel,
   cutoff geometry, G (cartesian bohr^-1, and in crystal coordinates), q (IBZ, cartesian bohr^-1), nq, q0, frequencies
   (complex, eV; from fragment 1, check that all fragments agree only on request).
2. `get_X(iq)` -> complex array (n_freq, ng, ng), order [w, G, G'], read lazily, netCDF file CLOSED after reading
   (prototype pitfall: a file left open makes a later 'r+' open in the same process fail with an HDF error).
3. Derived, no v needed: `eps_inv(iq)` = 1 + X; `chi0(iq)` = chi0~ = 1 - (1 + X)^-1 (symmetrized; raise if the
   kernel is not HARTREE). With v: `coulomb(iq)` -> v(q+G) for the ng G (from ndb.cutoff or 3D); `W(iq)` =
   v^1/2 (1 + X) v^1/2 (W of the system the database describes, a.u.); `chi(iq)` = v^-1/2 X v^-1/2.
4. `transform(X, iq, rotation=None, translation=None)`: X of the system moved by r -> R r + tau (R cartesian
   orthogonal, tau in bohr): chi'(r, r') = chi(R^-1(r - tau), R^-1(r' - tau)) gives
   X'_GG'(q) = e^{-i(G-G').tau} X_{R^-1 G, R^-1 G'}(R^-1 q). Implemented cases: R^-1 q = q (direct) and R^-1 q = -q
   (time reversal / Onsager, X_GG'(-q) = X_{-G',-G}(q): X'_GG' = e^{-i(G-G').tau} X_{-R^-1 G', -R^-1 G}(q)); otherwise
   NotImplementedError (X outside the IBZ needed: from the symmetries of the database, later). At iq = 1 use the q0
   DIRECTION, not q = 0 (C2z maps q0 -> -q0: the time-reversal branch; treating iq = 1 as q = 0 gave 5e-3 on hBN AA').
   Checks: the G set closed under R (error naming the missing G), R maps the cutoff geometry onto itself (slab: R z =
   +-z; then v(q+G) is invariant and the symmetrized X transforms as chi). A translation along z is a pure phase.
   Prototype `move_layer` (Embedding_Analysis Sec. 2): hBN ML_bot -> ML_top_AB (R = 1) and ML_top_AA1 (R = C2z, i.e.
   X'_GG' = e^{-i(G-G').tau} X_{sG', sG}, s = (Gx, Gy, -Gz)) vs the computed ML_top databases: relative error
   1.5-2.1e-5 (w = 0) and 1.0-1.3e-5 (i 27 eV) at q != 0, 1.5-1.7e-4 at q -> 0, all 37 q; controls: layer not moved
   1.35-1.41, AA' with the translation only 0.70-0.76.
5. Writer `write(dst, X)` (method of the TEMPLATE instance; X = dict iq -> array (n_freq, ng, ng), or a callable):
   copies the header and all the fragments of the template (only ndb.<kind>* files, NEVER ndb.RIM_W) into `dst`,
   replaces `X_Q_<iq>` (float32, G axes transposed back), refuses dst = the template folder, refuses an existing dst
   with different content unless asked; rewrites a fragment only if it differs by more than a tolerance (prototype:
   1e-7); reads back and checks (float32: ~1.5e-8). The header (bands, ng, cutoff, PPA energy) stays the template's:
   document that yambo then believes the template parameters. Return a summary (fragments written / unchanged).
6. 2D tools (a1, a2 in the xy plane, a3 along +z, as in WfcParser): `inplane_blocks(iq)` -> dict G_par -> indices of
   the G with that in-plane part (sorted by Gz); `planar(M, iq, G_par, z, z2)` -> M(Q; z, z') =
   (1/c) sum_{Gz, Gz'} e^{i Gz z} M_{GG'} e^{-i Gz' z'} on the G_par = G_par' block, Q = q + G_par (convention:
   M(r, r') = (1/V) sum_q sum_GG' e^{i(q+G).r} M_GG' e^{-i(q+G').r'} and in-plane (1/A) sum_Q e^{iQ.rho} M(Q; z, z')).
   Check: the slab-cut v gives (2pi/Q) e^{-Q|z - z'|} for |z - z'| < c/2 and 0 beyond. The prototype `along_z(M, iq,
   idx, z, z2, left)` and `on_top_plane` (Embedding_Analysis Sec. 3) are the same with M = v^1/2 X v^1/2 etc.

Tests and tutorial:
- MPPI tutorial data (ismhpc): `sphinx_source/tutorials/Yambo_tutorial/gw_ppa/ndb.pp` (+ 31 fragments; 27 G, 32 q,
  bands 1-8, kernel HARTREE, cutoff 'none', q0 = (1e-5, 0, 0)) with `Yambo_tutorial/SAVE` (12 symmetry matrices).
  Reading, X Hermitian at both frequencies (2-3e-8), eps_inv, chi0 round trip (Dyson back to X), writer round trip in
  a temporary folder. Symmetry test without new runs (checked 2026-11-18 with the prototype logic): for every q and
  every operation of `ns.db1` SYMMETRY with R^-1 q = +-q, transform(X, iq, R) = X: 64 direct + 64 time-reversal cases,
  max relative deviation 6.4e-7; a lattice translation must give X back too.
- Unit tests in `tests/`: synthetic only (a tiny database written in the test with the same netCDF variables; X with
  known structure for transform and planar; analytic slab v for `planar`).
- Tutorial notebook (name after the class), MPPI tutorial rules (execution-date cell, executed in place on ismhpc,
  rst entry).
- Later, in the project session: rewrite Embedding_Analysis on the class and check that the numbers do not change
  (move_layer errors above; new databases `..._emb_ML_bot` identical; vt/v on the top plane 0.757 at Q = 0.085
  bohr^-1); slab part: ndb.cutoff of the project vs the analytic formula.

Where to work: as for WfcParser, on the cluster in `~/Applications/MPPI` (test data only there): check first that
laptop, cluster and origin are at the same commit; commit and push from the cluster; `git pull` on the laptop at the end.

Status (2026-11-18): DONE on the cluster, name `YamboScreeningParser` (proposed name kept; the user may still rename it).
`mppi/Parsers/YamboScreeningParser.py`, plain class like WfcParser (not a dict). Constructor as in point 1 (`save` = SAVE
folder or its ns.db1; `cutoff` = ndb.cutoff or its folder); extra attributes `syms`, `time_reversal` (from the SAVE),
`nfreqs`, `G_indexes(vectors)`, `get_frequencies(iq)`, `check_frequencies()`, `bare_qpg(iq)`. The derived methods
(`eps_inv`, `chi0`, `W`, `chi`) take an optional `X=None` (a matrix of the user instead of the database). No cutoff:
bare_qpg = |q+G|, at iq = 1 |G| and |q0| for G = 0 (as `col_driver.F`: q_pt(1) = 0, bare_qpg(1,1) = q0_def_norm).
`transform` checks the invariance of v numerically (bare_qpg[R^-1 G] = bare_qpg[G]) when v is available, instead of a
check on the geometry string (warning if a cutoff has no ndb.cutoff). `write(dst, X, overwrite=False, tol=1e-7)`
returns {'written', 'copied', 'unchanged'}; only fragmented databases (NotImplementedError otherwise); with
overwrite=True a different header in dst is replaced (its old fragments removed) and a ndb.RIM_W is removed.
2D tools: `inplane_blocks()` (no iq: the G set is the same at every q; keys = in-plane Miller (m1, m2)),
`inplane_Q(iq, G_par)`, `planar(M, iq, G_par, z, z2)` -> (..., len(z), len(z2)). Tests `tests/test_screeningparser.py`
(synthetic ns.db1, ndb.pp + fragments, ndb.cutoff; slab-cut v in `planar` vs (2pi/Q) e^{-Q|z-z'|} within 2e-3 with
1201 Gz; note: the representation is periodic in z, so "0 beyond c/2" is not testable). `Tutorial_YamboScreeningParser`
(executed on ismhpc, gw_ppa of Yambo_tutorial, writes in `YamboScreeningParser_tutorial/`): X Hermitian 3.6e-7,
eps_M(q->0) = 15.6 (27 G, 8 bands), Dyson round trip 3e-16, 66 direct + 66 time-reversal operations of the SAVE give X
back within 6.4e-7 (max|dX|/max|X|; the earlier count 64 + 64 of the prototype check is not reproduced, irrelevant);
the Si atoms are at +-p (origin = bond center): inversion 1.3e-6, translation atom 1 -> 2 and inversion about an atom
0.54; writer round trip. The 2D tools are only tested on synthetic data (no 2D screening database among the tutorial
outputs: a small graphene em1s run would be needed, ask the user). NEXT: in the project session, rewrite
Embedding_Analysis on the class (numbers listed above).

## Next work: Coulomb parser, and the screening parser without v (from the 2d-layer-substrate-effective-screening project, 2026-11-18)
Decisions of the user (2026-11-18), after the first version of `YamboScreeningParser` (commit 17d1271):
- `YamboScreeningParser` must NOT handle the Coulomb interaction: it mixes two different databases. Everything
  intrinsic to the screening database needs no v, because the matrices are symmetrized: eps^-1 = 1 + X,
  chi0~ = 1 - (1 + X)^-1, transform, write, the 2D tools (`planar` works on any matrix), also eps_M = 1/(eps^-1)_00
  (diagonal elements are the same in the symmetrized and in the physical form). v is needed only to go to physical
  units (chi = v^-1/2 X v^-1/2, W = v^1/2 (1 + X) v^1/2, v_emb in Ha bohr^3): that is a combination of two parsers and
  is done in the user's code, e.g. `s = np.sqrt(coul.coulomb(iq, scr.ng)); W = s[:, None] * (1 + X) * s[None, :]`
  (the v^1/2 inside X are exactly those of ndb.cutoff: X_redux.F multiplies by 4pi/(bare_qpg(G) bare_qpg(G'))).
  => REMOVE from YamboScreeningParser: the `cutoff` argument, `_read_cutoff`, `bare_qpg`, `coulomb`, `W`, `chi`
  (breaking change of a class used so far only by the project notebook `hBN_bilayer/Embedding_Analysis.ipynb`, which
  is updated afterwards in the project session; tutorial and tests to update here).
  `transform` keeps the check that R leaves v invariant, but from the geometry in the header of the screening
  database (`cutoff_geometry`, the CUTOFF of ndb.<kind>) instead of the numbers: 'none' -> any orthogonal R that maps
  the G set onto itself; 'slab <axis>' -> R must map the axis onto +-itself (slab z: R[2, :2] = R[:2, 2] = 0); other
  geometries (box, cylinder, sphere, ws) -> raise NotImplementedError unless the user disables the check explicitly
  (to agree in the MPPI session).
- NEW class for the Coulomb interaction of a run: `YamboCoulombParser` (`mppi/Parsers/YamboCoulombParser.py`, name
  agreed). It reads `ndb.cutoff` and/or `ndb.RIM` (both describe the cut Coulomb interaction of the run, depend only
  on the geometry: cell, k grid, cutoff, RandQpts/RandGvec; a run may have one, the other or none). In the project
  both are in the shared folder `ndb_rim_cutoff`, written by the same `yambo -r` run, but in general they MAY BE IN
  TWO DIFFERENT FOLDERS (user, 2026-11-18): the constructor takes their paths independently (below).

Formats (checked 2026-11-18 on the Lumen source and on the hBN ML_bot run of the project, Lumen 5.4 r27515):
- `ndb.cutoff` (written by `coulomb/cutoff_driver.F`; slab: `coulomb/cutoff_slab.F`): `GEOMETRY`, `DESCRIPTION` (char,
  'slab z'), `PARS` ([5] = ng_cut, 12311 = FFTGvecs for hBN), `HEAD_QPT` (3, nq) iku, `HEAD_R_LATT` ([nqibz, nqbz, ...]
  = [37, 324, 37, 324]), `CUT_BARE_QPG` netCDF (ng_cut, nq, 2) float32 = Fortran bare_qpg(iq, ig) complex (imaginary
  part 0 for the slab): v_G(q) = 4pi / bare_qpg^2. G = the first ng_cut G of the global list of `SAVE/ns.db1`
  (`G-VECTORS`, cartesian bohr^-1 = iku * 2pi/alat_i) -- ndb.cutoff has NO G list: the correspondence relies on the
  global order (cutoff_slab.F loops ig = 1..wf_ng over g_vec). Slab formula: bare_qpg = |q+G| / sqrt(1 - e^{-Q_par lcut}
  cos(q_z lcut)), lcut = c/2, Q_par regularized to >= q0_def_norm (1e-5); at iq = 1 the code uses q = 0 for G != 0
  (|q+G| = |G|) and |q0| for G = 0; without RIM, bare_qpg(1,1) is replaced by an average over a circle (RIM_qpg not
  allocated). Checked: file = analytic formula within 2e-7 (q != 0), 2e-5 at iq = 1.
- `ndb.RIM` (`io/io_RIM.F`, `coulomb/rim.F`, `coulomb/rim_integrate_v.F`): `PARS` (10) = [RIM_ng (101 = RandGvec),
  RIM_is_diagonal (1 for the slab), RIM_n_rand_pts (5e6), RIM_RL_vol (Monte Carlo BZ volume, 0.31984), RL_vol (true BZ
  volume, 0.31980), RIM_id_epsm1_reference, RIM_epsm1 (3), RIM_anisotropy] (eps^-1 anisotropy entries 0 when Em1Anys
  is not set); `HEAD_QPT`, `HEAD_R_LATT` as above; `CUTOFF` (char): UNRELIABLE, it says 'none' in the hBN file although
  the content is slab-cut (checked numerically, below): do not use it; `RIM_qpg` netCDF (RIM_ng, RIM_ng, nq) float32 =
  Fortran RIM_qpg(iq, G, G'): transpose (2, 1, 0) -> [iq, G, G']; with RIM_is_diagonal only the diagonal is filled
  (3737 = 37 x 101 nonzero entries). Meaning: Monte Carlo average over the small cell of the BZ around q (volume
  d3q_factor = RL_vol / nqbz) of 1/b_G(q + q')^2 (slab-cut), instead of the value at q. Normalization
  (`wf_and_fft/scatter_Gamp.F`): with RIM, gamp = RIM_qpg/2 replaces q_weight/(bare_qpg(G) bare_qpg(G')), q_weight =
  d3q_factor/(2pi)^3, so the cell-averaged Coulomb interaction is <v>_cell(q, G) = 4pi RIM_qpg(q, G, G) / (2 q_weight).
  Head at q = 0 (iq = 1, G = 0): the integrable divergence (analytic integral on a small circle for the slab + Monte
  Carlo outside). Checked on hBN (<1/b^2>_cell / (1/b^2) with b from ndb.cutoff): 1.000-1.009 for iq >= 5 or large G,
  1.085 at G = 0 of the smallest q (iq = 2); at iq = 1 the point value is not usable (G = 0: at q0; G along z with even
  n: the slab-cut point value at q0 almost vanishes, the cell average is finite, ratio ~2160): that is what RIM is for.
  Used by Yambo in Sigma_x (scatter_Gamp), in the head of the BSE kernel (K_screened_interaction) and to build
  ndb.RIM_W (QP_interpolate_W); NOT in the screening (X_redux uses bare_qpg). ndb.RIM depends only on the geometry, so
  it can be shared with a modified screening database (unlike ndb.RIM_W, built from X).

API (proposal, to refine in the MPPI session; same style as YamboScreeningParser):
- `YamboCoulombParser(save, cutoff=None, rim=None, verbose=True)`: `save` = the SAVE (or its ns.db1), for alat and the
  G list; `cutoff` = the ndb.cutoff file or the folder that contains it; `rim` = the ndb.RIM file or the folder that
  contains it. The two paths are INDEPENDENT (the files may be in different folders; passing the same folder to both
  is the common case). At least one of the two is required (error otherwise); a path given but without the file ->
  error naming the missing file. If both are given, check that their q points (HEAD_QPT, HEAD_R_LATT) agree with each
  other (same run/geometry). Attributes: geometry (from
  ndb.cutoff; 'none' if only ndb.RIM), q (cartesian bohr^-1), nq, nqbz, G (first ng_cut G of ns.db1, cartesian), ng_cut,
  and for the RIM: rim_ng, rim_is_diagonal, rim_points, RL volumes.
- `bare_qpg(iq, ng=None)` and `coulomb(iq, ng=None)`: point values on the first ng G (v = 4pi/b^2). Without
  ndb.cutoff (no cutoff): 4pi/|q+G|^2 computed from G and q, with |q0| for G = 0 at iq = 1 (q0 norm as an argument,
  default 1e-5 = q0_def_norm of Yambo; as `col_driver.F`).
- `coulomb_rim(iq)`: the cell average <v>_cell for the rim_ng G (vector if diagonal, matrix otherwise).
- `check_slab()`: comparison of bare_qpg with the analytic slab formula above (this also checks the correspondence
  with the G list of the SAVE, which the file does not store).

Work order (user, 2026-11-18): (1) remove v from YamboScreeningParser (and update its tests and tutorial); (2) write
YamboCoulombParser; (3) the user then reviews the new class (API and checks) before the project notebook is updated.

Tests: no ndb.cutoff / ndb.RIM among the MPPI tutorial outputs (Yambo_tutorial is 3D without cutoff and without RIM):
synthetic files in `tests/` (same netCDF variables); tutorial: a small 2D run with slab cutoff and RIM would be needed
(ask the user), otherwise the hBN data of the project in the project session (values above).
Then in the project session: Embedding_Analysis with YamboScreeningParser (no cutoff) + YamboCoulombParser for the
quantities in physical units (Sec. 1 table, Sec. 3 v_emb, Sec. 4 W/v on the top plane); numbers must not change.
