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
  the resonances of chi2 (Miller's rule) and on the signed frequencies were removed at the user's request. When a file is added, add its rst page / notebooks.rst entry to sphinx_source.
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
- [ ] `MergeQPndb.merge_qp` does not close its input Datasets (do it together with the YamboQPParser work)
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
- [ ] `NLanalysisYamboPy.py`: many latent NameErrors; decide whether to fix it or remove it

Need the cluster (QE/Yambo/slurm):
- [x] `Dataset.run_the_calculations` and `Utilities.Parallel.loop` used nested functions with multiprocessing →
      failed under `spawn`. Now module-level `calculator_run`/`func_loop`; a failing run returns None; results
      collected before join. Tested with spawn (laptop) and fork (cluster). Note: the `verbose` (and the other
      global options) of a Dataset are passed to the runs, so they override the verbose of the calculator
- [x] `PostProcessing.pw_get_energy`/`pw_get_gap` crashed on a failed run → None
- [ ] slurm `run_ended`: waits forever if the job dies before `JOB_DONE`; IndexError on an empty .out; no
      squeue/sacct check
- [ ] `run_job`: `job` is undefined for an unknown scheduler; `os.system('rm ...')` breaks on paths with spaces
- [ ] `RunRules`: default `omp_num_threads` is read from the environment at import time

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
