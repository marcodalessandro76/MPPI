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
- Tests: `python -m pytest tests` (fixtures in `tests/conftest.py`). Tests that need the executables are marked
  `@pytest.mark.requires_qe` / `@pytest.mark.requires_yambo` and are skipped when `pw.x`/`yambo` are not in PATH, so
  the same suite runs on both machines. `Reference_data/nl_results` (~1.3 GB) is not in git: the tests that use it
  are skipped when it is missing. Every bug fix gets a test.

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
