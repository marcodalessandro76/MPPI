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
- Tests: `pytest` suite under `tests/` (to be created). Tests that need `pw.x`/`yambo` must be marked and skipped
  automatically when the executables are not in PATH, so the same suite runs on both machines.

## Conventions
- Match the existing style: classes that inherit from `dict`, Sphinx-style docstrings with `:py:class:` types,
  `verbose` flags.
- Keep the public API backward compatible (the user's notebooks and scripts on the clusters depend on it). Flag it
  explicitly when a change breaks it.

## Status: bug-fix checklist
Found in the analysis of 2026-09-30. Tick an item once it is fixed AND tested. "verified" = reproduced by running
the code; the other items come from reading it.

Runnable on the laptop:
- [ ] `PhInput.store` calls `self._slicefile`, but the method is `slicefile` → parsing any file fails (verified)
- [ ] `Dataset.fetch_results` matches ids by substring: `{'ecut':4}` also selects `ecut_40` (verified). It also
      re-runs `post_processing()` for every index (O(n²))
- [ ] `Utils.dict_merge` uses `collections.Mapping` (removed in py3.10) (verified)
- [ ] `YamboOutputParser`: Yambo writes `o-*.orbt_magnetization`, but the column dict key is `orb_magnetization` →
      columns come out as col1..col7 (verified)
- [ ] `Dos`: default `broad_kind=lorentzian` is the function, not the string → no broadening (Dos.py:121,130,155,180)
- [ ] `Dos`: `data.get_evals(set_gap,set_direct_gap)` is positional → the gap is used as a scissor (Dos.py:150,201)
- [ ] `Xn_single_frequency/Xn_frequency_mixing.from_file`: `cls(data,verbose)` passes verbose as X_order
- [ ] `Xn_frequency_mixing.py:313` uses the stale `ifreq` instead of `plot_ifreq`; debug prints at lines 312 and 319
- [ ] `Utils.file_parser` ignores its `skip` argument
- [ ] `YamboInput.read_file`: `filename` is undefined in the except branch; `exit()` on error
- [ ] `PwInput`: `exit()` on bad k-points; ATOMIC_POSITIONS with if_pos flags and `K_POINTS gamma` are not parsed;
      files are not closed
- [ ] `PwParser`: only catches TypeError; with lsda the XML has nbnd_up/nbnd_dw instead of nbnd; nbands_full is
      meaningless for metals
- [ ] `YamboDipolesParser.get_info` prints dip_v shape for dip_p; netCDF files are not closed (QP, DB utils, merge_qp)
- [ ] `BandStructure.py:178` passes allclose atol/rtol swapped; `Constants.high_sym_fcc['K']` is equivalent to X
- [ ] `Optics/Utils.fit_sum_frequencies`: `np.sqrt(residuals)[0]` → IndexError when lstsq returns an empty residual
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
