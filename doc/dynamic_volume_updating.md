# Corrected Summer Project Report

## Dynamic Volume Updating of Convex Polytopes via Random-Walk Sampling

## 1. Correct project scope

The project concerns a sequence of nested polytopes formed by adding linear
constraints. In a metabolic model the initial feasible body is

`P = {v : S v = 0, lb <= v <= ub}`.

No biomass objective floor or FBA-optimality restriction is inserted. Boundary
sampling is a separate project and has no role here.

The numerical foundation is the supplied static `relative_volume.m` method:

`rho = vol(P intersect {c^T v <= c0}) / vol(P)`.

That method builds a simulated-annealing schedule of parallel halfspaces,
estimates a telescoping product of volume ratios, and obtains the necessary
samples with PolytopeSamplerMatlab CRHMC.

This integration does not replace that logic. For a stream of cuts it invokes
the same static routine repeatedly, adding each accepted cut to the parent
before the next call.

## 2. What was wrong in the previous Dingo wrapper

The earlier folder contained three incompatible ideas:

1. `volume_reference.py` manually projected `S v = 0` into a SciPy nullspace
   and divided two independent volesti absolute-volume estimates. This bypassed
   the supplied method completely. It also did not robustly cover fixed-width reactions,
   infinite bounds, additional affine equalities, or dimension-changing cuts.
2. `MatlabVolumeEngine` attempted to import `dingo.volume_updating`, but that
   module was absent from the delivered bundle.
3. The engine expected diagnostics such as `cut_reused_pilot` and
   `cut_fresh_samples` that the supplied MATLAB code does not return.

Consequently, the wrapper was neither self-contained nor a faithful wrapper of
the static method.

## 3. Corrected Dingo architecture

The corrected integration contains two Python modules:

- `dingo.dynamic_volume`: the real implementation;
- `dingo.volume_updating`: a backward-compatible import path.

The public objects are:

- `HalfspaceCut`;
- `cuts_from_bounds`;
- `StaticVolumeUpdater`;
- backward-compatible alias `MatlabVolumeUpdater`;
- `StaticVolumeResult`.

The wrapper validates dimensions, nested bounds, options, executable paths and
all MATLAB outputs. It supports Linux MATLAB and Windows MATLAB invoked from
WSL. Temporary request/result files are isolated for each call.

## 4. Sequential static algorithm

Let `P^(0)` be the supplied parent and let cut `j` be

`C_j = {v : c_j^T v <= t_j}`.

For each `j`:

1. call `relative_volume(P^(j-1), c_j, t_j)`;
2. obtain `rho_j = vol(P^(j-1) intersect C_j)/vol(P^(j-1))`;
3. set `P^j = P^(j-1) intersect C_j`;
4. continue with the next cut.

The final estimate is

`rho = product_j rho_j`,

computed in log space. This is only the sequential composition of the static
volume estimator.

## 5. Relationship to the supplied cutting-planes ZIP

The installed MATLAB directory contains the eight original source files.

Byte-for-byte unchanged:

- `default_params.m`;
- `new_window.m`;
- `run_relative_volume.m`;
- `sliding_window.m`;
- `ul_test.m`;
- `update_window.m`.

Integration-only modification:

- `init.m` retains its compilation/setup behavior but additionally searches
  `DINGO_POLYTOPE_SAMPLER_MATLAB` and common Dingo-relative locations.
- `relative_volume.m` retains the supplied estimator and adds the two interface
  initializations normally performed by PolytopeSamplerMatlab `sample.m`:
  the shared `startTime` timeout clock and the seeded `simdTwister` stream.

New bridge:

- `volume_update_batch.m` reads a Python-generated MAT request, invokes
  `relative_volume` once per cut, appends the cut to the next parent, strips
  nonportable MATLAB handle objects from saved diagnostics, and writes a v7 MAT
  result readable by SciPy.

No annealing, ratio-estimation, convergence, CRHMC-transition, or sequential
update mathematics in `relative_volume.m` was rewritten.

The ZIP's model data, standalone random-cut runner documentation, HTML/Markdown
exposition, and Git metadata are not installed into the Dingo Python package
because they are not runtime dependencies. The supplied TeX theory is
preserved under `docs/volume_updating/theory.tex`, the original MIT license is
preserved under `docs/volume_updating/MetabolicPolytopesCuttingPlanes-LICENSE`,
and exact source hashes and both complete integration patches are under `docs/`.

## 6. Bound and GPR handling

The volume wrapper accepts linear cuts. A biological application may first use
a GPR rule and severity model to produce child reaction bounds. The helper
`cuts_from_bounds` converts every nested bound change into the correct scalar
halfspace:

- upper bound `v_k <= u'_k` becomes normal `+e_k` and threshold `u'_k`;
- lower bound `v_k >= l'_k` becomes normal `-e_k` and threshold `-l'_k`.

Both sides are emitted for a tightened reversible reaction. Any attempted bound
relaxation is rejected because the static telescoping update requires a nested
child.

GPR semantics, severity sweeps, angles and RL outputs remain application
layers. They do not change the volume estimator.

## 7. Parameter policy

The wrapper accepts only parameters present in the supplied static API:

`r`, `delta`, `alpha`, `nu`, `N_utest`, `window_size`, `epsilon`,
`max_annealing_iter`, `max_bisect_iter`, `simdLen`, `t_max_safety`, and `verb`.

Wrapper-only options are `algorithm` and `wrapper_timeout_sec`.
`algorithm` must be `static`, `baseline`, or `relative_volume`. A requested
`reuse` algorithm raises an error rather than silently running a different
method. Unsupported fields such as `ratio_ess_min` are rejected.

This means parameter changes are explicit while the underlying static logic
remains the supplied logic.

When no override is provided, the wrapper uses the ZIP's exact defaults:
`r=0.10`, `delta=0.05`, `alpha=0.20`, `nu=10`, `N_utest=320`,
`window_size=5000`, `epsilon=0.10`, `max_annealing_iter=200`,
`max_bisect_iter=100`, `simdLen=8`, `t_max_safety=1.01`, and `verb=1`.
The prior `alpha=0.10`, `N_utest=125`, `epsilon=0.05` profile is therefore not
silently imposed by this wrapper; it must be requested as a named experiment.

## 8. Result contract

The wrapper returns:

- cumulative `rho` and `log_rho`;
- per-cut ratios and log ratios;
- number of annealing phases per cut;
- total CRHMC step counts per cut;
- per-cut convergence flags;
- cut labels;
- `algorithm='static_relative_volume_sequential'`;
- `objective_floor=False`.

It does not invent fresh-sample, reused-pilot or cross-update-state statistics
that the static implementation does not expose.

## 9. Testing completed

Nine Python wrapper tests pass:

- two-cut sequential result parsing;
- identity behavior for an empty cut stream;
- bound-to-halfspace conversion including reversible two-sided tightening;
- rejection of bound relaxation;
- rejection of dimension mismatch;
- rejection of a fake `reuse` mode;
- rejection of unsupported options;
- rejection of malformed source numeric options.
- verification that the required sampler clock and random stream are initialized.

The Python modules compile successfully. A real sequential integration test
was executed through WSL using Windows MATLAB and PolytopeSamplerMatlab. For
the exact two-cut toy problem, the estimator returned rho=0.0721726401 versus
11/150=0.0733333333, corresponding to 1.58% relative error. Both cuts
converged, the phase counts were [1, 2], the CRHMC step counts were
[2488, 4072], and objective_floor=False. Total wall-clock time was 15.25 s.

## 10. Correct experimental comparison

The three scientific columns should remain distinct:

1. objective-free Dingo/volesti reference;
2. supplied static CRHMC, freshly applied to each sequential cut;
3. any future optimized reuse implementation, if built, tested against column
   2 and clearly separated from the static wrapper.

The wrapper delivered here implements column 2 only. This is intentional and
matches the requested architecture.

## 11. Deliverable

The `dynamic_volume_updating.zip` archive can be unzipped directly into the
repository. It places the wrapper files into the existing Python package and
includes the report, tests, doctor, theory, source hashes, and exact integration
diff.
