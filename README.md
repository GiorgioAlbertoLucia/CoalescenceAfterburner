# CoalescenceAfterburner

A ROOT/C++ toy Monte Carlo for computing light (anti)nuclei coalescence yields
via the Wigner-density formalism. Nucleons (protons and neutrons) are sampled
event-by-event from an input transverse-momentum spectrum and a Gaussian
emission source, transformed into relative (Jacobi) coordinates, and combined
into a nucleus with a probability given by the overlap of their phase-space
density with the nucleus' internal wavefunction.

Currently implemented nuclei:
- **He3** (A = 3, 2p + 1n)
- **He4** (A = 4, 2p + 2n)

## How it works

For each event:
1. `PositionSampler` and `MomentumSampler` draw nucleon positions (Gaussian
   source of a configurable radius) and momenta (from an input p_T histogram,
   flat in rapidity within `±fYMax`).
2. All valid nucleon combinations (e.g. 2 protons + 1 neutron for He3) are
   formed and boosted to their center-of-mass frame (`CMFrameBooster`).
3. `JacobiTransform` converts the boosted 4-momenta/positions into relative
   Jacobi coordinates.
4. `WignerDensity` evaluates the coalescence probability from the Wigner
   density of the nucleus' internal wavefunction (Gaussian wavefunction via
   `GaussianWigner`), following the parameterization in
   [arXiv:2302.12696](https://doi.org/10.48550/arXiv.2302.12696).
5. Per-nucleus and per-nucleon kinematic distributions (p_T, rapidity, phi)
   plus the yield and its statistical uncertainty are accumulated in
   `BookKeeping` and written to the output ROOT file.

Event processing is parallelized across `nThreads` worker threads
(`CoalescenceEngine::run`), each with its own cloned histograms and RNG seed
offset, merged at the end.

## Repository layout

```
include/    Header files (engine classes, samplers, physics utilities)
src/        Implementation files
main.cxx    Entry point / run configuration
input/      Input spectra and a HEPData extraction helper script
CMakeLists.txt
```

Key classes:

| Class | Role |
|---|---|
| `CoalescenceEngine` | Abstract base: event generation, multithreaded run loop, rapidity acceptance cut |
| `CoalescenceEngineHe3` / `CoalescenceEngineHe4` | Nucleus-specific combinatorics and Wigner-density evaluation |
| `CoalescenceEngineFactory` | Builds the right engine from `Config::nucleusName` |
| `ToyMcEngine` | Driver: owns the `Config`, runs the engine, writes output |
| `PositionSampler` / `MomentumSampler` | Sample nucleon position / momentum |
| `JacobiTransform` | Converts A-body momenta/positions to relative Jacobi coordinates |
| `WignerDensity` | Abstract Wigner-density interface; `GaussianWigner` implements a Gaussian wavefunction |
| `CMFrameBooster` | Boosts a set of particles to their common CM frame |
| `BookKeeping` | Owns/fills/writes the output histograms |
| `Config` | Run configuration, loadable from a YAML file |

## Requirements

- CMake ≥ 3.14
- A C++17 compiler
- [ROOT](https://root.cern) (discovered via `root-config`, so `ROOTSYS` or a
  `root-config` on `PATH` is required)
- [yaml-cpp](https://github.com/jbeder/yaml-cpp) (path currently hardcoded in
  `CMakeLists.txt` — adjust `target_include_directories` /
  `target_link_libraries` for `coalescence_core` to your local install)

## Building

```bash
mkdir build && cd build
cmake ..
make
```

This produces the static library `libcoalescence_core.a` and the executable
`coalescence`.

## Running

The run configuration (nucleus type, source radius, input spectrum, number of
events/threads, random seed, output file) is currently set in `main.cxx` via
a `Config` struct, e.g.:

```cpp
Config config = {
    .nucleusName = "He4",
    .sourceRadius = 4.6,               // fm
    .inputPtHistogramFile = "../input/spectra_0_10.root",
    .inputPtHistogramName = "hProton_0_10",
    .outputFile = "../output/output_he4.root",
    .nEvents = 100,
    .nThreads = 20,
    .randomSeed = 42
};
```

Edit `main.cxx` to switch nucleus/parameters and rebuild, then run:

```bash
./coalescence
```

`Config` also supports `loadFromFile(path)` to load these fields from a YAML
file instead of hardcoding them in `main.cxx`.

### Output

A ROOT file containing, per run:
- Sampled nucleon and reconstructed-nucleus p_T, rapidity, and phi
  distributions
- The average coalescence yield, its uncertainty, and its per-event
  distribution

### Preparing input spectra

`input/extract_hists_from_hepdata.py` converts HEPData ROOT exports (proton
and He3 p_T spectra) into the histogram format expected as
`inputPtHistogramFile` / `inputPtHistogramName`.

## Notes / known limitations

- Nucleus type and run parameters must currently be edited directly in
  `main.cxx` before rebuilding; there is no command-line interface yet.
- `CMakeLists.txt` has a hardcoded absolute path to a local `yaml-cpp`
  install — update it for your machine before building elsewhere.