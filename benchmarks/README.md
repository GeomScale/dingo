# Dingo Benchmarking Suite

This directory contains performance tests to measure the overhead of the Python interface relative to the C++ volesti core.

## How to Run
1. Install dependencies: `pip install pytest-benchmark numpy`
2. Execute benchmarks: `pytest benchmarks/`

## Metrics Tracked
- **Wall-clock time:** Total execution time per sample.
- **Interface Overhead:** Time delta between Python entry and C++ execution.
- **ESS (Planned):** Statistical quality of samples.