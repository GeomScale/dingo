import pytest
import numpy as np
from dingo import PolytopeSampler
from .metrics_utils import calculate_ess

@pytest.fixture
def unit_cube():
    """Generates a 10-dimensional unit cube for benchmarking."""
    dimension = 10
    A = np.vstack([np.eye(dimension), -np.eye(dimension)])
    b = np.concatenate([np.ones(dimension), np.zeros(dimension)])
    return A, b

def test_sampling_performance_and_quality(benchmark, unit_cube):
    """
    Measures both the execution time and the statistical quality (ESS)
     of the dingo/volesti interface.
    """
    A, b = unit_cube
    sampler = PolytopeSampler()
    num_samples = 2000

    def run_benchmark():
        # 1. Sample from the polytope
        samples = sampler.sample_from_polytope(A, b, n=num_samples, method='cdhr')
        
        # 2. Convert to numpy and calculate quality metric
        samples_array = np.array(samples)
        ess_score = calculate_ess(samples_array)
        
        return samples_array, ess_score

    # The 'benchmark' fixture runs 'run_benchmark' multiple times to get stable timing
    samples, ess_score = benchmark(run_benchmark)
    
    # Assertions ensure the sampler is actually working
    assert samples.shape[0] == num_samples
    assert ess_score > 0
    
    # This will print in the 'Captured stdout' section of your pytest results
    print(f"\n[Bench Results] ESS: {ess_score:.2f} for {num_samples} samples")