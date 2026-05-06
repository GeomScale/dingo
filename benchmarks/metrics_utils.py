import numpy as np
import arviz as az

def calculate_ess(samples):
    """
    Calculates the Effective Sample Size (ESS) using ArviZ.
    This provides a robust measure of how well the dingo/volesti 
    sampler is exploring the space.
    """
    # Convert numpy samples to ArviZ InferenceData format
    # samples shape: (num_samples, dimensions)
    idata = az.convert_to_dataset(samples[np.newaxis, :]) 
    
    # Calculate ESS for all dimensions
    ess_values = az.ess(idata)
    
    # Return the mean ESS across all dimensions as a single metric
    return float(np.mean(ess_values.x.values))

def report_convergence_quality(ess_value, total_samples):
    """Provides a qualitative rating of the sampling quality."""
    ratio = ess_value / total_samples
    if ratio > 0.1:
        return "Excellent"
    elif ratio > 0.01:
        return "Acceptable"
    else:
        return "Poor (High Autocorrelation)"