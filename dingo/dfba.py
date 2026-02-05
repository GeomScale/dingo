"""Dynamic flux balance analysis utilities.

This module implements a sampling-assisted dynamic FBA (dFBA) workflow inspired by
Mahadevan et al. (2002), Sánchez et al. (2014) and COMETS (Harcombe et al. 2014,
Borenstein et al. 2021). The non-uniqueness of the steady-state FBA solutions is
mitigated by sampling feasible flux distributions at every time step and picking
(at steady state) the flux profile that is closest to the one selected during the
previous step. During the first step the mean value of the sampled fluxes is
used, mirroring the suggestion in DFBAlab.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import math
import warnings

import numpy as np

from dingo.MetabolicNetwork import MetabolicNetwork

try:
    from dingo.PolytopeSampler import PolytopeSampler
except ImportError:  # pragma: no cover - optional dependency (volestipy)
    PolytopeSampler = None


@dataclass
class DynamicFBAResult:
    """Container collecting the time series produced by a dFBA simulation."""

    time: np.ndarray
    biomass: np.ndarray
    concentrations: Dict[str, np.ndarray]
    fluxes: np.ndarray


class DynamicFBA:
    """Dynamic FBA simulator that breaks ties via sampling.

    Parameters
    ----------
    model : MetabolicNetwork
        Metabolic network that will be simulated.
    time_step : float
        Time step in hours.
    total_time : float
        Total length of the simulation horizon in hours.
    initial_biomass : float
        Initial biomass concentration (gDW/L or in units compatible with the
        biomass reaction).
    initial_concentrations : dict
        Mapping from exchange reaction id to extracellular concentration
        (mmol/L or any unit consistent with the flux units).
    sample_size : int, optional
        Number of flux samples generated at every step (default 128).
    sampler_method : str, optional
        One of the sampling routines supported by
        :func:`PolytopeSampler.generate_steady_states_no_multiphase`.
    sampler_kwargs : dict, optional
        Extra keyword arguments forwarded to the sampling routine.
    biomass_tol : float, optional
        Threshold below which the biomass is considered extinct and substrate
        uptake is disabled.
    """

    def __init__(
        self,
        model: MetabolicNetwork,
        time_step: float,
        total_time: float,
        initial_biomass: float,
        initial_concentrations: Dict[str, float],
        sample_size: int = 128,
        sampler_method: str = "billiard_walk",
        sampler_kwargs: Optional[Dict[str, float]] = None,
        biomass_tol: float = 1e-12,
    ) -> None:
        if time_step <= 0:
            raise ValueError("time_step has to be strictly positive.")
        if total_time <= 0:
            raise ValueError("total_time has to be strictly positive.")
        if initial_biomass < 0:
            raise ValueError("initial_biomass must be non-negative.")
        if sample_size < 1:
            raise ValueError("sample_size must be positive.")
        if model.biomass_index is None:
            raise ValueError("The metabolic network must expose a biomass reaction.")

        self._model = model
        self._dt = float(time_step)
        self._total_time = float(total_time)
        self._initial_biomass = float(initial_biomass)
        self._sample_size = max(1, int(sample_size))
        self._sampler_method = sampler_method
        self._biomass_tol = float(biomass_tol)
        self._n_steps = int(math.ceil(self._total_time / self._dt))

        self._rxn_index = {rid: idx for idx, rid in enumerate(self._model.reactions)}
        self._tracked_exchanges = self._validate_concentrations(initial_concentrations)
        self._initial_concentrations = {
            rxn_id: float(initial_concentrations[rxn_id])
            for rxn_id in self._tracked_exchanges
        }

        self._sampler_kwargs = {
            "burn_in": 0,
            "thinning": 1,
            "variance": 1.0,
            "bias_vector": None,
            "ess": max(100, self._sample_size),
        }
        if sampler_kwargs:
            self._sampler_kwargs.update(sampler_kwargs)

        self._base_lb = np.copy(self._model.lb)

    def _validate_concentrations(self, concentrations: Dict[str, float]) -> Dict[str, int]:
        if not concentrations:
            raise ValueError("At least one extracellular concentration must be provided.")
        result = {}
        for rxn_id, conc in concentrations.items():
            if rxn_id not in self._rxn_index:
                raise KeyError(f"Reaction {rxn_id} not present in the model.")
            if conc < 0:
                raise ValueError(f"Concentration for {rxn_id} must be non-negative.")
            result[rxn_id] = self._rxn_index[rxn_id]
        return result

    def run(self) -> DynamicFBAResult:
        """Simulate the model dynamics with sampling-guided flux selection."""

        biomass = self._initial_biomass
        concentrations = {rid: val for rid, val in self._initial_concentrations.items()}
        time_trace = [0.0]
        biomass_trace = [biomass]
        conc_trace: Dict[str, list] = {
            rid: [val] for rid, val in concentrations.items()
        }
        flux_trace = []
        prev_flux: Optional[np.ndarray] = None

        try:
            for step in range(self._n_steps):
                self._apply_dynamic_bounds(concentrations, biomass)
                flux = self._select_flux(prev_flux)
                flux_trace.append(flux)
                prev_flux = flux.copy()

                mu = flux[self._model.biomass_index]
                biomass_prev = biomass
                biomass = self._update_biomass(biomass, mu)
                biomass_trace.append(biomass)

                for rxn_id, idx in self._tracked_exchanges.items():
                    concentrations[rxn_id] = self._update_concentration(
                        concentrations[rxn_id], flux[idx], biomass_prev
                    )
                    conc_trace[rxn_id].append(concentrations[rxn_id])

                time_trace.append(min(self._total_time, (step + 1) * self._dt))
        finally:
            # Restore the base bounds so that the model can be reused afterwards.
            self._model.lb = np.copy(self._base_lb)

        time_arr = np.asarray(time_trace)
        biomass_arr = np.asarray(biomass_trace)
        conc_history = {rid: np.asarray(values) for rid, values in conc_trace.items()}
        flux_arr = np.asarray(flux_trace)

        return DynamicFBAResult(time_arr, biomass_arr, conc_history, flux_arr)

    def _apply_dynamic_bounds(self, concentrations: Dict[str, float], biomass: float) -> None:
        lb = np.copy(self._base_lb)
        if biomass <= self._biomass_tol:
            for idx in self._tracked_exchanges.values():
                lb[idx] = 0.0
        else:
            for rxn_id, idx in self._tracked_exchanges.items():
                conc = concentrations[rxn_id]
                limit = -conc / (self._dt * biomass)
                lb[idx] = max(lb[idx], limit)
        self._model.lb = lb

    def _select_flux(self, prev_flux: Optional[np.ndarray]) -> np.ndarray:
        try:
            samples = self._sample_fluxes()
        except Exception as exc:  # pragma: no cover - defensive path
            warnings.warn(
                f"Sampling failed ({exc}). Falling back to a single FBA solution.",
                RuntimeWarning,
            )
            return self._fallback_flux()

        if samples.size == 0:
            return self._fallback_flux()

        if prev_flux is None:
            selected = samples.mean(axis=0)
        else:
            distances = np.linalg.norm(samples - prev_flux, axis=1)
            selected = samples[np.argmin(distances)]

        return np.asarray(selected, dtype=float)

    def _sample_fluxes(self) -> np.ndarray:
        if PolytopeSampler is None:
            raise ImportError("volestipy is required for sampling in DynamicFBA.")

        sampler = PolytopeSampler(self._model)
        steady_states = sampler.generate_steady_states_no_multiphase(
            method=self._sampler_method,
            n=self._sample_size,
            burn_in=int(self._sampler_kwargs.get("burn_in", 0)),
            thinning=int(self._sampler_kwargs.get("thinning", 1)),
            variance=float(self._sampler_kwargs.get("variance", 1.0)),
            bias_vector=self._sampler_kwargs.get("bias_vector"),
            ess=int(self._sampler_kwargs.get("ess", self._sample_size)),
        )
        return np.asarray(steady_states.T, dtype=float)

    def _fallback_flux(self) -> np.ndarray:
        flux, _ = self._model.fba()
        return np.asarray(flux, dtype=float)

    def _update_biomass(self, biomass: float, growth_rate: float) -> float:
        return max(0.0, biomass + self._dt * growth_rate * biomass)

    def _update_concentration(self, conc: float, exchange_flux: float, biomass: float) -> float:
        new_conc = conc + self._dt * exchange_flux * biomass
        return max(0.0, new_conc)
