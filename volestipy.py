# volestipy mock/stub module
# =============================================================================
# This provides a minimal mock of the volestipy C++ extension when the
# native build is not available (e.g. on Windows without Boost/lp_solve).
# The mock HPolytope uses cobra's sampling when available, or falls back
# to a simple hit-and-run sampler implemented in pure Python.
# =============================================================================

import numpy as np
from scipy.optimize import linprog


class HPolytope:
    """Mock HPolytope when volestipy C++ extension is not available."""

    def __init__(self, A, b):
        self.A = np.array(A, dtype=np.float64)
        self.b = np.array(b, dtype=np.float64)
        self.dim = A.shape[1]

    def _find_interior_point(self):
        """Find an interior point using multiple strategies."""
        m, n = self.A.shape

        # Strategy 1: Chebyshev center (max inscribed ball)
        try:
            norms = np.linalg.norm(self.A, axis=1)
            c_lp = np.zeros(n + 1)
            c_lp[-1] = -1  # maximize r => minimize -r

            A_lp = np.hstack([self.A, norms.reshape(-1, 1)])
            bounds = [(None, None)] * n + [(0, None)]

            res = linprog(c_lp, A_ub=A_lp, b_ub=self.b, bounds=bounds, method='highs')
            if res.success and res.x[-1] > 1e-12:
                return res.x[:n]
        except Exception:
            pass

        # Strategy 2: Find a feasible point via linprog
        try:
            c_feas = np.zeros(n)
            res = linprog(c_feas, A_ub=self.A, b_ub=self.b, method='highs')
            if res.success:
                x = res.x
                # Check strict feasibility
                if np.all(self.A @ x <= self.b + 1e-10):
                    return x
        except Exception:
            pass

        # Strategy 3: Try origin or midpoint of bounds
        if np.all(self.A @ np.zeros(n) <= self.b + 1e-10):
            return np.zeros(n)

        return None

    def _hit_and_run_step(self, x):
        """One step of coordinate hit-and-run from point x."""
        n = self.dim
        # Random direction (coordinate direction for CDHR)
        coord = np.random.randint(n)
        d = np.zeros(n)
        d[coord] = 1.0

        # Compute step range: A*(x + t*d) <= b => t*(A*d) <= b - A*x
        Ad = self.A @ d
        residual = self.b - self.A @ x

        t_min = -np.inf
        t_max = np.inf

        for i in range(len(Ad)):
            if Ad[i] > 1e-12:
                t_max = min(t_max, residual[i] / Ad[i])
            elif Ad[i] < -1e-12:
                t_min = max(t_min, residual[i] / Ad[i])

        if t_min >= t_max:
            return x  # infeasible direction, stay

        t = np.random.uniform(t_min, t_max)
        return x + t * d

    def generate_samples(self, method, n, burn_in, thinning, variance, bias_vector, solver=None):
        """Generate samples using pure-Python CDHR (coordinate hit-and-run).
        Returns shape (n_samples, dim) to match real volestipy API.
        """
        x = self._find_interior_point()
        if x is None:
            raise RuntimeError("Could not find interior point for polytope")

        # Burn-in
        for _ in range(max(burn_in, 100)):
            x = self._hit_and_run_step(x)

        # Sample — real volestipy returns (n_samples, dim)
        samples = np.zeros((n, self.dim))
        for i in range(n):
            for _ in range(max(thinning, 1)):
                x = self._hit_and_run_step(x)
            samples[i, :] = x

        return samples

    def mmcs(self, ess, psrf, parallel_mmcs, num_threads, solver=None):
        """Mock MMCS: falls back to hit-and-run with identity rounding."""
        n_samples = max(ess * 2, 1000)
        x = self._find_interior_point()
        if x is None:
            raise RuntimeError("Could not find interior point for polytope")

        # Burn-in
        for _ in range(200):
            x = self._hit_and_run_step(x)

        # Sample
        samples = np.zeros((self.dim, n_samples))
        for i in range(n_samples):
            for _ in range(2):
                x = self._hit_and_run_step(x)
            samples[:, i] = x

        Tr = np.eye(self.dim)
        Tr_shift = np.zeros(self.dim)

        return self.A, self.b, Tr, Tr_shift, samples

    def rounding(self, method="john_position", solver=None):
        """Mock rounding: returns identity transformation."""
        Tr = np.eye(self.dim)
        Tr_shift = np.zeros(self.dim)
        return self.A, self.b, Tr, Tr_shift, 1.0
