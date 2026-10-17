"""Erasure correction for circular contour averages in a Laurent RKHS.

All error norms are deterministic, conditional on the stated q-norm bound.
Floating-point evaluations are not outward-rounded interval certificates.
"""
from dataclasses import dataclass
import numpy as np


def spectrum(n, q):
    if not isinstance(n, (int, np.integer)) or n < 2 or not 0 < q < 1:
        raise ValueError("Require integer n >= 2 and 0 < q < 1")
    k = np.arange(n)
    lam = n * (q ** k + q ** (n-k)) / (-np.expm1(n*np.log(q)))
    if not np.all(np.isfinite(lam)) or min(lam) <= 0:
        raise ValueError("Spectrum underflows; use higher precision")
    return lam


def kernel(theta, q):
    """Poisson kernel, evaluated without subtraction of nearly equal numbers."""
    return (1-q*q)/((1-q)**2 + 4*q*np.sin(np.asarray(theta)/2)**2)


def error_norm(weights, q):
    """Sharp unit-ball error of ANY real grid weights, via nonnegative aliases."""
    w = np.asarray(weights, dtype=float)
    n = len(w)
    lam = spectrum(n, q)
    ft = np.fft.fft(w)
    s = w.sum()
    t = q**n
    e2 = (1-s)**2 + 2*t/(1-t)*s*s
    e2 += np.dot(lam[1:]/n, np.abs(ft[1:])**2)
    return float(np.sqrt(e2))


@dataclass
class Rule:
    weights: np.ndarray
    error: float
    loss: float
    gains: dict
    normal_residual: float
    condition: float


def recover(n, q, missing=(), condition_limit=1e12):
    """FFT plus an e-by-e solve; no dense n-by-n matrix is formed.

    gains[j] is the exact-arithmetic reduction in SQUARED unit-ball error
    from restoring node j. Return a numeric diagnostic, not an interval proof.
    """
    lam = spectrum(n, q)
    cond = float(max(lam)/min(lam))
    if cond > condition_limit:
        raise ValueError("Full-grid condition exceeds configured limit")
    raw = list(missing)
    if any(not isinstance(j, (int, np.integer)) or j < 0 or j >= n for j in raw):
        raise ValueError("Missing indices must be integers in [0,n)")
    if len(raw) != len(set(raw)):
        raise ValueError("Duplicate missing index")
    e = np.array(sorted(raw), dtype=int)
    u0 = 1/lam[0]
    w = np.full(n, u0)
    loss = 0.
    gains = {}
    if len(e):
        c = np.fft.ifft(1/lam).real
        b = c[(e[:, None]-e[None, :]) % n]
        rhs = np.full(len(e), u0)
        x = np.linalg.solve(b, rhs)
        v = np.zeros(n)
        v[e] = x
        w -= np.fft.ifft(np.fft.fft(v)/lam).real
        w[e] = 0.  # Enforce erasures exactly, including numerical roundoff.
        loss = float(rhs @ x)
        invb = np.linalg.solve(b, np.eye(len(e)))
        gains = {int(j): float(x[k]**2/invb[k,k]) for k,j in enumerate(e)}
    available = np.ones(n, bool)
    available[e] = False
    kw = np.fft.ifft(lam*np.fft.fft(w)).real
    residual = float(np.max(np.abs(kw[available]-1))) if available.any() else 0.
    return Rule(w, error_norm(w,q), loss, gains, residual, cond)


def dense_rule(n, q, missing=()):
    a = np.array([j for j in range(n) if j not in set(missing)], dtype=int)
    w = np.zeros(n)
    if len(a):
        theta = 2*np.pi*a/n
        k = kernel(theta[:,None]-theta[None,:],q)
        w[a] = np.linalg.solve(k,np.ones(len(a)))
    return w


def observed_error_bound(weights, q, norm_bound, pointwise_errors=None):
    if norm_bound < 0:
        raise ValueError("Norm bound must be nonnegative")
    bound = norm_bound*error_norm(weights,q)
    if pointwise_errors is not None:
        d = np.asarray(pointwise_errors)
        if d.shape != np.asarray(weights).shape or np.any(d<0):
            raise ValueError("Require a nonnegative error for each grid node")
        bound += float(np.dot(np.abs(weights),d))
    return bound


def greedy_retries(n,q,missing,budget):
    pending = set(missing)
    order=[]
    errors=[recover(n,q,pending).error]
    for _ in range(min(budget,len(pending))):
        r=recover(n,q,pending)
        j=max(sorted(pending), key=lambda j:r.gains[j])
        pending.remove(j)
        order.append(j)
        errors.append(recover(n,q,pending).error)
    return order, errors
