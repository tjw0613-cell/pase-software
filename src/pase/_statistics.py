"""Patient-level PASE summaries. The r and Fisher-z summaries have distinct targets."""
from __future__ import annotations

import math
from functools import lru_cache

import numpy as np
from scipy import optimize, special


def _number(value):
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        result = float(value)
    except (ValueError, TypeError, OverflowError):
        return None
    return result if math.isfinite(result) else None


def _count(value):
    result = _number(value)
    if result is None or not result.is_integer() or not 0 <= result < 2**53:
        return None
    return int(result)


def _patient(row):
    pid = row.get("patient_id")
    if pid is None or not str(pid).strip():
        raise ValueError("Each record requires a patient_id.")
    eligible = str(row.get("eligible", True)).strip().lower()
    if eligible not in ("true", "false", "1", "0"):
        raise ValueError(f"Unrecognized upstream eligibility for patient {pid}.")
    r = _number(row.get("r", row.get("partial_r")))
    n = _count(row.get("n", row.get("n_cells")))
    rank = _count(row.get("design_rank"))
    k = _count(row.get("k")) if "k" in row else rank - 1 if rank else None
    equal_reason = None
    if eligible in ("false", "0"):
        equal_reason = row.get("failure_reason") or "upstream_ineligible"
    elif r is None or abs(r) > 1:
        equal_reason = "missing_or_invalid_r"
    elif n is not None and k is not None and n - k - 1 < 2:
        equal_reason = "insufficient_residual_df"
    fisher_reason = equal_reason
    if fisher_reason is None:
        if abs(r) == 1:
            fisher_reason = "boundary_r"
        elif n is None:
            fisher_reason = "missing_or_invalid_complete_case_n"
        elif k is None:
            fisher_reason = "missing_or_invalid_covariate_rank"
        elif n - k - 3 <= 0:
            fisher_reason = "nonpositive_fisher_denominator"
        elif "design_rank" in row and "k" in row and rank != k + 1:
            fisher_reason = "inconsistent_design_rank"
    return {
        "patient_id": str(pid), "r": r, "n": n, "k": k,
        "equal_eligible": equal_reason is None, "equal_reason": equal_reason,
        "fisher_eligible": fisher_reason is None, "fisher_reason": fisher_reason,
        "z": math.atanh(r) if fisher_reason is None else None,
        "v": 1.0 / (n - k - 3) if fisher_reason is None else None,
    }


@lru_cache(maxsize=256)
def _critical(alpha, df):
    return -float(special.ndtri(alpha / 2) if df is None
                  else special.stdtrit(df, alpha / 2))


def _interval(mean, se, alpha, df=None, fisher=False):
    if se is None:
        return dict(ci_low=None, ci_high=None, se=None, p_value=None,
                    df=df, status="insufficient_patients")
    if se == 0:
        return dict(ci_low=None, ci_high=None, se=0.0, p_value=None,
                    df=df, status="degenerate_zero_se")
    margin = _critical(alpha, df) * se
    lo, hi = mean - margin, mean + margin
    if fisher:
        lo, hi = math.tanh(lo), math.tanh(hi)
    if lo >= hi:
        return dict(ci_low=None, ci_high=None, se=se, p_value=None, df=df,
                    status="degenerate_zero_width")
    t = -abs(mean / se)
    p = 2 * float(special.ndtr(t) if df is None else special.stdtr(df, t))
    return dict(ci_low=lo, ci_high=hi, se=se, p_value=p, df=df,
                status="ok")


def _equal(rows, alpha):
    r = np.array([p["r"] for p in rows if p["equal_eligible"]])
    j = len(r)
    mean = float(r.mean()) if j else None
    if j and np.all(r == r[0]):
        mean = float(r[0])
        se = 0.0 if j > 1 else None
    else:
        se = float(r.std(ddof=1) / math.sqrt(j)) if j > 1 else None
    return dict(J=j, estimate=mean, target="equal_patient_arithmetic_mean_of_finite_cell_r",
                **_interval(mean, se, alpha, df=(j - 1 if j > 1 else None)))


def _weighted_mean(z, w):
    # Center first, giving exact zero residuals for identical z values.
    return float(z[0] + np.dot(w, z - z[0]) / w.sum())


def _reml_tau2(z, v):
    """Minimize the intercept-only restricted negative log-likelihood over tau² >= 0."""
    if len(z) == 2:
        return max(0.0, float(((z[0] - z[1])**2 - v.sum()) / 2))
    if np.all(v == v[0]):
        return max(0.0, float(z.var(ddof=1) - v[0]))
    # In orthonormal intercept-free contrasts the upper bound on any squared
    # contrast coefficient is ||z-zbar||². Above this bound minus min(v),
    # every contrast's negative log-likelihood is increasing.
    ss = float(np.sum((z - z.mean())**2))
    upper = ss - float(v.min())
    if upper <= 0:
        return 0.0
    # Score stationary points can include multiple optima. Sample a log1p
    # scale anchored to the smallest variance and refine each visible basin.
    scale = float(v.min())
    limit = math.log1p(upper / scale)
    grid = np.linspace(0.0, limit, max(41, math.ceil(limit / 0.08) + 1))
    t = scale * np.expm1(grid)
    def objective(x):
        var = v + scale * math.expm1(x)
        w = 1 / var
        resid = z - _weighted_mean(z, w)
        return float(np.log(var).sum() + math.log(w.sum()) + np.dot(w, resid**2))
    var = v[None, :] + t[:, None]
    w = 1 / var
    mu = z[0] + (w * (z - z[0])).sum(axis=1) / w.sum(axis=1)
    loss = (np.log(var).sum(axis=1) + np.log(w.sum(axis=1))
            + (w * (z - mu[:, None])**2).sum(axis=1))
    best_x, best_loss = 0.0, float(loss[0])
    basins = list(np.where((loss[1:-1] <= loss[:-2]) & (loss[1:-1] <= loss[2:]))[0] + 1)
    # A shallow interior minimum can lie between the boundary and first grid point.
    basins.extend([0, len(grid) - 1])
    for i in basins:
        lo, hi = grid[max(0, i - 1)], grid[min(len(grid) - 1, i + 1)]
        fit = optimize.minimize_scalar(objective, bounds=(lo, hi), method="bounded",
                                       options={"xatol": 1e-13, "maxiter": 150})
        if not fit.success or not np.isfinite(fit.fun):
            raise ArithmeticError("restricted likelihood minimization failed")
        # Prefer the exact boundary when differences are at round-off scale.
        tol = 8 * np.finfo(float).eps * max(1.0, abs(best_loss))
        if fit.fun < best_loss - tol:
            best_x, best_loss = float(fit.x), float(fit.fun)
    return scale * math.expm1(best_x)


def _empty_weighted(j, reason):
    skipped = dict(ci_low=None, ci_high=None, p_value=None, se=None, df=None, status=reason)
    return dict(J=j, estimate=None, estimate_z=None, tau2=None, weights={},
                q=None, wald=skipped.copy(), hk=skipped.copy(), mkh=skipped.copy(),
                ci_low=None, ci_high=None, p_value=None, status=reason)


def _fit_z(pids, z, v, tau2, alpha, random):
    w = 1 / (v + tau2)
    mu = _weighted_mean(z, w)
    j = len(z)
    se = math.sqrt(1 / float(w.sum()))
    weights = dict(zip(pids, (float(x) for x in w / w.sum())))
    target = ("backtransformed_random_effects_mean_fisher_z" if random
              else "backtransformed_common_effect_or_precision_weighted_fisher_z")
    result = dict(J=j, estimate=math.tanh(mu), estimate_z=mu, weights=weights,
                  target=target, tau2=tau2, status="ok")
    wald = _interval(mu, se, alpha, fisher=True)
    if not random:
        return {**result, **wald}
    df = j - 1
    q = float(np.dot(w, (z - mu)**2) / df)
    return dict(result, q=q, wald=wald,
                hk=_interval(mu, math.sqrt(q) * se, alpha, df, True),
                mkh=_interval(mu, math.sqrt(max(1.0, q)) * se, alpha, df, True))


def summarize_patients(records, alpha=0.05, shrinkage=False):
    """Summarize one record per independent patient. See README.md for schema."""
    if not 0 < alpha < 1:
        raise ValueError("alpha must lie strictly between zero and one.")
    rows = [_patient(row) for row in records]
    if len(set(p["patient_id"] for p in rows)) != len(rows):
        raise ValueError("Duplicate patient_id. Supply one effect per independent patient.")
    fit_rows = [p for p in rows if p["fisher_eligible"]]
    j = len(fit_rows)
    ids = [p["patient_id"] for p in fit_rows]
    z = np.array([p["z"] for p in fit_rows])
    v = np.array([p["v"] for p in fit_rows])
    out = dict(equal_r=_equal(rows, alpha), patients=rows, alpha=alpha, warnings=[],
               shrinkage=[], shrinkage_status="not_requested")
    if j < len([p for p in rows if p["equal_eligible"]]):
        out["warnings"].append("Fisher analyses use a smaller eligible patient set than equal-r.")
    if j < 2:
        for name in ("fixed_z", "dl", "reml"):
            out[name] = _empty_weighted(j, "insufficient_patients")
        if shrinkage:
            out["shrinkage_status"] = "insufficient_patients"
        return out
    if j <= 5:
        out["warnings"].append("Few patients: heterogeneity and t intervals have limited precision.")
    out["fixed_z"] = _fit_z(ids, z, v, 0.0, alpha, False)
    w = 1 / v
    q_fe = float(np.dot(w, (z - out["fixed_z"]["estimate_z"])**2))
    out["fixed_z"]["Q"] = q_fe
    # 2 sum_{i<j} w_i w_j / sum(w), stable under extreme imbalance.
    sorted_w = np.sort(w)
    c_dl = 2 * float(np.dot(sorted_w[1:], np.cumsum(sorted_w)[:-1])) / w.sum()
    tau2_dl = max(0.0, float((q_fe - j + 1) / c_dl))
    out["dl"] = _fit_z(ids, z, v, tau2_dl, alpha, True)
    try:
        tau2_reml = _reml_tau2(z, v)
    except (ArithmeticError, RuntimeError, ValueError) as error:
        out["reml"] = _empty_weighted(j, "reml_optimization_failed")
        out["warnings"].append(str(error))
        if shrinkage:
            out["shrinkage_status"] = "reml_optimization_failed"
        return out
    reml = out["reml"] = _fit_z(ids, z, v, tau2_reml, alpha, True)
    if shrinkage:
        status = "conditional_plugin" if tau2_reml > 0 else "zero_tau2_complete_pooling"
        out["shrinkage_status"] = status
        for p in fit_rows:
            b = tau2_reml / (tau2_reml + p["v"])
            mean = reml["estimate_z"] + b * (p["z"] - reml["estimate_z"])
            se = math.sqrt(b * p["v"])
            delta = _critical(alpha, None) * se
            out["shrinkage"].append(dict(patient_id=p["patient_id"], status=status,
                observed_data_weight=b, estimate_z=mean, estimate=math.tanh(mean),
                conditional_se_z=se, ci_low=math.tanh(mean - delta) if se else None,
                ci_high=math.tanh(mean + delta) if se else None))
    return out


def summarize_arrays(r, n, k, alpha=0.05, shrinkage=False):
    """Array interface using positional IDs. n and k may be scalar or patient-aligned."""
    r = np.asarray(r, dtype=object)
    if r.ndim != 1:
        raise ValueError("r must be one-dimensional.")
    def patient_values(value):
        arr = np.asarray(value, dtype=object)
        if arr.ndim == 0:
            return np.full(r.shape, arr.item(), dtype=object)
        if arr.shape != r.shape:
            raise ValueError("n and k must be scalars or the same length as r.")
        return arr
    rows = (dict(patient_id=i, r=value, n=count, k=rank)
            for i, (value, count, rank) in enumerate(zip(r, patient_values(n), patient_values(k))))
    return summarize_patients(rows, alpha=alpha, shrinkage=shrinkage)
