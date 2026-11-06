"""Inspectable results and a forest plot for a single association analysis."""
from __future__ import annotations

from copy import deepcopy
from textwrap import fill

import numpy as np
import pandas as pd


class AssociationResult:
    """Patient effects, cohort summary, and the settings used to obtain them.

    ``patients`` includes excluded patients and their reasons. ``metadata``
    describes preprocessing, the statistical target, and relevant limitations.
    Individual raw correlations have no attached confidence interval. For
    ``aggregate='eb'``, individual intervals condition on fitted population
    parameters and do not include their estimation uncertainty.
    """

    def __init__(self, patients: pd.DataFrame, summary: dict, metadata: dict):
        self.patients = patients.copy(deep=True)
        self._summary = deepcopy(summary)
        self.metadata = deepcopy(metadata)

    def summary(self) -> pd.DataFrame:
        """Return a one-row table; withheld inference is represented by NaN.

        ``status`` and ``reason`` distinguish technical failure, a user-selected
        patient threshold, and an available interval. An EB analysis has no
        cohort estimate, confidence interval or P value in this table.
        ``inference_available`` refers to cohort inference; EB uncertainty is
        conditional. Descriptive fields refer to method-eligible raw effects.
        ``warning_codes`` and ``reporting_status`` qualify available estimates;
        no status certifies calibration. MDE fields are model-based planning.
        """
        return pd.DataFrame([deepcopy(self._summary)])

    def plot_forest(self, ax=None, *, title=None, include_excluded=True):
        """Plot patient correlations and, when available, the cohort interval.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to draw on. A new figure is created if omitted.
        title : str, optional
            Override the main title. Withholding and caution notes remain.
        include_excluded : bool, default True
            Show excluded patient names with their exclusion reasons.

        Returns
        -------
        matplotlib.axes.Axes
            The axes containing the forest plot.
        """
        import matplotlib.pyplot as plt

        table = self.patients
        if not include_excluded:
            table = table.loc[table["eligible"]]
        table = table.reset_index(drop=True)
        is_eb = self._summary["method"] == "eb"
        has_summary = self._summary["status"] == "reported" and not is_eb
        total_rows = len(table) + int(has_summary)
        if ax is None:
            _, ax = plt.subplots(figsize=(7.0, max(3.0, 0.3 * total_rows + 1.5)))
        labels = []
        for position, row in table.iterrows():
            label = str(row["patient_id"])
            if not row["eligible"]:
                label += " (" + str(row["reason"]).replace("_", " ") + ")"
            labels.append(label)
            if pd.notna(row["r"]):
                ax.plot(row["r"], position, "o", color="#999999" if is_eb or not row["eligible"]
                        else "#176B91", markersize=4, zorder=3)
            if is_eb and pd.notna(row["eb_estimate"]):
                if pd.notna(row["eb_ci_low"]) and pd.notna(row["eb_ci_high"]):
                    ax.plot([row["eb_ci_low"], row["eb_ci_high"]], [position, position],
                            color="#B26419", linewidth=1.3)
                ax.plot(row["eb_estimate"], position, "D", color="#B26419", markersize=4)
        if has_summary:
            row = self._summary
            position = len(table)
            ax.plot([row["ci_low"], row["ci_high"]], [position, position], color="#252525", lw=1.8)
            ax.plot(row["estimate"], position, "D", color="#252525", markersize=6)
            labels.append("Equal-patient mean" if row["method"] == "equal" else "REML + HKSJ")
        ax.axvline(0, color="#777777", lw=0.8, ls="--", zorder=0)
        ax.set_yticks(np.arange(len(labels)), labels=labels)
        ax.set_ylim(len(labels) - 0.4, -0.8)
        ax.set_xlabel("Patient-level correlation")
        plot_title = title if title is not None else f"{self.metadata['anchor']} and {self.metadata['target_label']}"
        has_individual = is_eb and self._summary["status"] == "individual_estimates_only"
        if has_individual:
            plot_title += "\nEB diamonds; conditional intervals; grey observed effects"
        elif not has_summary:
            plot_title += "\nInference withheld: " + self._summary["reason"].replace("_", " ")
        if has_summary or has_individual:
            codes = self._summary.get("warning_codes", ())
            labels_by_code = {
                "below_default_min_patients": "expert override",
                "small_cohort_calibration_warning": "small cohort",
                "wide_interval": "wide interval",
                "model_based_power_warning": "planning power below target",
                "mde_exceeds_correlation_range": "MDE exceeds one",
            }
            cautions = [labels_by_code[code] for code in codes if code in labels_by_code]
            if cautions:
                plot_title += "\nCaution: " + "; ".join(cautions)
            elif codes:
                plot_title += "\nCaution: see result metadata for assumptions"
        ax.set_title("\n".join(fill(line, width=65) for line in plot_title.splitlines()))
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        ax.grid(axis="x", color="#dddddd", lw=0.5)
        ax.figure.tight_layout()
        return ax

    def __repr__(self):
        return (f"AssociationResult(method={self._summary['method']!r}, "
                f"patients={len(self.patients)}, "
                f"eligible={self._summary['n_patients_eligible']}, "
                f"status={self._summary['status']!r})")
