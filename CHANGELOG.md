# Changelog

## 1.1.0

- Changed the default reporting floor from two to six method-eligible patients.
  Explicit lower values remain possible and are flagged. This is a reporting
  policy change, not a change to the equal-r, REML–HKSJ or EB estimators.
- Preserved descriptive patient effects, exclusions and direction counts when
  inference is withheld. Added machine-readable reporting and warning fields.
- Added warnings for fewer than twenty eligible patients and for intervals
  wider than a configurable correlation-scale tolerance (default 0.2).
- Added an equal-r, two-sided noncentral-t planning MDE using an external SD;
  Fisher and EB targets do not inherit that planning model.
- Documented the simulation evidence, its result-informed timing and limits.
  No patient-count boundary is presented as a calibration guarantee.

## 1.0.0

Initial public software release with an AnnData interface, patient-level
association calculations, equal-r, REML–HKSJ and individual EB summaries,
a real-data tutorial and cross-platform tests.
