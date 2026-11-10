# Changelog

## 1.1.1

- Corrected the environment commands for Windows, macOS and Linux.
- Made README links work on PyPI and added a saved-plot example.
- Linked the software, simulation protocol, core results and manuscript
  reproduction workflow, with the later supplementary analyses listed separately.
- Clarified categorical inputs and the four-gene teaching example versus the
  manuscript's 15-gene workflow.
- This documentation and packaging update leaves the statistical calculations,
  reporting defaults and bundled data unchanged from 1.1.0. The manuscript
  analyses used 1.1.0.

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
