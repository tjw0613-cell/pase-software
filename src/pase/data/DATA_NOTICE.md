# GSE234933 tutorial data

These small files contain selected, transformed observations from publicly
deposited third-party data. They are **not covered by PASE's MIT software license**,
and PASE's authors do not claim ownership of, or grant a new license to, the source
data. No patient names or inferred clinical attributes are included. The retained
specimen labels and cell barcodes are the identifiers deposited by the study.

Source: [GSE234933](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE234933),
public since 24 July 2023. Please cite:

Bill R, Wirapati P, Messemaker M, et al. CXCL9:SPP1 macrophage polarity identifies
a network of cellular programs that control human cancers. Science.
2023;381:515-524. https://doi.org/10.1126/science.ade2292

GEO places no repository restrictions on use or distribution, while noting that
submitters may hold intellectual-property rights. See the
[GEO data-use terms](https://www.ncbi.nlm.nih.gov/geo/info/disclaimer.html).
Respect any applicable source-holder terms and acknowledge the source study in
subsequent use. Inclusion here does not imply endorsement by the source authors.

`gse234933_demo_manifest.json` lists the original input checksums, extraction
rules, dimensions and output checksums. `examples/build_gse234933_subset.py`
reconstructs these files from the deposited count archive and annotation table.
The sample uses only the author's broad myeloid annotation, positive library size
and deterministic cell sampling; it is not selected for a favorable association.

`X` is log1p(CP10K), with the denominator calculated from **all original genes**
before retaining the tutorial's genes. The reduced `counts` layer cannot recover
the full library by summing its columns. Use the supplied full library sizes.

The 15-gene SPP1-positive macrophage list is recorded in the manifest and can be
traced to https://doi.org/10.3389/fimmu.2021.694801. Using it in this software
tutorial does not recreate the source publication's scoring procedure.

This is a compact teaching dataset, not a replication of the full PASE paper or
an independent external validation cohort.
