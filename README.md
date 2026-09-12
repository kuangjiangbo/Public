# IL-6 and keloid risk — Mendelian randomization analysis code

Analysis code for the manuscript:

> **Genetic Evidence Linking Interleukin-6 Levels to Keloid Risk: A Mendelian Randomization Study**
> Yong-fei Wang, Jiang-bo Kuang

## Public data sources (GWAS summary statistics)

| Trait | Accession | Source | N |
|---|---|---|---|
| IL-6 levels | GCST90012005 | Folkersen et al., *Nat Metab* 2020 | 21,758 |
| IL-1RA levels | GCST90012004 | Folkersen et al., *Nat Metab* 2020 | 21,758 |
| TNFR1 levels | GCST90012015 | Folkersen et al., *Nat Metab* 2020 | 21,758 |
| CRP levels | GCST90029070 | Said et al., *Nat Commun* 2022 | 575,531 |
| Keloid | GCST90018874 | Sakaue et al., *Nat Genet* 2021 (European-ancestry stratum) | 481,912 (668 cases) |

All datasets are freely available from the **GWAS Catalog** (<https://www.ebi.ac.uk/gwas/>)
and mirrored in **IEU OpenGWAS** (<https://gwas.mrcieu.ac.uk/>).
VCF inputs are not redistributed here — download by accession.

## Software

- R 4.4.3
- TwoSampleMR 0.6.29, MRPRESSO 1.0, coloc 5.2.3, VariantAnnotation 1.52.0, ieugwasr, biomaRt
- Python 3 (numpy / pandas) for instrument extraction

## Contents

- `01_univariate_MR/` — instrument selection, harmonisation, IVW plus five sensitivity
  methods, Cochran's Q, MR-Egger intercept, Steiger, leave-one-out.
- `02_MVMR/` — instrument extraction and multivariable MR (IL-6 ± CRP / IL-1RA / TNFR1),
  conditional F-statistics, Q-statistics, forest plot.
- `03_MR-PRESSO/` — MR-PRESSO global and outlier tests.
- `04_colocalization/` — Bayesian colocalisation (coloc 5.2.3) at the IL6R and HLA loci,
  with prior and MAF sensitivity analyses.

## Reproduction notes

- All coordinates are GRCh37/hg19.
- The 23 IL-6 instruments are confined to the IL6R (chr1) and HLA (chr6) regions; this is
  discussed in the manuscript as a key methodological limitation.
- Outputs quoted in the manuscript (tables, conditional F-statistics, colocalisation
  posteriors) are written by these scripts to `.csv`/`.txt` in the same directories.

## Licence

Code released under the MIT licence. Data remain subject to the terms of the original
GWAS providers.
