# Data use, provenance, and license

## Provenance and disclosure boundary

The two data tables are a processed extraction from the Ravel-led University
of Maryland Baltimore Human Microbiome Project (UMB-HMP) longitudinal vaginal
cohort used to construct the graph in the `grip` R Journal paper. This HMP
demonstration-project cohort is distinct from the NIH HMP healthy-reference
cohort. The parent study enrolled 135 nonpregnant women of reproductive age
for daily vaginal self-sampling over ten weeks (Ravel et al., 2013, DOI:
10.1186/2049-2618-1-29). Later work identifies the cohort as UMB-HMP and its
profiles as V3--V4 16S rRNA data (Lee et al., 2023, DOI:
10.1371/journal.pcbi.1011295).

The tables contain microbial 16S rRNA feature counts, sequencing-sample
identifiers, project, platform, and sequencing-phase fields. They do not
contain clinical or diary variables, direct participant identifiers, host
sequence, or records from the related U01 cohort present in the historical
working source. The original cohort publication notes that additional
sensitive metadata are available only from the principal investigators;
those data are outside this deposit.

## License

To the extent that copyright or database rights attach to the authors'
selection and arrangement of these processed tables, they are made available
under the [Creative Commons Attribution 4.0 International license](https://creativecommons.org/licenses/by/4.0/).
The underlying UMB-HMP data remain subject to the original study's applicable
data-use conditions and normal scientific attribution practices. The
accompanying R scripts remain under the `grip` repository's GPL-3-or-later
software license.

Users should cite Ravel et al. (2013), Lee et al. (2023), and the versioned
archive record once its DOI has been assigned.
