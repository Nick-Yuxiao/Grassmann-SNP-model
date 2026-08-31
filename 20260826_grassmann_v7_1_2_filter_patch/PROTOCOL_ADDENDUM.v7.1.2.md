# Grassmann V7 protocol addendum v7.1.2

Status: **FROZEN BEFORE ANY NONEMPTY MAF RESULT WAS OBSERVED**  
Amends: v7.1.1 record-FILTER interpretation for the frozen chr22 source only  
Preserves: all v7.1.0 and v7.1.1 sample, model, masking, metric, GPU and claim decisions

## Observed source contract

The frozen phased HGDP+1KGP v2 chr22 BCF contains 1,093,149 records. Every
record has VCF `FILTER=.`; no record has the literal `PASS` value. The source
object name declares it filtered, but the VCF header/records do not encode a
per-record PASS decision. Applying the v7.1.1 literal `PASS` rule therefore
produced zero records. That empty result and its audit files are retained under
`site_selection/v7.1.1` as provenance.

## Amended executable rule

For this exact source object only, the record gate is:

1. source object SHA-256 must equal
   `09e50c698522d19bbc2c43a2ea2d29b290d63cf8c1d6983c5198dcb4289ff3fa`;
2. record `FILTER` must equal the unset value `.`;
3. record must be a biallelic SNP;
4. donor-pool MAF, recomputed from FORMAT/GT in all 2,496 frozen 1KGP release
   individuals, must be strictly greater than 0.01;
5. duplicate genomic positions are audited and later rejected rather than
   silently choosing one record;
6. HGDP evaluability is measured in a later frozen step before final L is signed.

This is named `SOURCE_PREFILTERED_FILTER_UNSET`, not `PASS`. It is not evidence
that upstream filters were applied identically to SNPBag. The compatibility arm
remains NOT_COMPARABLE. No MAF, missingness or length threshold may be adjusted
to approach 81,920.

The malformed source `FORMAT/PP` header warnings are recorded but do not affect
the GT-only AC/AN/AF/MAF/F_MISSING calculation.
