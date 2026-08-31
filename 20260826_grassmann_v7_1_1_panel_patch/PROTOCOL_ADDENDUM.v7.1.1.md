# Grassmann V7 protocol addendum v7.1.1

Status: **FROZEN FOR T01 PANEL CONSTRUCTION**  
Amends: v7.1.0 sample-release, donor-split and SNPBag-calibration clauses only  
Preserves: all v7.1.0 architecture, masking, metric, GPU and claim-ceiling decisions

## 1. Frozen source-specific sample rule

For the phased HGDP+1KGP v2 chr22 BCF and its matching gnomAD metadata, an
eligible individual is one with `release=true` in the metadata and membership in
the indexed BCF. The `release` field is the source's high-quality, relatedness-
pruned release decision. `sample_filters.all_samples_related=true` is not an
additional exclusion when `sample_filters.release_related=false`: it can mark the
retained representative of a relationship whose counterpart was excluded.

The audited source contains 2,496 eligible 1KGP individuals and 768 eligible HGDP
individuals. All HGDP individuals remain forbidden as generator donors or
preprocessing-fit data.

## 2. Donor train/validation split

The 1KGP release set is split before synthesis with seed 71001. Within every 1KGP
population, individuals are ranked by SHA-256 of
`seed + TAB + population + TAB + sample_id`. The first 10% (nearest integer,
minimum one and maximum `n-1`) form donor validation; the remainder form donor
train. Population stratification also preserves nested superpopulation coverage.
No donor may occur in both partitions. Synthetic train and validation corpora
must be generated exclusively from their corresponding donor partition.

## 3. HGDP primary holdout

With exact SNPBag calibration membership unavailable for this source release, all
768 eligible HGDP individuals form the primary external holdout. The primary
estimator remains the equal-population estimator frozen in v7.1.0. Population
block bootstrap and individual-weighted sensitivity reporting are unchanged.

## 4. SNPBag calibration availability

The publication reports 216 HGDP individuals sampled as four individuals from
each of 54 populations, from a 929-individual HGDP source. The audited joint v2
BCF contains 925 HGDP individuals and 52 population labels; after its frozen
release filter only 44 populations contain at least four eligible individuals.
Consequently this source cannot construct the exact published subset.

The calibration status for this release is `NOT_COMPARABLE_SOURCE_MISMATCH`.
The calibration sample manifest must be empty except for its header. A 180-person
or 176-person balanced subset must not be substituted or described as an exact
reproduction. If exact published sample IDs and a matching source release later
become available, activating them requires a new signed protocol and manifest
before calibration outcomes are inspected.

This availability decision does not alter the 8M matched-parameter scientific
comparison, which uses the frozen 768-person HGDP primary holdout.
