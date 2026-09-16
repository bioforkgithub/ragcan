# RaGCAn

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22791638.svg)](https://doi.org/10.5281/zenodo.22791638)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![version](https://img.shields.io/badge/version-3.0.0-informational)

**Rapid Genome Coherence Analyzer.** Ask one question about a named prokaryotic genus in minutes,
on a laptop. Do these genomes hold together?

Given 389 *Pseudomonas* genomes and no taxonomic labels, RaGCAn pulled out exactly seven. All seven
have since been reassigned to *Halopseudomonas*. The same seven came back on two further genome
sets chosen a different way.

It needs one protein search. No tree, no alignment, no reference database to download.

## What it does

RaGCAn takes a directory of proteomes, one file per genome. It builds the core proteome of the
set, computes average amino-acid identity (AAI) between every pair of genomes, and partitions the
set by complete-linkage clustering at a chosen AAI threshold (default 65%).

One bin means the genus is coherent. More than one bin is a proposed re-categorisation. The
default is complete linkage and not single linkage, so two bins merge only when every genome in
one is at or above the threshold to every genome in the other. A chain of intermediates cannot
drag two unrelated genomes into the same group.

The aim of this program was to give a fast first check on a genus, not a final taxonomic answer.
A flag from RaGCAn is a reason to run ANI, dDDH or a phylogeny. It is not a replacement for them.

## It runs on a laptop

Measured, not estimated. One genus of 17 genomes (78,933 proteins), on **2 CPU cores**:

| | |
|---|---|
| Wall clock | **8 min 27 s** |
| Peak memory | **3.6 GB** |
| Cores used | 2 |

So an ordinary laptop with 8 GB of RAM runs a typical genus. Nothing here needs a cluster. The
1,160-genus survey used many cores because it ran 1,160 genera, not because one genus is expensive.

Measured on one AMD EPYC 7763 core pair with DIAMOND 2.1.8. Lowering `--block-size` does **not**
reduce peak memory, because the memory is in the identity matrix rather than in DIAMOND, and the
results are byte-identical either way.

Typical wall time by genus size, from the 1,160-genus survey:

| Genomes in the genus | Genera | Median |
|---|---|---|
| 3 to 9 | 736 | 6 s |
| 10 to 19 | 243 | 36 s |
| 20 to 49 | 127 | 54 s |
| 50 to 99 | 33 | 3.3 min |
| 100 to 299 | 17 | 11 min |
| 300 or more | 3 | 1.6 h |

Those are at 140 threads. Divide by your core count and you will not be far off.

## Install

Only DIAMOND (2.1.8) and numpy are needed.

```bash
micromamba env create -f environment.yml
micromamba activate ragcan
```

## Run

```bash
RaGCAn.py -i proteomes/ -o results/ -t 16
```

`-i` is a directory of protein FASTA files, one per genome. It is never modified. Results go to
`-o` (default `<input>_Result`).

The thresholds that matter:

| flag | default | what it sets |
|---|---|---|
| `--bin-aai` | 65.0 | AAI at which two bins may merge, the genus call |
| `--pident` | 85.0 | percent identity for orthology |
| `--evalue` | 1e-4 | DIAMOND e-value |
| `--min-coverage` | 0.0 | alignment coverage filter, off by default |

`--bin-aai` is a parameter and not a constant of nature. Across 1,147 genera the threshold that
best reproduces GTDB moves by 41.8 points. No single value is right everywhere. 65% was chosen
because it is within 1.5 points of the best fixed value over the whole survey. Raise it if the
genus is tight, lower it if the genus is broad, and say in your methods which value you used.

## Output

`reports/` holds the answer. `summary.txt` is the short version. `aai_matrix.tsv` is every
pairwise AAI. `proposed_bins.tsv` is the partition. `merge_trace.tsv` records every merge and the
worst-case AAI that justified it, so a bin that exists only because of one pair sitting barely
above threshold can be seen for what it is.

## What it has been tested on

We ran it across 1,160 prokaryotic genera (17,979 genomes, bacteria and archaea) in 30 h on one
machine. Against GTDB, 93.78% of genomes (15,984 of 17,045 with a GTDB assignment) fall in a bin
whose majority label matches their GTDB genus. When RaGCAn does flag a genus, GTDB agrees that the
genus should be divided in 70.0% of cases (95% CI 63.3-75.9).

The worked case is *Pseudomonas*. Out of 389 genomes it pulled out 7, and all 7 were
*Halopseudomonas*. The same seven came back on 432 genomes and again on 391 genomes selected a
different way.

Its one limitation is that it lumps. Across 15 genera examined in detail it never split a genus
that should have stayed whole, but it missed splits it could have made. Thus a flag is strong
evidence and a clean result is weaker evidence.

The same measurement carries species-level information. At a 94.10% AAI cut-off all 751
same-species pairs in the survey stayed together and 99.6% of 788,943 different-species pairs
separated (99.81% balanced accuracy). RaGCAn does not perform species delimitation and no such
claim is made here.

## How it compares

| | What it answers | Cost |
|---|---|---|
| **RaGCAn** | Does this genus hold together, and if not, where does it divide | One protein search. Minutes on a laptop |
| POCP | Do two genomes belong to one genus | Comparable answers to RaGCAn in our hands, at higher cost |
| ANI | Do two genomes belong to one species | Cheap, but it is a species criterion, not a genus one |
| dDDH | Do two genomes belong to one species | The reference standard for species. Slow, usually a web submission |
| Marker-gene phylogeny | How the organisms are related | The real answer. Slow and costly, so it is run where a problem is already suspected |

Use RaGCAn to find out where to point the expensive methods. It is the first step, not the last.

## History


The original program was written in 2021 and is still online at
**https://github.com/ManishVictor/PY-GEMINI**. It was called PY-GEMINI then and it did taxonomic
analysis, core gene finding, and recombinant and non-recombinant gene detection. The method, the
parameters and the complete-linkage design are from that version.

RaGCAn v3.0.0 is the same method, rewritten for scale. It is identical in computation to the
version that produced every published number. Six lines differ and all six are display strings
(the logger name, the argparse program name and three usage examples). `RaGCAn_rename.diff` and
`RENAME_NOTE.md` in this repository record exactly what changed.

## What is in here

```
RaGCAn.py              the program. Run this one
fast_aai.py            vectorised AAI for large genera
run_fast.py            driver that calls it
environment.yml        DIAMOND and numpy, nothing else
survey_scripts/        how the 1,160-genus survey selected, fetched and ran
analysis_scripts/      the GTDB comparison, the statistics, the internal verification
figures/               every figure in the paper, drawn from the result tables
versions/              the 2021 original and the full pipeline, for checking the lineage
SURVEY_PLAN.md         the protocol, written before the survey ran
```

Nothing in `versions/` is maintained. It is there so the history can be checked rather than taken
on trust. See `versions/README.md`.

## Citation

`CITATION.cff` is in this repository. The manuscript is in preparation and this section will be
updated with the reference on acceptance.

## Licence

MIT. See `LICENSE`.

## AI assistance

This disclosure clarifies that the work is independently developed, including the research
question, methodology, pipeline design, parameters and thresholds, output design, tool stack
selection, and interpretation of results. The original 2021 script, `PY-GEMINI.py`, serves as the
foundation. An AI assistant, Anthropic Claude, was utilized for engineering tasks such as
debugging the 2021 script, restructuring it into a command-line tool, enhancing portability,
packaging, and developing the survey and analysis scripts within this repository. Large-scale runs
were also conducted under supervision.

Comparison of the two script versions shows approximately 2% textual similarity at the character
level and nearly 100% logical similarity in algorithm design. Analysis parameters remain
identical, and 29 out of 31 output-path names are retained, demonstrating a faithful
re-implementation. The updated program represents an evolved and refined version of the initial
side project.

Errors introduced by the assistant were identified and corrected, including the reinstatement of a
core-anchoring rule previously removed in error and addressing a data-loss defect. Both issues are
documented in the comprehensive assistance record accompanying the manuscript.
