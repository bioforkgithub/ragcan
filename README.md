# RaGCAn

**Rapid Genome Coherence Analyzer.** A core-genome screen that asks one question about a named
prokaryotic genus. Do these genomes hold together?

Manish Prakash Victor, Institute of Marine Research, Bergen, Norway.

---

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

## Install

Only DIAMOND (v2.1.9) and numpy are needed.

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

## History

The original program was written in 2021 and is still online at
**https://github.com/ManishVictor/PY-GEMINI**. It was called PY-GEMINI then and it did taxonomic
analysis, core gene finding, and recombinant and non-recombinant gene detection. The method, the
parameters and the complete-linkage design are from that version.

RaGCAn v3.0.0 is the same method, rewritten for scale. `RaGCAn.py` has md5
`d3cfa0d1fe2ae197686ef4cfc42f115d`. It is identical in computation to the file that produced every
published number, which carries md5 `9d55089baa73ca459ad0634808535616`. Six lines differ and all
six are display strings (the logger name, the argparse program name and three usage examples).
`RaGCAn_rename.diff` in this repository is the proof.

If you edit `RaGCAn.py` its checksum changes and that claim no longer holds. Commit any such
change on its own, with a message saying what it was.

## Citation

`CITATION.cff` is in this repository. The manuscript is in preparation and this section will be
updated with the reference on acceptance.

## Licence

MIT. See `LICENSE`.

## AI assistance

This is disclosed because the work should be readable on its own terms.

The science is mine. The question, the method, the pipeline design, the parameters and
thresholds, the output design, the choice of tool stack and the interpretation of every result
are my own. So is the 2021 original, `PY-GEMINI.py`, which I wrote by hand.

An AI assistant (Anthropic Claude) was used for the engineering. It debugged the 2021 script,
restructured it into a command-line tool, made it portable and packaged it, wrote the survey and
analysis scripts in this repository, and did the large-scale runs under my direction.

The two versions were measured against each other. Textual similarity between `PY-GEMINI.py` and
the rewritten `pygemini.py` is about 2% at character level. Logical similarity, meaning the
algorithm design, is about 100%. The analysis parameters are byte-identical and 29 of 31
output-path names are preserved. A 2% textual and 100% logical result is what a faithful
re-implementation looks like. The assistant wrote the code and I designed the program.

The assistant was also wrong, and I corrected it. It removed my core-anchoring rule and called it
actively harmful. The rule is sound and it is back in the program. A separate review caught a
data-loss defect the assistant had introduced. Both are recorded in the full assistance record
kept with the manuscript.
