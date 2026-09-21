# Results

The tables behind the paper, as the program and the comparison scripts wrote them. Every number in
the paper can be traced to a file here. Tab-separated, one header line, no hidden formatting.

| File | What it holds |
|---|---|
| `survey_per_genus.tsv` | One row per genus for the 1,160-genus survey: core genes, bins at 65% AAI, minimum and mean AAI, GTDB genera, the outcome and the run time |
| `gtdb_agreement_per_genus.tsv` | The comparison with GTDB for the 1,147 genera that could be scored: bins, GTDB genera, adjusted Rand index and the verdict |
| `flagged_genera.tsv` | The genera the program divided, with the genomes it set apart |
| `unscored_flags_ncbi_check.tsv` | The 25 genomes the program set apart that GTDB does not classify, each checked against NCBI's own ANI taxonomy check |

GTDB release R232 (April 2026) is the reference throughout. The settings were `--pident 60
--bin-aai 65` with DIAMOND 2.1.8 `--very-sensitive`.
