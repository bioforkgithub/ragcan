# Results

The tables behind the paper, as the program and the comparison scripts wrote them. Every number in
the paper can be traced to a file here. Tab-separated, one header line, no hidden formatting.

| File | What it holds |
|---|---|
| `survey_per_genus.tsv` | One row per genus for the 1,160-genus survey: core genes, bins at 65% AAI, minimum and mean AAI, GTDB genera, the outcome and the run time |
| `gtdb_agreement_per_genus.tsv` | The comparison with GTDB for the 1,147 genera that could be scored: bins, GTDB genera, adjusted Rand index and the verdict |
| `flagged_genera.tsv` | The genera the program divided, with the genomes it set apart |
| `unscored_flags_ncbi_check.tsv` | The 25 genomes the program set apart that GTDB does not classify, each checked against NCBI's own ANI taxonomy check |
| `two_to_three_species_genera.tsv` | The extension: the 1,046 genera with two or three named species, one row each, and the genomes set apart that GTDB does not classify |
| `genera_tested.txt` | Every genus tested, survey and extension, in one readable list with its result |
| `enterobacteriaceae_family_screen.tsv` | The 37 genera GTDB places in Enterobacteriaceae, screened as one set of 456 genomes, with the three predictions made before the run |
| `eukaryote_sets.tsv` | Yeasts, *Caenorhabditis*, *Drosophila* and primates: each genome with its bin, the sweep against NCBI genera, and the six predictions made before the runs |
| `threshold_sweep_per_genus.tsv` | For every scored genus, how the agreement with GTDB changes as the cut-off moves from 50% to 95% |
| `RaGCAn_results.xlsx` | The same results as a workbook: one sheet per organism group, one row per database genus, the bins side by side, and a sheet linking all 20,522 genomes to NCBI |
| `make_workbook.py` | The script that built the workbook, kept as a record; it reads the survey working folders, which are too large to publish |

The survey covered every genus with four or more named species. The extension covered those with
two or three. Together that is every bacterial and archaeal genus with at least two named species
that have a usable RefSeq genome: 2,206 genera, 20,434 genomes screened.

GTDB release R232 (April 2026) is the reference throughout. The settings were `--pident 60
--bin-aai 65` with DIAMOND 2.1.8 `--very-sensitive`.
