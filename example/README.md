# Worked example: five genomes, about a minute

Five genomes that NCBI files under *Pseudomonas*. RaGCAn is given no taxonomy. It returns two groups
and puts one genome on its own. That genome is *Pseudomonas abyssi*, which GTDB now places in
*Halopseudomonas*.

## Run it

```bash
RaGCAn.py -i example/proteomes -o example/result -t 2
python3 example/check_example.py
```

The second command compares your result with `expected/` and prints `PASS` or `FAIL`.

## What it costs

Measured on 2 CPU cores:

| | |
|---|---|
| Wall clock | 1 min 10 s |
| Peak memory | 376 MB |
| Input | 5 genomes, 25,697 proteins, 11 MB |

Installing the environment from `environment.yml` took 30 s with an empty package cache on a fast
connection. Expect longer on a slow one.

## What you should see

| Group | Genomes |
|---|---|
| 1 | *P. aeruginosa* PAO1, *P. fluorescens*, *P. putida*, *P. syringae* |
| 2 | *P. abyssi* on its own |

*P. abyssi* is below 65% AAI to all four of the others (62.4 to 64.7). The four are all above 71% to
each other. The core genome of all five is 95 genes. Without *P. abyssi* it is 432.

## The genomes

| File | Accession | GTDB (April 2026) |
|---|---|---|
| `Pseudomonas_aeruginosa_PAO1__GCF_000006765.1.faa` | GCF_000006765.1 | *Pseudomonas aeruginosa* |
| `Pseudomonas_fluorescens__GCF_900215245.1.faa` | GCF_900215245.1 | *Aquipseudomonas fluorescens* |
| `Pseudomonas_putida_NBRC_14164__GCF_000412675.1.faa` | GCF_000412675.1 | *Aquipseudomonas putida* |
| `Pseudomonas_syringae__GCF_018394375.1.faa` | GCF_018394375.1 | *Aquipseudomonas syringae* |
| `Pseudomonas_abyssi__GCF_003444685.1.faa` | GCF_003444685.1 | ***Halopseudomonas abyssi*** |

Protein sequences are the NCBI RefSeq annotations for each accession, unchanged.

## What this example also shows

GTDB splits the first four across two genera, *Pseudomonas* and *Aquipseudomonas*. RaGCAn keeps them
together at 65%. That is the program's known limitation: it lumps more often than it splits. A flag
from RaGCAn is strong evidence. A single group is weaker evidence.

The example also shows why the genus as a whole matters. *P. abyssi* is above 65% AAI to
*P. aeruginosa* taken alone. It separates because complete linkage asks every pair to clear the
threshold, and the broader genus supplies pairs that do not. Run the same test on a narrow slice of
a genus and a real boundary can disappear.

## Tested on

| Python | numpy | DIAMOND | Groups | Largest AAI difference |
|---|---|---|---|---|
| 3.11 | 1.24 | 2.1.8 | as above | reference |
| 3.14 | 2.5 | 2.2.6 | as above | 0.13 points |

The same input run twice on the same software gave byte-identical output.
