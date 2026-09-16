# Figure scripts

Every figure in the paper is drawn by one of these. No value is typed in by hand. Each script reads
the result tables at run time and prints the numbers it used, with the file each came from, so a
figure can be checked against its source.

| Script | Draws |
|---|---|
| `make_figures.py` | Figures 1 to 3 |
| `make_figure_nature_methods.py` | Two-panel summary figure, Nature Portfolio specifications |
| `make_graphical_abstract.py` | Graphical abstract, full version |
| `make_graphical_abstract_short.py` | Graphical abstract, short version |

They need the result tables, which are not in this repository because of their size. Point the
scripts at a copy:

```bash
python3 make_figures.py --root /path/to/results
```

The directory given to `--root` must contain `LPSN_Genus_Survey/`. `RAGCAN_DATA_ROOT` works as an
environment variable instead. With neither, the script looks upward from its own location.

No generative image model was used. These are plots of measured values.
