# The program's name

The tool is **RaGCAn — Rapid Genome Coherence Analyzer**. The file is **`RaGCAn.py`**.

## History, from the session record

| When | What Manish said |
|---|---|
| 2026-08-24 17:27 | *"How about Rapid Analysis of Gene Coherence"* |
| 2026-08-24 17:30 | *"How about RaGCan"* → *"No it should be RaGCAn"* |
| **2026-08-24 17:31** | ***"Yes change it to that as well as rename the program"*** |
| 2026-08-25 10:14 | offered a third expansion, *"Rapid Genomes Comprehensive Analysis"* |
| **2026-08-29** | **settled: "RaGCAn — Rapid Genome Coherence Analyzer, RaGCAn.py"** |

⚠ The file rename was authorised on **2026-08-24** and was not carried out until **2026-08-29**,
when he pointed it out. `PROJECT_STATE.md` had recorded the name as "never explicitly agreed",
which was false — the transcript above shows otherwise. Both errors are now corrected.

## What the rename changed

**Nothing computational.** Six lines, all display strings:

| line | was | is |
|---|---|---|
| 2 | `PY-GEMINI - rapid core-genome screen` | `RaGCAn (Rapid Genome Coherence Analyzer) - ...` |
| 51 | `logging.getLogger('pygemini')` | `logging.getLogger('RaGCAn')` |
| 756 | `'PY-GEMINI %s - core genome screen'` | `'RaGCAn %s - core genome screen'` |
| 1146 | `prog='pygemini'` | `prog='RaGCAn'` |
| 1150 | usage example | usage example |
| 1152 | conda env name | conda env name |

`RaGCAn_rename.diff` in this directory is the full diff. Verify it yourself:

```bash
diff <(git show HEAD:Code/RaGCAn.py) Code/RaGCAn.py
```

## The two checksums, and which one the paper cites

| file | md5 | meaning |
|---|---|---|
| `pygemini.py` (paper 1, `../../pygemini_work/PY-GEMINI/`) | `9d55089baa73ca459ad0634808535616` | **the version that produced every published number** |
| `RaGCAn.py` (here) | `d3cfa0d1fe2ae197686ef4cfc42f115d` | same code, renamed, six display strings |

Both are stated in the manuscript Methods. The original is deliberately left untouched under
paper 1 so the checksum cited there still resolves.

`run_fast.py` was updated to `import RaGCAn as pg` and now resolves the module from its own
directory, so this folder is self-contained and can be cloned and run anywhere.
