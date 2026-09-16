# Earlier versions

These files are here so the history of the program can be checked rather than taken on trust.
**None of them is maintained. Run `RaGCAn.py` in the root of this repository.**

| File | md5 | Date | What it is |
|---|---|---|---|
| `v1_PY-GEMINI_2021.py` | `36c14d927d13b24f8e2e7a1e38428189` | 2021-08-02 | The original program, 396 lines, written by hand. It answered one question about a set of genomes sequenced in the author's own laboratory |
| `v2_pygemini_full.py` | `4f2cb772384f533548faef3178e7dc1f` | 2026 | The full pipeline, which aligns the ortholog groups and infers a tree. Preserved because it is where the alignment and tree code lives. **It produced no result reported in the paper** |
| `../RaGCAn.py` | `d3cfa0d1fe2ae197686ef4cfc42f115d` | 2026 | **The version the paper describes.** Run this one |

## The one that produced the published numbers

Every number in the paper came from `pygemini.py`, md5 `9d55089baa73ca459ad0634808535616`.
`RaGCAn.py` in the root of this repository is that same file renamed. Six lines differ and all six
are display strings: the logger name, the `argparse` program name and three usage examples. No
computational line differs. `RaGCAn_rename.diff` is the diff and `RENAME_NOTE.md` explains it.

## Checking the 2021 lineage yourself

The date is the part that matters, and a copy of a file proves nothing about when it was written.
The evidence is the commit history of the original repository:

**https://github.com/ManishVictor/PY-GEMINI**, commit `3efaac0`, 2021-08-02.

`v1_PY-GEMINI_2021.py` here is byte-identical to that commit.

The paper states that textual similarity between the 2021 original and the rewritten program is
about 2% at character level while logical similarity is about 100%. Both files are in this
repository, so that claim can be checked directly:

```bash
diff versions/v1_PY-GEMINI_2021.py RaGCAn.py
```

The two figures are not in tension. A low textual and a high logical similarity is what a faithful
re-implementation looks like. The analysis parameters are byte-identical between the two.

## Why v2 is here at all

It produced nothing in the paper. It is included because the alignment-and-tree pipeline is part of
how the program developed, and leaving it out would make the lineage look tidier than it was.
