> ## ⚠ THIS PLAN IS STALE ON ONE POINT. READ THIS FIRST.
>
> **§3 mitigation 1 — "drop any proteome whose protein count is more than 1 SD below
> the genus median" — is SUPERSEDED and is deliberately NOT implemented.**
>
> It was superseded on 2026-08-25/26 by measurement, not by opinion. See
> `PROJECT_STATE.md` §4.1, §12.2 and §12.3:
>
> - That filter removes **all seven *Halopseudomonas*** from the 389-genome
>   *Pseudomonas* set. *Halopseudomonas* is a genuinely reduced-genome lineage, so a
>   protein-count filter cannot tell *small because fragmentary* from *small because
>   reduced* — it preferentially deletes exactly the lineages the screen exists to
>   find, and it destroys the headline result.
> - Filtering by assembly level is not a safe substitute either: "complete" does not
>   mean "gene-rich" (*C. taklimakanense* is a complete assembly with 2,553 proteins).
>
> **The implemented design is one run per genus, drafts INCLUDED, no pre-filter at
> all**, with assembly level and protein count recorded per genome as covariates and
> analysed afterwards. Phase 2's "QC filter" step in §2 is superseded the same way, and
> the `genomes dropped by QC` column named in §4 is not produced.
>
> If you find yourself writing a filter that removes genomes, stop and read
> `PROJECT_STATE.md` §12.3 first.
>
> Everything else in this plan stands. Phase 1 is built and running — see
> `HOW_TO_RESUME.md` and `reports/PROGRESS.md`.

# Domain-wide genus coherence survey — plan

**Question:** run RaGCAn genus by genus across every validly published bacterial genus
with enough sequenced species, and ask for each one: *does this genus form a single
coherent bin?* Then test the disagreements against GTDB.

**Status 2026-08-26:** **Phase 0 COMPLETE** — all three definitive runs landed
(*Chryseobacterium* 135 → 607 core; *Pseudomonas* 389 → 226 core, 2 bins, 7/7
*Halopseudomonas*; *Pseudomonadaceae* 432 → 205 core, 2 bins). **Phase 1 STARTED** —
pipeline being built and the 16 pilot genera (1,453 genomes) queued.

⚠ **§3 mitigation 1 has been REVISED — the original QC filter must not be implemented.**
See the boxed warning in §3 before writing any code.

---

## 1. Feasibility — measured, not estimated

From `reference/assembly_summary_bacteria.txt` (NCBI RefSeq, downloaded 2026-08-25):

- 22,170 reference/representative assemblies carrying a proper binomial
- 4,122 distinct genera

One genome per named species, which is the selection rule already used for
*Chryseobacterium* and *Pseudomonas*:

| Inclusion rule | Genera | Genomes to fetch |
|---|---:|---:|
| ≥2 species | 2,086 | 19,951 |
| ≥3 species | 1,442 | 18,663 |
| **≥4 species (recommended)** | **1,095** | **17,622** |
| ≥5 species | 887 | 16,790 |
| ≥10 species | 404 | 13,717 |

Below 4 species a coherence verdict carries almost no information — two genomes
either clear 65% AAI or they don't, and there is no internal structure to find.

**Largest genera** (these dominate the compute): *Streptomyces* 831, *Pseudomonas* 391,
*Paenibacillus* 333, *Flavobacterium* 327, *Sphingomonas* 189, *Microbacterium* 185,
*Corynebacterium* 183, *Vibrio* 179, *Nocardioides* 177, *Clostridium* 164.

### Cost

- **Storage:** ~30 GB of proteomes; DIAMOND hit files are transient and deleted per genus.
- **Download:** ~1,200 batched API calls, roughly 1–2 hours with rate limiting.
- **Compute:** cost per genus scales with (proteins)², so the total is dominated by the
  ~20 largest genera. Measured reference points on this machine (220 threads):
  74 genomes / 340 K proteins = 2 min; *Streptomyces* at ~6.2 M proteins is the worst
  case at several hours. **Estimated total: 1–2 days of background wall time.**

This is comfortably within reach — 256 cores, 1 TB RAM, 2.7 PB free.

---

## 2. Phases

### Phase 0 — finish what is already running *(in progress)*
The three definitive runs: *Chryseobacterium* 135, *Pseudomonas* 389, *Pseudomonadaceae* 432.
These also calibrate throughput for the sweep. Nothing else starts until they land.

### Phase 1 — build and pilot
Build the pipeline, then run ~30 genera chosen so the answers are already known:

- **Should be coherent:** *Chryseobacterium*, *Stutzerimonas*, *Halopseudomonas*,
  *Bordetella*, *Brucella*, *Yersinia*
- **Known to be problematic / split by GTDB:** *Pseudomonas*, *Clostridium*,
  *Bacillus*, *Lactobacillus*, *Streptococcus*, *Mycobacterium*
- **Recently split, so the answer is documented:** *Aeromonas*, *Ralstonia*,
  *Paraburkholderia* vs *Burkholderia*

If the pilot does not reproduce the known answers, the sweep is not worth running.
This is the go/no-go gate.

### Phase 2 — full sweep
1,095 genera, one at a time, resumable, background. Per genus: fetch proteomes →
run (no pre-filter — see §3) → record summary → mark done → **then** delete hit files.
Deleting the hit files before a genus is marked `done` destroys the mid-genus resume
point: pygemini reuses a finished DIAMOND search, so a killed run resumes in seconds
rather than hours, but only if `working/` survives.

### Phase 3 — the decisive test *(this is what makes it a paper)*
A survey that just reports "N genera look incoherent" is descriptive. The test is
whether **RaGCAn's disagreements with LPSN coincide with GTDB's documented
disagreements.** `reference/bac120_taxonomy.tsv` (GTDB, April 2026) gives, for every
genome, both its GTDB genus and its NCBI genus, so for each LPSN genus we can compute:

- does GTDB also split it, and into how many groups?
- do RaGCAn's bins align with GTDB's groups? (adjusted Rand index / cluster agreement)

Four outcomes per genus:

| RaGCAn | GTDB | Reading |
|---|---|---|
| coherent | coherent | agreement, no action |
| split | split | **the informative case** — independent corroboration |
| split | coherent | candidate false positive, or a real finding GTDB has not made |
| coherent | split | RaGCAn's 65% threshold is too permissive — the known limitation |

Row 4 is expected to be common. We have already quantified why: between-genus AAI
(*Pseudomonas* vs *Stutzerimonas*, up to 75.73%) can exceed within-genus AAI
(73.89%), so no single threshold separates them. The survey measures **how often**
that happens across the whole domain, which is a real and reportable result.

### Phase 4 — analysis and writing

---

## 3. Controls that decide whether the results mean anything

**Assembly completeness is the main confound.** The core genome is an intersection,
so one incomplete proteome shrinks it and looks exactly like a divergent organism
(the project notes warn about this; *C. antarcticum* is the worked example).
At this scale it would generate false "incoherent" verdicts wholesale.

> ### ⚠ REVISED 2026-08-26 — mitigation 1 below was WRONG and must not be implemented
>
> The original mitigation 1 said: *"Drop any proteome whose protein count is more than
> 1 SD below the genus median."* **Do not do this.** It was tested on 2026-08-25/26 and
> it is actively harmful. See `RESULTS_RECORD.md` §7 and `PROJECT_STATE.md` §4.1, §12.2.
>
> Applied to the 389-genome *Pseudomonas* set, that exact rule **removes all seven
> *Halopseudomonas*** — the entire headline result. They sit in the bottom 5.1% of the
> set by protein count (mean 3,633 vs median 5,327) because *Halopseudomonas* is a
> genuinely **reduced-genome lineage**, not because the assemblies are poor.
>
> A protein-count filter cannot distinguish *small because the assembly is fragmentary*
> from *small because the organism has a reduced genome* — and genome reduction tracks
> divergence. **It preferentially deletes exactly the lineages this survey exists to
> find.**
>
> Filtering by assembly *level* instead does not work either: the complete-genomes-only
> tier gave a **lower** core count (810) than the size-filtered tier (959), because
> "complete" does not mean "gene-rich" (*C. taklimakanense*, complete, 2,553 proteins).
>
> **Neither filter is safe. The corrected design is: no pre-filter at all. Run every
> genus with drafts included, record quality as a covariate, and test its effect
> afterwards across all 1,095 genera** (`PROJECT_STATE.md` §12.3).

Mitigations, all applied per genus — **as revised**:
1. **No pre-filter.** Run every named species in the genus, drafts included, exactly as
   a real user would. Quality is measured, not filtered on.
2. Record protein count, assembly level and N50 for every genome — as a **covariate**
   for later analysis, not as a basis for exclusion.
3. Flag any genus whose verdict is driven by a single genome (leave-one-out already
   reports this) — coherence that is restored by dropping one genome is a QC finding,
   not a taxonomic one. **This is now the primary quality control**, and it works
   *after* the run rather than corrupting the input to it.
4. Report assembly level composition alongside every verdict.
5. Re-run **complete-genomes-only as a check** for the ~50–150 genera called incoherent,
   not as a filter on all 1,095. A complete-only survey would cover only 31% of genera
   (338 of 1,097 have ≥4 complete genomes) and would be biased toward the well-studied.

**Fixed parameters across all genera**, so results are comparable:
`--pident 60 --bin-aai 65 --sensitivity very-sensitive`. Additionally, sweep
`--bin-aai` over 65–85 per genus — this is nearly free because the search is
checkpointed, and it turns a binary verdict into a curve (*at what threshold does
this genus fall apart?*), which is far more informative.

---

## 4. Tracking

Everything is recorded so the run is auditable and resumable.

```
LPSN_Genus_Survey/
├── PLAN.md                       this file
├── reference/                    NCBI assembly summaries, GTDB taxonomy
├── scripts/                      pipeline
├── proteomes/<Genus>/            fetched proteomes (one per species)
├── results/<Genus>/              full RaGCAn output per genus
├── logs/<Genus>.log              per-genus log
└── reports/
    ├── survey_status.tsv         one row per genus: state, timings, retries
    ├── survey_results.tsv        one row per genus: the verdict + all statistics
    ├── threshold_curves.tsv      bins vs --bin-aai, per genus
    ├── gtdb_comparison.tsv       Phase 3 output
    ├── failures.tsv              what broke and why
    └── PROGRESS.md               regenerated continuously, human-readable
```

`survey_status.tsv` is the resume point: a genus is `pending`, `fetching`, `running`,
`done`, `skipped` (fewer than 4 named species) or `failed`. Re-running the driver
picks up wherever it stopped. Nothing is recomputed that is already `done`, and a genus
left in `running` by an abrupt kill resumes from its checkpointed DIAMOND search.
There is no `failed QC` state — nothing is filtered out (§3).

`survey_results.tsv` columns: genus, species count, genomes used (= species count;
nothing is dropped), assembly-level composition and min/median/max protein count as
**quality covariates**, total proteins, core genes, bins at 65%, min/mean AAI, pairs ≥95% AAI, largest
leave-one-out core gain and which genome, threshold at which the genus first splits,
GTDB genus count, agreement score, wall time.

---

## 5. Honest assessment of scope

This is a substantially different paper from the *Chryseobacterium* manuscript —
stronger, but bigger. Two things worth being clear about before committing:

1. **The compute is the easy part.** 1–2 days. The analysis and writing for a
   1,095-genus survey is where the real time goes, and the interpretation burden is
   heavy: every genus RaGCAn calls incoherent needs a check before it can be
   asserted, because a wrong call there is a claim about someone's taxonomy.

2. **The expected headline is about the method's limits as much as about taxonomy.**
   Given the AAI overlap we measured, a fixed 65% threshold will merge genera that
   are genuinely distinct. The survey will quantify that across the domain. That is a
   good, honest result — but it is a paper about *when a fast screen can and cannot be
   trusted*, not a paper that redraws bacterial taxonomy.

**Recommended structure**, matching what was proposed: open with *Chryseobacterium*
as the fully worked, three-line-validated example (AAI → ANI → dDDH, with the audit
trail), use *Stutzerimonas*/*Halopseudomonas* as the independent replication against
a published genus proposal, then present the domain-wide survey as the scale test —
with GTDB concordance as the quantitative result rather than the raw verdict counts.
