# Analogue expansion analysis

Scripts for preparing analogue libraries, inspecting their chemical space, comparing Glide docking and PELE metrics with experimental parent hits, and selecting analogues from the expansion round.

The revised scripts accept researcher-supplied paths, hit IDs, constraints, score columns, and thresholds. No target or compound is selected implicitly. Glide and PELE selection comparisons assume **lower is better**, and their selection thresholds are **strict**: equality does not pass. Molecular-similarity cutoffs used by the library helpers have separate semantics defined by `MolecularAnalysis`.

## Installation

The CSV/SDF/MolDB wrappers use these helper scripts, now supplied with the project:

- `get_sdf_from_csv.py`
- `get_moldb_from_sdf.py`
- `intersect_moldbs.py`

Keep the helpers in the script folder and pass their paths with `--script`, as shown below. If they need a separate environment, select its interpreter with `--python`. The wrappers accept helper paths as arguments, so renaming the files does not require changing the wrapper code.

All three helpers require the `MolecularAnalysis` package, specifically `MolecularAnalysis.moldb`. Install that project and its dependencies in the interpreter used to run the helpers. The MolDB format and molecular-similarity behaviour are defined by the installed `MolecularAnalysis` implementation.

## Workflow

1. Convert supplier CSVs to SDF and MolDB as needed; optionally intersect matching libraries.
2. Summarise molecular properties and inspect chemical space with t-SNE.
3. Supply best-pose Glide CSVs, PELE metrics, and manually accepted visual-inspection poses.
4. Compare PELE metrics, plot Glide/PELE relationships, and run final selection.

These scripts analyse existing outputs. They do not run docking, PELE simulations, or visual inspection, and they do not generate the union of analogue libraries.

## Generality of the original scripts

| Script | Original limitation | Revision |
| --- | --- | --- |
| `sdf_creator.py` | Fixed directory tree, target folders, helper path, and CSV column names | Input folders, helper, interpreter, columns, delimiter, and exported properties are parameters |
| `moldb_creator.py` | Fixed directory tree, target folders, and helper path | Folders, helper, and interpreter are parameters |
| `intersect_analogues.py` | Fixed folder trees, target folders, helper, and similarity | Input trees, helper, output directory, similarity, and SDF export are parameters |
| `get_metrics.py` | Fixed input directory; SMILES property calculation was unused | Directory, file pattern, properties, and output are parameters; calculation from SMILES is optional |
| `tsne_plot.py` | Fixed files, hit-ID patterns, fingerprint settings, and t-SNE settings | Explicit input and optional hit; configurable columns, fingerprint, seed, metric, perplexity, and iterations |
| `pele_analysis.py` | Already accepted a file, hit, and reference constraint, but split names at the final underscore | Explicit constraint vocabulary handles compound suffixes; configurable metric columns and saved results |
| `pele_plots.py` | Fixed target, hit, files, constraint branches, colours, and reference comparisons | Manifest-driven files and arbitrary constraint labels; exact joins; configurable metrics and reference constraint |
| `final_selection.py` | Fixed target, paths, per-hit cutoffs, fingerprint settings, and target-specific filename parsing | Manifest-driven inputs and thresholds; configurable metrics and fingerprints; exact ID/constraint matching |

The revisions generalise the supplied workflow within the documented data formats. They still require best-pose tables, consistent identifiers, valid SMILES, and numeric scores; they do not infer arbitrary schemas or automatically resolve duplicates. Finite values such as a failed-docking sentinel are not removed automatically: clean failed simulations/dockings before analysis, especially before correlations and plots.

## 1. Prepare analogue libraries

### Convert CSVs to SDF

```bash
python sdf_creator.py \
  --directories data/analogue_hunter/target_A data/analogue_hunter/target_B \
  --script ./get_sdf_from_csv.py \
  --id-column name --smiles-column smiles \
  --properties space smiles --delimiter ','
```

Processes each `*.csv` in each supplied folder, using its stem as the output prefix alongside the input. Defaults match the original supplier schema: `name`, `smiles`, and properties `space smiles`.

### Convert SDFs to MolDB

```bash
python moldb_creator.py \
  --directories data/analogue_hunter/target_A data/analogue_hunter/target_B \
  --script ./get_moldb_from_sdf.py
```

Passes all `*.sdf` files in each folder to the helper in one call, with matching output stems alongside the inputs. The helper saves pickle files first. Its current rename step searches only the process working directory, so outputs saved in the input folders may retain `.pickle` rather than become `.moldb`. See the helper behaviour notes below before running the intersection step.

### Intersect analogue sources

```bash
python intersect_analogues.py \
  --directories data/analogue_hunter/target_A data/analogue_hunter/target_B \
  --other-base data/motif_matcher \
  --script ./intersect_moldbs.py \
  --similarity 1 --output-dir results/intersections
```

For each supplied folder, looks for a folder with the same basename under `--other-base`, and pairs identical `*.moldb` filenames. Outputs use `results/intersections/<folder>/<file>_intersect` as their prefix. The helper calls `intersectMolDBs(..., simt=similarity)`, then also calls `filterSimilarity(simt=similarity)` when the value is nonzero. Thus the cutoff affects both intersection and subsequent similarity filtering. SDF export is requested by default; use `--no-sdf` to disable it. Unmatched files are reported; a folder with no matched pairs raises an error.

All three wrappers support `--dry-run`. They pass arguments without a shell, including paths containing spaces. Output files may be overwritten according to the external helper's behaviour. Folders are not searched recursively, and missing or empty input folders raise errors.

### Run the helpers directly

The batch wrappers are optional. For one input or one pair of libraries, run the helpers directly from the folder containing them:

```bash
python get_sdf_from_csv.py \
  --csv data/analogues.csv --id name --smi smiles \
  --params space smiles -d ',' -o results/analogues

python get_moldb_from_sdf.py \
  --sdfs results/analogues.sdf -o results/analogues

python intersect_moldbs.py \
  --moldb1 results/source_A.moldb --moldb2 results/source_B.moldb \
  --sim 1 --sdf -o results/analogues_intersect
```

These are separate usage examples; the intersection requires two existing libraries. Create the output directories before running a helper directly. Use a consistent `MolecularAnalysis` environment when producing and loading the serialized libraries.

| Helper | Required arguments | Optional arguments | Behaviour |
| --- | --- | --- | --- |
| `get_sdf_from_csv.py` | `--csv`, `--id`, `--smi` | `--params`, `-d` (default `,`), `--sim_filt`, `-o` | Constructs one `Mol` per CSV row, sets requested SDF properties, creates a `MolDB`, optionally filters similarity, and saves SDF |
| `get_moldb_from_sdf.py` | `--sdfs` and `-o` / `--outname`, both accepting lists of equal length | None | Loads each SDF using `MolDB(sdfDB=..., paramaters=True)`, saves each library with `saveToPickle`, then attempts to rename working-directory `.pickle` files to `.moldb` |
| `intersect_moldbs.py` | `--moldb1`, `--moldb2`, `-o` / `--outname` | `--sim` (default `1`), `--sdf` | Loads two serialized libraries, intersects and optionally filters them using the supplied similarity cutoff, prints the resulting size, optionally saves SDF, and saves a pickle only if the result is nonempty |

`--params` copies the named CSV columns into SDF properties. It can be omitted for direct CSV-to-SDF conversion; the batch wrapper requests `space` and `smiles` by default. The CSV helper checks that ID, SMILES, and requested property columns exist before processing rows.

### Behaviour of the supplied helpers

The helper code has been documented as supplied; the README update does not change it. Its inputs are configurable, with the following implementation constraints:

- **CSV output prefix:** always provide `-o`. Without it, the current CSV helper constructs a list from the input filename rather than a string output prefix.
- **CSV similarity filtering:** omit `--sim_filt` to leave it disabled. The option takes a value rather than acting as a boolean flag; any supplied value, including the text `False`, enables the call to `filterSimilarity()`. That call uses the library's default settings, and the helper does not expose a numeric similarity cutoff. The batch CSV wrapper does not expose `--sim_filt`.
- **CSV error reporting:** a row-processing exception prints an error and calls `exit()` without a nonzero status. A batch wrapper can therefore fail to recognise that conversion failed; check the expected SDF and console output.
- **Pickle renaming:** both the SDF-to-MolDB and intersection helpers rename **every** `.pickle` file in the current working directory. They do not restrict renaming to files produced by the current run and do not search the output-prefix directory. This can rename unrelated files while leaving the intended outputs with `.pickle` extensions. The `.moldb` name represents a renamed pickle, not a separate conversion step. Check generated filenames and rename only the intended outputs before using wrappers that search for `*.moldb`.
- **Intersection filtering:** `--sim 0` skips the post-intersection `filterSimilarity()` call, but still passes zero as the intersection cutoff. It does not disable intersection matching. The helper requests SDF output whenever `--sdf` is present, including for an empty result, while pickle output is saved only for nonempty results. Empty-SDF handling depends on `MolecularAnalysis`.
- **Library API compatibility:** the SDF helper passes the keyword `paramaters=True` exactly as written. The installed `MolecularAnalysis` version must accept that keyword; it has not been verified against the library implementation here.

## 2. Summarise properties and chemical space

### Molecular properties

For CSVs already containing property columns:

```bash
python get_metrics.py --directory data/analogues \
  --columns mw logp tpsa -o results/property_summary.csv
```

To calculate molecular weight, RDKit LogP, and TPSA directly from SMILES:

```bash
python get_metrics.py --directory data/analogues \
  --compute-from-smiles --smiles-column smiles \
  -o results/property_summary.csv
```

Use `--pattern` to select input files (default `*.csv`). The summary contains file, property, count, mean, maximum, minimum, and sample standard deviation (`ddof=1`). Standard deviation is undefined for fewer than two observations. Invalid/missing SMILES and invalid supplied numeric values raise errors. In calculation mode, the calculated properties are always `mw`, `logp`, and `tpsa`; `--columns` applies only to supplied properties. The output summary itself is excluded if it lies in the input directory.

### t-SNE

```bash
python tsne_plot.py --input data/union/HIT.csv \
  --smiles-column smiles --group-column space \
  --experimental-csv data/experimental_hits.csv \
  --experimental-delimiter ';' --hit HIT \
  --radius 4 --n-bits 2048 --perplexity 50 --iterations 5000 --seed 42 \
  -o results/HIT_tsne
```

The optional experimental table defaults to `ID` and `SMILES`; change these with `--experimental-id-column` and `--experimental-smiles-column`. Supply both `--experimental-csv` and `--hit`, or omit both. Hit IDs are explicit and need not follow a supplier naming pattern.

The input can also be an SDF. For SDF input, molecular structures supply the SMILES and the grouping property defaults to `space`. For CSV input, the default columns are `smiles` and `space`, with `--delimiter ','`.

Outputs are `<prefix>.pdf` and `<prefix>_embeddings.csv`. Group percentages use analogue rows only. When supplied, the reference hit participates in the embedding and is marked separately. Duplicates are retained as input observations. Invalid structures raise an error rather than becoming zero fingerprints. Perplexity must be positive and smaller than the total sample count, including the optional hit.

The default fingerprint is a 2048-bit Morgan fingerprint with radius 4. `--metric euclidean` preserves the original distance choice; `--metric jaccard` uses binary fingerprint Jaccard distance. The implementation accepts both scikit-learn iteration parameter spellings. t-SNE is exploratory: coordinates depend on inputs and settings, and separate runs do not share a common coordinate system.

## 3. Define docking and PELE inputs

`pele_plots.py` and `final_selection.py` use a CSV manifest. Copy and edit `manifest.example.csv`, with one row per target, parent hit, and constraint. Paths inside the manifest are resolved relative to the manifest's location; absolute paths are also accepted.

| Field | Meaning |
| --- | --- |
| `target` | Researcher-defined target label; also used in output paths |
| `hit` | Parent hit's exact compound ID |
| `constraint` | Researcher-defined constraint label, including labels such as `HP_HB` |
| `analogue_csv` | Best-pose analogue docking table for this row |
| `hit_csv` | Reference docking table containing exactly one row for the parent |
| `pele_csv` | PELE metrics table containing candidates and the parent |
| `visual_dir` | Required for final selection: folder of manually accepted `*.pdb` poses |
| `pose_dir` | Required when copying final PELE poses: folder containing source PDBs |
| `gscore_threshold` | Optional additional absolute Glide cutoff; blank means no extra cutoff |
| `pele_threshold` | Optional additional absolute PELE cutoff; blank means no extra cutoff |
| `pele_suffix` | Optional PELE naming suffix; defaults to `constraint` |
| `visual_suffix` | Optional accepted-pose suffix; defaults to `constraint` |
| `pose_suffix` | Optional source PELE-pose suffix; defaults to `constraint` |

By default, docking tables use `ID`, `SMILES`, and `glide_gscore`; plotting needs only `ID` and `glide_gscore`. PELE tables use `molecule`, `all_BFE`, and `lowBFE_minBE`. Column names are case-sensitive. Supply column-name options to accommodate other schemas. IDs are read as strings and preserved, including underscores and leading zeros. Align identifiers between sources yourself; commas are no longer silently replaced with underscores.

The expected PELE molecule key is `<ID>_<pele_suffix>`, and pose filenames follow `<ID>_<suffix>.pdb`. Constructing complete names avoids confusing `HP` with `HP_HB` or stripping underscores from compound IDs. Suffix overrides allow labels to differ between Glide, PELE, and pose files without target-specific code branches. PELE keys must be unique within each supplied metrics file; docking IDs must be unique within each supplied docking table. Select or aggregate poses before running these analyses.

`original_thresholds.csv` preserves every target/hit/constraint cutoff from the pasted `final_selection.py`. Copy applicable values into the manifest's `gscore_threshold` field to reproduce those settings. This reference table is not loaded automatically. Blank thresholds mean no additional cutoff; a numeric zero is a valid cutoff.

## 4. Analyse PELE results

### Compare against one chosen reference constraint

```bash
python pele_analysis.py --pele_df data/pelemetrics.csv \
  --hit HIT --const HP_HB --constraints unbias HP HB HP_HB \
  --metrics all_BFE lowBFE_minBE --output-prefix results/HIT_pele
```

Supply every constraint suffix present in the PELE file with `--constraints`. Names are parsed using the longest declared matching suffix, preserving both underscore-containing IDs and compound constraints. A missing or ambiguous reference row raises an error.

Every candidate row is compared against the parent in the explicitly selected `--const` condition. All rows belonging to that parent ID are excluded from candidate counts. The script writes:

- `<prefix>_comparisons.csv`: candidate metrics, parsed IDs/constraints, and a boolean flag per metric.
- `<prefix>_summary.csv`: passing row counts and unique-ID counts by constraint and across all constraints (`__ALL__`).

`--molecule-column` changes the PELE identifier column. `--pele-csv` and `--reference-constraint` are aliases for the original `--pele_df` and `--const` options.

### Plot Glide versus PELE

```bash
python pele_plots.py --manifest manifest.csv \
  --metrics all_BFE lowBFE_minBE \
  --reference-constraint HP_HB --output-dir results/pele_plots
```

Joins analogue docking and PELE values by exact compound/constraint keys. The parent and every reference ID listed in that row's `hit_csv` are excluded from candidate points. Docking rows without PELE metrics are reported and omitted. There is no visual-inspection or Glide-threshold filtering in this plotting step.

For each target and hit, the script writes one `<metric>_vs_glide.pdf` per metric and `joined_metrics.csv`, with boolean improvement flags. Plots use separate colours for constraints and mark the corresponding parent scores. At the output root it also writes:

- `correlations.csv`: per-constraint Pearson coefficients, p-values, and sample counts. Fewer than two points or constant values produce undefined correlations.
- `improvement_counts.csv`: per-constraint and overall counts of candidates below their matched parent metric. Optional `--reference-constraint` adds an overall comparison against that particular parent condition.
- `reference_metrics.csv`: the parent values used in each condition.

Overall counts deduplicate IDs across constraints; row counts retain repeated testing under different constraints. The optional reference constraint must exist for every target/hit included in that run. Use `--id-column`, `--score-column`, and `--pele-id-column` to override input names.

## 5. Select final analogues

```bash
python final_selection.py --manifest manifest.csv \
  --pele-metric lowBFE_minBE --radius 4 --n-bits 2048 \
  --output-dir results/final_selection --copy-poses
```

For each manifest row, candidates must satisfy **all** of the following:

1. Have a Glide score below the parent hit in the same constraint.
2. Pass any additional `gscore_threshold` supplied in that row.
3. Have an exact matching PDB filename in the visual-inspection folder.
4. Have a PELE metric below the parent in the same constraint.
5. Pass any additional `pele_threshold` supplied in that row.

The parent ID itself is excluded. Missing parent rows, missing PELE metrics for visually accepted Glide candidates, duplicate keys, invalid SMILES, and nonfinite scores raise errors. Missing source poses raise an error when `--copy-poses` is requested. Empty selections produce CSVs with headers.

Outputs under `<output-dir>/<target>/<hit>/` are:

- `<target>_<hit>_analogues.csv`: all selected constraints combined.
- `<target>_<hit>_<constraint>_analogues.csv`: selection for an individual constraint.
- `PELE_poses/<constraint>/<ID>.pdb`: copied source poses, only with `--copy-poses`.

The output root also contains `selection_summary.csv`, reporting counts before filtering and after Glide, visual-inspection, and PELE filtering.

Selection tables contain `ID`, `SMILES`, `similarity_to_hit`, `glide_gscore`, `glide_constraint`, `PELE_metric`, and `PELE_value`. `PELE_metric` records which metric was selected; `PELE_value` replaces the original fixed `PELE_lowBFE_minBE` output column. Similarity is the Tanimoto coefficient of Morgan fingerprints, calculated against the parent structure in that constraint's hit table; it is reported without applying a similarity cutoff.

Use `--id-column`, `--smiles-column`, `--score-column`, and `--pele-id-column` to change input names; `--pele-metric` selects the PELE score. Pose copying is opt-in and uses Python file operations. The whole selection is validated before output writing begins. Re-runs overwrite matching outputs; old files from earlier runs are not cleared, so use a fresh output directory when changing the manifest or selection criteria.
