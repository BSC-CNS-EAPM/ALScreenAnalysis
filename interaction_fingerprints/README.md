# Interaction fingerprints

To compute the fingerprints you need an **sdf file with all the molecules** you want to analyse,
with their docked poses, and a **pdb file with the prepared protein receptor**.

In the paper, the sdf files were the compounds obtained per target from ALScreen, selected with less strict Glide gscores.

Computing the fingerprints is costly, so it is done once by `compute_fingerprints.py`
and stored in a pkl file. `analyse_fingerprints.py` then only has to load that pkl, which means the analysis can be repeated, for example with different clustering parameters, without ever computing the fingerprints again. 

That is why there are two separate scripts. Both work on one target at a time and need the target name and the input paths to be given explicitly.

Example:

```bash
python compute_fingerprints.py --target IL1b \
    --pdb ../data/interaction_fingerprints/IL1b/receptor.pdb \
    --sdf ../data/interaction_fingerprints/IL1b/poses.sdf

python analyse_fingerprints.py --target IL1b --ifp results/IL1b_ifp.pkl
```
