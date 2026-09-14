# Ligand-local Fpocket volume calculation

`calcvol.py` recalculates the volume of the part of one or more Fpocket
pockets associated with a reference ligand. It first selects alpha spheres
located near the ligand and then estimates the **union volume** occupied by
those spheres using the same bounding-box Monte Carlo principle used by
Fpocket. Consequently, overlapping regions are counted only once.

## Required inputs

- `--fpocket-dir`: an Fpocket result directory, such as
  `PPS_Protein_out`. The script searches its `pockets` subdirectory for the
  requested `pocketN_vert.pqr` files.
- `--target`: comma-separated Fpocket pocket numbers, for example
  `2,38,51`.
- `--ligand`: a PDB file containing the reference ligand coordinates. Only
  heavy atoms are used.
- `--ligand-resname`: optional residue-name filter (for example `2OW`) when
  `--ligand` is a complete protein-ligand complex rather than a ligand-only
  PDB. Without this option, every heavy `ATOM`/`HETATM` coordinate in the input
  PDB is treated as part of the ligand.
- `--all-spheres`: bypasses the ligand-distance filter and includes every
  alpha sphere in each requested target pocket. This represents the
  no-threshold Fpocket pocket volume for the specified targets.

## Calculation method

### 1. Alpha-sphere parsing

For every requested pocket, the script reads the alpha-sphere centre

$$
\mathbf{c}_i=(x_i,y_i,z_i)
$$

and raw radius $r_i$ from `pocketN_vert.pqr`.

### 2. Ligand-proximity selection

Let $\mathbf{l}_j$ be the coordinate of ligand heavy atom $j$. Alpha
sphere $i$ is retained when its centre lies within the selected threshold
$d$ of at least one ligand heavy atom:

$$
S_{\mathrm{local}}
=
\left\{
i:\min_j\left\|\mathbf{c}_i-\mathbf{l}_j\right\|\le d
\right\}.
$$

The default is $d=3.0$ Å. This is a **centre-to-heavy-atom** criterion; the
distance is not measured from the alpha-sphere surface. Sphere selection uses
the original centre coordinates and is not affected by the radius correction
described below.

### 3. Fpocket radius correction

The current Fpocket volume implementation calls
`get_verts_volume_ptr(..., -1.6)`. To reproduce this behaviour, the radius
used for volume integration is

$$
r_i^{\mathrm{eff}}=r_i+\Delta r,
\qquad \Delta r=-1.6\ \text{Å}.
$$

The default value of `--radius-offset` is therefore `-1.6`. Using the raw PQR
radii directly would substantially overestimate the value reported by the
current Fpocket output. `--radius-offset 0` should only be used when an
uncorrected alpha-sphere union is intentionally required.

### 4. Bounding box

The smallest axis-aligned box enclosing all selected, radius-corrected spheres
is constructed as

$$
x_{\min}=\min_i(c_{i,x}-r_i^{\mathrm{eff}}),\qquad
x_{\max}=\max_i(c_{i,x}+r_i^{\mathrm{eff}}),
$$

with equivalent expressions for $y$ and $z$. Its volume is

$$
V_{\mathrm{box}}
=(x_{\max}-x_{\min})(y_{\max}-y_{\min})(z_{\max}-z_{\min}).
$$

### 5. Monte Carlo union-volume estimate

The script samples $N$ uniformly distributed random points inside the
bounding box. A sampled point $\mathbf{p}_k$ is classified as inside the
pocket when it lies inside at least one selected sphere:

$$
I_k=
\begin{cases}
1, & \exists i:\left\|\mathbf{p}_k-\mathbf{c}_i\right\|^2
\le (r_i^{\mathrm{eff}})^2,\\
0, & \text{otherwise}.
\end{cases}
$$

The alpha-sphere union volume is then estimated by

$$
\widehat{V}_{\mathrm{union}}
=V_{\mathrm{box}}\frac{N_{\mathrm{inside}}}{N},
\qquad
N_{\mathrm{inside}}=\sum_{k=1}^{N}I_k.
$$

Because a point is counted once even if it lies in several spheres, sphere
overlap is not double-counted. This differs from the simple sum

$$
\sum_i\frac{4}{3}\pi(r_i^{\mathrm{eff}})^3,
$$

which overestimates the pocket volume whenever spheres overlap. The latter is
written to `naive_sphere_volume_sum_A3` for diagnostics only and is **not** the
reported pocket volume.

### 6. Monte Carlo uncertainty

Writing $\hat{p}=N_{\mathrm{inside}}/N$, the reported sampling standard
error is

$$
\operatorname{SE}(\widehat{V})
=V_{\mathrm{box}}
\sqrt{\frac{\hat{p}(1-\hat{p})}{N}}.
$$

Increasing `--iterations` reduces this random sampling error. A fixed
`--seed` makes repeated runs reproducible. The default is 1,000,000 samples
with seed 42; per-pocket seeds are deterministically offset so that different
pockets do not reuse identical random streams.

### 7. Multiple target pockets

Each requested pocket is calculated independently. The script also produces a
`combined` result by pooling the ligand-local spheres from every requested
pocket, removing exactly duplicated spheres, and calculating their union in a
single bounding box. The combined result is therefore not the arithmetic sum
of the individual pocket volumes: spatial overlap between different Fpocket
pockets is counted only once.

If a target contains no sphere satisfying the ligand-distance threshold, its
reported volume is zero.

## Example

```powershell
python .\calcvol.py `
  --fpocket-dir "Q:\coding_dir\701_Project\Resultsbin\fpocket_res\PPS_Protein_out" `
  --target 2,38,51 `
  --ligand "Q:\path\to\ligand.pdb" `
  --distance-threshold 3.0 `
  --iterations 1000000 `
  --seed 42 `
  --output-dir ".\results" `
  --output-prefix PPS_reference_3A
```

## Outputs

The script writes matching CSV and JSON files. Both include:

- requested pocket and combined scopes;
- total and ligand-selected alpha-sphere counts;
- Monte Carlo union volume in Å³;
- sampling standard error in Å³;
- diagnostic uncorrected sum of individual effective-sphere volumes;
- bounding-box limits and volume;
- number of sampled and inside points;
- threshold, radius offset, iteration count and random seed in the JSON
  metadata.

## Interpretation and limitations

The reported value is a **ligand-associated alpha-sphere volume**, not the
physical volume occupied by the ligand and not the strict geometric
intersection between the pocket and a 3 Å shell around the ligand. Once an
alpha-sphere centre passes the 3 Å selection criterion, its complete
radius-corrected sphere contributes to the union calculation.

The estimate inherits the geometric assumptions of Fpocket and should be used
primarily for consistently processed comparisons. Differences in receptor
preparation, ligand coordinates, Fpocket parameters, distance threshold,
radius offset or Monte Carlo sampling settings can change the result.
