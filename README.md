# Structural connectivity module

Compute individual and group structural covariance matrices, and the
corresponding cortical gradients, from cortical thickness TSV files.

The script is designed for left and right hemisphere thickness files such as:

```
thickness_68_lh_stats.tsv
thickness_68_rh_stats.tsv
```

## Installation

Tested with Python 3.11.

```bash
git clone https://github.com/nicfre1/structural-connectivity-module.git
cd structural-connectivity-module
python -m pip install -r requirements.txt
```

## Usage

Provide the left and right hemisphere thickness files and an output directory:

```bash
python compute_covariance.py \
    path/to/lh_thickness.tsv \
    path/to/rh_thickness.tsv \
    --out_dir path/to/output
```

A small example is included in `data/`:

```bash
python compute_covariance.py \
    data/BrainnetomeChild_thickness_lh_stats.tsv \
    data/BrainnetomeChild_thickness_rh_stats.tsv \
    --out_dir results
```

Use `-v` for more verbose output.

## Input format

One tab-separated file per hemisphere, with a `sample` column (subject
identifier) and one column per cortical region containing the thickness value.
Hemisphere files are merged on `sample`.

## What does the script do?

For the input left and right hemisphere thickness TSV files, the script:

1. Reads the left and right hemisphere TSV files.
2. Selects the `sample` column and all columns containing cortical thickness values.
3. Merges the left and right hemisphere data into a single DataFrame.
4. Computes a z-score across cortical regions for each subject.
5. Averages homologous left and right hemisphere regions into a single regional value.
6. Computes one individual region-by-region covariance matrix per subject using:

   ```python
   diff = z[:, None] - z[None, :]
   matrix = np.exp(-(diff ** 2))
   ```

7. Stores all individual covariance matrices in a single DataFrame.
8. Computes a group covariance matrix by averaging the individual covariance matrices across subjects.
9. Computes group cortical gradients from the group covariance matrix using BrainSpace diffusion map embedding.
10. Extracts the corresponding group eigenvalues (lambdas), representing the variance explained by each gradient.
11. Computes individual cortical gradients for each subject and aligns them to the group gradients using Procrustes alignment.
12. Computes the corresponding individual eigenvalues (lambdas).
13. Saves the output files and figures listed below.

## Outputs

Written to `--out_dir`:

* `*_individual_matrices.csv`: individual covariance matrices
* `*_group_matrix.csv`: group covariance matrix
* `*_group_gradients.csv`: group cortical gradients
* `*_group_lambdas.csv`: group eigenvalues (lambdas)
* `*_individual_gradients.csv`: individual cortical gradients
* `*_individual_lambdas.csv`: individual eigenvalues (lambdas)
* `figures/`: visualizations of the group covariance matrix, cortical
  gradients and eigenvalues, plus one figure per subject (`group/`,
  `individual_matrices/`, `individual_gradients/`)

## Citation

If you use this software, please cite it using the metadata in
[`CITATION.cff`](CITATION.cff) (GitHub also offers a "Cite this repository"
button on the repository page).

## License

Released under the [MIT License](LICENSE).
