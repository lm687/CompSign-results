## Ignored folders
The data folder `data/` is restricted and therefore not available

## Environment
To set up the environment, run
```
module load miniconda3-4.5.4-gcc-5.4.0-hivczbz
conda create -n snakemake-globalDA -c conda-forge bioconda::snakemake bioconda::snakemake-minimal -c bioconda

```
to enter the environment, type

```
source activate snakemake-globalDA
...
conda deactivate
```

## Synthetic datasets
- Generation A: 20200625. There is a beta intercept of zero



