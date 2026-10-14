# Binding Parameter Monte Carlo Simulation

R code for the Monte Carlo simulations in:

> Y. Du, *Estimation and Distribution of Parameters in Protein–Ligand Binding and Competitive Binding by Monte Carlo Simulation*.

The code simulates equilibrium binding, binding kinetics and competitive binding experiments, estimates the binding parameters by nonlinear regression, and tests the sampling distributions of the estimates with a parametric bootstrap Kolmogorov–Smirnov (KS) test. Section, table and figure numbers below refer to the manuscript and its Supplementary Materials.

## Requirements

- A Linux cluster with the SLURM job scheduler
- R 4.5.3 with the packages listed in Section 2.1 of the manuscript (minpack.lm, drc, ggplot2, plyr, RColorBrewer, extrafont, svglite, scales, parallel), plus MASS and KSgeneral for the KS tests
- `dos2unix`

The SLURM headers in the `submit-*.sh` scripts (partition, number of CPUs, time limit) are set for the author's cluster. Adjust them to your system before submitting jobs.

## Setup

```bash
cd Binding_Para_MC_Simulation

# Convert Windows line endings in scripts (needed if files were edited on Windows)
find . -type f \( -name "*.R" -o -name "*.sh" \) -exec dos2unix {} +

# Load the R environment (module names are cluster-specific)
module purge
module load apps/anaconda3/2023.09
source activate r_env
```

## Workflow

### Directory layout

| Directory | Content | Manuscript |
|---|---|---|
| `1-kd-main` | Estimation of *K*<sub>d</sub> | Section 3.1.1 |
| `1-kd-ks` | Sampling distribution of *K*<sub>d</sub>, measurement error only | Section 3.1.2; Table 1; Figure S6A, B |
| `1-kd-ks-plnoise` | Sampling distribution of *K*<sub>d</sub>, with errors in [P]<sub>0</sub> and [L]<sub>0</sub> | Section 3.1.2; Table 1; Figure S6C, D |
| `2-kin-main` | Estimation of *k*<sub>off</sub>, *k*<sub>on</sub> and kinetic *K*<sub>d</sub> | Section 3.2.1 |
| `2-kin-ks` | Sampling distributions of *k*<sub>off</sub>, *k*<sub>on</sub> and kinetic *K*<sub>d</sub> | Section 3.2.2; Table 2; Figure S9 |
| `3-comp-main` | Estimation of apparent IC<sub>50</sub>, apparent *K*<sub>i</sub> and *K*<sub>i</sub> | Section 3.3.1 |
| `3-comp-main-plinoise` | As above, with errors in [P]<sub>0</sub>, [L]<sub>0</sub> and [I]<sub>0</sub> | Section 3.3.1 |
| `3-comp-ks.ki{1,10,100}` | Sampling distributions of apparent IC<sub>50</sub> and *K*<sub>i</sub>, measurement error only | Section 3.3.2; Table S15; Figures S10, S11 |
| `3-comp-ks-plinoise.ki{1,10,100}` | As above, with errors in [P]<sub>0</sub>, [L]<sub>0</sub> and [I]<sub>0</sub> | Section 3.3.2; Table S15; Figures S12, S13 |
| `qqplot_n100000.R` | Normal Q–Q plots of the first two runs at *n* = 100,000 | Figures S6, S9–S13 |

In directory names, `plnoise` denotes random errors in [P]<sub>0</sub> and [L]<sub>0</sub>, and `plinoise` denotes random errors in [P]<sub>0</sub>, [L]<sub>0</sub> and [I]<sub>0</sub>.

### Distribution analysis steps

Every `*-ks*` directory follows the same steps:

1. **Stage 1.** `submit-stage1.sh` simulates and fits the experiments. Each job writes its results (e.g., `res_kd_expN1*_part1.rds`) to its own `job_*` subdirectory.
2. **Collect.** Move the stage-1 results into the analysis directory with `mv job_*/res* .`.
3. **Merge.** The `merge-stage1*.R` script combines the parts into one `.rda` file per parameter.
4. **Stage 2.** The `normal` scripts (e.g., `submit-stage2.normal.sh`, `submit-stage2.koff-normal.sh`) run the parametric bootstrap KS test on the raw estimates. The `lognormal` scripts (e.g., `submit-stage2.lognormal.sh`, `submit-stage2.koff-lognormal.sh`) run it on the log<sub>10</sub>-transformed estimates, which is equivalent to testing lognormality of the raw estimates.
5. **Q–Q plots.** `qqplot_n100000.R` reads a merged `.rda` file and writes `<name>_n100000_run1.png` and `<name>_n100000_run2.png` to the same directory.

`sbatch` returns as soon as a job is queued. Wait until all stage-1 jobs have finished (check with `squeue -u $USER`) before collecting and merging the results. Stage 2 and the Q–Q plots use only the merged `.rda` files, so they can be run in either order.

## Part 1. Equilibrium binding

```bash
# Estimation of Kd
cd 1-kd-main
sbatch submit-kd.sh

# Sampling distribution of Kd, measurement error only
cd ../1-kd-ks
sbatch submit-stage1.sh
# --- wait until all stage-1 jobs have finished ---
mv job_*/res* .                 # e.g., res_kd_expN1*_part1.rds
Rscript merge-stage1.R
sbatch submit-stage2.normal.sh
sbatch submit-stage2.lognormal.sh
Rscript ../qqplot_n100000.R res_kd.rda

# Sampling distribution of Kd, with errors in [P]0 and [L]0
cd ../1-kd-ks-plnoise
sbatch submit-stage1.sh
# --- wait until all stage-1 jobs have finished ---
mv job_*/res* .
Rscript merge-stage1.R
sbatch submit-stage2.normal.sh
sbatch submit-stage2.lognormal.sh
Rscript ../qqplot_n100000.R res_kd_plnoise_sd-medium_tp-medium_tl-medium.rda
```

## Part 2. Binding kinetics

```bash
# Estimation of koff, kon and kinetic Kd
cd ../2-kin-main
sbatch submit-kin.sh

# Sampling distributions of koff, kon and kinetic Kd
cd ../2-kin-ks
sbatch submit-stage1.sh
# --- wait until all stage-1 jobs have finished ---
mv job_*/res* .
Rscript merge-stage1-kin.R      # writes res_koff.rda, res_kon.rda and res_kinkd.rda
for p in koff kon kinkd; do
    sbatch submit-stage2.${p}-normal.sh
    sbatch submit-stage2.${p}-lognormal.sh
    Rscript ../qqplot_n100000.R res_${p}.rda
done
```

## Part 3. Competitive binding

```bash
# Estimation of apparent IC50, apparent Ki and Ki
cd ../3-comp-main
sbatch submit-comp.sh

# As above, with errors in [P]0, [L]0 and [I]0
cd ../3-comp-main-plinoise
sbatch submit-comp.sh
```

The six distribution directories, three *K*<sub>i</sub> values (1, 10 and 100 nM) each without and with concentration errors, follow the same steps. For example, for *K*<sub>i</sub> = 1 nM without concentration errors:

```bash
cd ../3-comp-ks.ki1
sbatch submit-stage1.sh
# --- wait until all stage-1 jobs have finished ---
mv job_*/res* .
Rscript merge-stage1-comp.R     # writes res_ic50_ki1.rda and res_ki_ki1.rda
sbatch submit-stage2.ic50-normal.sh
sbatch submit-stage2.ic50-lognormal.sh
sbatch submit-stage2.ki-normal.sh
sbatch submit-stage2.ki-lognormal.sh
Rscript ../qqplot_n100000.R res_ic50_ki1.rda
Rscript ../qqplot_n100000.R res_ki_ki1.rda
```

To run all six directories from the project root, submit stage 1 first:

```bash
dirs="3-comp-ks.ki1 3-comp-ks.ki10 3-comp-ks.ki100 \
      3-comp-ks-plinoise.ki1 3-comp-ks-plinoise.ki10 3-comp-ks-plinoise.ki100"

for d in $dirs; do
    (cd "$d" && sbatch submit-stage1.sh)
done
```

When all stage-1 jobs have finished, run the remaining steps:

```bash
for d in $dirs; do
    ki=${d##*.}                 # ki1, ki10 or ki100
    (
        cd "$d" &&
        mv job_*/res* . &&
        Rscript merge-stage1-comp.R &&
        for p in ic50 ki; do
            sbatch submit-stage2.${p}-normal.sh
            sbatch submit-stage2.${p}-lognormal.sh
            Rscript ../qqplot_n100000.R res_${p}_${ki}.rda
        done
    )
done
```

## Related software

Theoretical binding curves (Figures S5 and S8) were obtained with Binding Curve Viewer (<https://www.eplatton.net/binding-curve-viewer/>):

> Y. Du, Binding Curve Viewer: Visualizing the Equilibrium and Kinetics of Protein–Ligand Binding and Competitive Binding, *J. Chem. Inf. Model.* 64 (2024) 4180–4192. <https://doi.org/10.1021/acs.jcim.4c00130>

## Contact

Yu Du, Department of Clinical Laboratory, The Second Affiliated Hospital of Jiaxing University (duyu@zjxu.edu.cn)
