# Spatial ctniche of multiple heart diseases

Public figure-reproduction code for the integrated single-nucleus and spatial transcriptomic study.

## Contents

- `plot/`: formatted scripts corresponding to the main and supplementary result figures.
- `docs/FIGURE_CODE_MAP.md`: mapping between manuscript figures and plotting directories.
- `config/`: portable R and Python path templates.
- `environment.yml`, `requirements.txt`, and `R_packages.txt`: software dependencies.

## Configuration

Copy the example path files and replace the placeholder roots with locations available on your system. The plotting scripts assume that the required processed tables and objects have already been generated.

## Scope

This repository intentionally excludes raw or processed expression matrices, clinical metadata, intermediate analysis objects and generated figures. The spatial deconvolution used in the reported analyses was RCTD; cell2location was not used.

## Usage

Use `docs/FIGURE_CODE_MAP.md` to locate the scripts for each figure. Run the relevant scripts under `plot/figure1` through `plot/figure6`; supplementary plotting code is under `plot/supplementary`.
