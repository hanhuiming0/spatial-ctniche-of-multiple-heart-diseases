# Spatial ctniche of multiple heart diseases

Minimal public code for reproducing manuscript figures from the integrated single-nucleus and spatial transcriptomic analyses.

## Contents

- `plot_only/`: scripts that load processed results and reproduce the main and supplementary figures.
- `config/`: portable R and Python path templates.
- `environment.yml`, `requirements.txt`, and `R_packages.txt`: software dependencies.

## Configuration

Copy the example path files and replace the placeholder roots with locations available on your system. The plotting scripts assume that the required processed tables and objects have already been generated.

## Scope

This repository intentionally excludes raw or processed expression matrices, clinical metadata, intermediate analysis objects, internal server paths, and generated figures. The spatial deconvolution used in the reported analyses was RCTD; cell2location was not used.

## Usage

Run the relevant scripts under `plot_only/figure1` through `plot_only/figure6`. Supplementary plotting code is provided under `plot_only/supplementary`.
