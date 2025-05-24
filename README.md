# Cell_Annotation_Tools_Benchmark

Welcome to the **Cell_Annotation_Tools_Benchmark** repository, a comprehensive benchmark study evaluating the performance of three automated cell type annotation tools—**SingleR**, **scmap**, and **scANVI**—for single-cell RNA-sequencing (scRNA-seq) data. This project systematically assesses these tools under varying reference data availability, batch heterogeneity, and computational constraints, providing practical guidance for tool selection in real-world scRNA-seq analysis.

This repository contains the code, data processing pipeline, and visualization scripts used to generate the results presented in the associated research report (Task 82074). The study focuses on pragmatic robustness, addressing challenges such as scarce reference data, batch effects, and computational scalability.

## Table of Contents
1. [Project Overview](#project-overview)
2. [Repository Structure](#repository-structure)
3. [Dependencies](#dependencies)
4. [Installation](#installation)
5. [Usage](#usage)
6. [Data](#data)
7. [Pipeline Details](#pipeline-details)
8. [Key Findings](#key-findings)
9. [Visualization](#visualization)
10. [Contributing](#contributing)
11. [License](#license)
12. [References](#references)

## Project Overview
Single-cell RNA-sequencing (scRNA-seq) enables high-resolution analysis of cellular heterogeneity, but accurate cell type annotation remains a critical challenge. Automated annotation tools streamline this process, yet their performance varies under real-world constraints such as limited reference data, batch effects, and computational resources. This project benchmarks three representative tools:

- **SingleR (v2.4.1)**: A correlation-based method using Spearman rank correlation on highly variable genes.
- **scmap (v1.24.0, scmap-cell)**: A projection-based method employing k-nearest neighbors or centroids with cosine similarity.
- **scANVI (v1.3.1)**: A deep learning method using a variational autoencoder (VAE) for probabilistic data integration and batch correction.

The benchmark evaluates these tools across:
- **Reference Scale**: From extremely scarce (1.25% of cells) to abundant (75% of cells) reference sets.
- **Batch Heterogeneity**: Using a single CEL-Seq2 dataset and a combined 10x Chromium + CEL-Seq2 dataset to simulate batch effects.
- **Computational Efficiency**: Measuring execution time, CPU RAM, and GPU VRAM (for scANVI).

The results provide a pragmatic decision framework for selecting the most suitable tool based on data characteristics and computational constraints.

## Repository Structure
