# Constructions of S-Hadamard Matrices

This repository contains source code and data for my Master's thesis titled "Constructions of S-Hadamard Matrices".
Citations in my thesis for this repository include a section number to specify the data and files being referred to.
The section numbers cited are given below.
For each section, we describe where the relevant data and files can be found in the repository.

Note that this code in this repository has the following software version requirements:

- Maple 2024 minimum required
- Magma V2.26-12
- Python 3.11

## Table of Contents

- [Section 1: Relevant links and sources](#section-1-relevant-links-and-sources)
- [Section 2: Butson type S-Hadamard matrices](#section-2-butson-type-s-hadamard-matrices)
  - [Section 2.1: Catalog of S-Hadamard BH matrices](#section-21-catalog-of-s-hadamard-bh-matrices)
  - [Section 2.2: Full SLA-parametrizations of BH matrices](#section-22-full-sla-parametrizations-of-bh-matrices)

## Section 1: Relevant links and sources

This section is not referred to in the thesis, but included for quick links to relevant sources:

- Aalto BH Matrix Repository: <https://wiki.aalto.fi/spaces/Butson/pages/120482304/Matrices+up+to+monomial+equivalence>
- Catalog of CHM: <https://chaos.if.uj.edu.pl/~karol/hadamard/>
- Defect of BH (for matrices in Aalto repository): <https://chaos.if.uj.edu.pl/~wojtek/CHM_BH_0/index.html#BH-defect>
- Source code for BH equivalence algorithm: <https://www.daneflannery.com/classifying-cocyclic-butson-hadamard-matrices>

## Section 2: Butson type S-Hadamard matrices

This section describes references in "Chapter 2: Background" of the thesis, particularly, "Section 2.6: Butson Hadamard Matrices".
Note that examples of Butson Hadamard matrices used throughout this thesis are taken from the [Aalto repository](https://wiki.aalto.fi/spaces/Butson/pages/120482304/Matrices+up+to+monomial+equivalence).

### Section 2.1: Catalog of S-Hadamard BH matrices

The specific BH matrices that we use for square defect computations in Section 2.6.2 are given in the [SBH-Matrices](SBH-Matrices/) directory.
See this directory for more details.
Notably, the position of each $BH(n,q)$ matrix in our data file containing a list of BH(n,q)$ matrices is what we refer to as the index in the thesis.

### Section 2.2: Full SLA-parametrizations of BH matrices

The source code for computing defect and square defect values of non-isolated $BH(n,q)$ matrices given in tables in Section 2.6.2 are given in the [Square-Defect](Square-Defect/) directory.
See this directory for more details.

All full SLA-parametrizations found during these computations are saved to files named `sbh-n-q-parametrization.m` in [Square-Defect/output](Square-Defect/output/).
See this directory for details on reading and using these files to view the parametrizations.
