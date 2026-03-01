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
- [Section 3: Computational methods and results for S-Hadamard matrices](#section-3-computational-methods-and-results-s-hadamard-matrices)
  - [Section 3.1: Order 12 S-Hadamard approximations](#section-31-order-12-s-hadamard-approximations)
  - [Section 3.2: PSLQ-based parametrization](#section-32-pslq-based-parametrization)
  - [Section 3.3: Grobner basis method for parametrizing](#section-33-grobner-basis-method-for-parametrizing)
  - [Section 3.4: Order 12 S-Hadamard parametrized matrices](#section-34-order-12-s-hadamard-parametrized-matrices)
- [Section 4: Identifying Butson Hadamard matrices](#section-4-identifying-butson-hadamard-matrices)
- [Section 5: Analysis of parametrized matrices](#section-5-analysis-of-parametrized-matrices)

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

## Section 3: Computational methods and results S-Hadamard matrices

This section describes references in "Chapter 3: Computational Methods for Constructing S-Hadamard Matrices" of the thesis.
Note that we do not include the source code for the optimizer from "Section 3.1: Finding Approximations of S-Hadamard Matrices" as it is relatively straightforward to reproduce.
We only include the relevant approximations found by the optimizer.

### Section 3.1: Order 12 S-Hadamard approximations

The 79 order 12 S-Hadamard matrices obtained by the optimizer are stored in [order-12-optimizer-solns.m](PSLQ/order-12-optimizer-solns.m).
See the associated [README.md](PSLQ/README.md) for instructions on reading matrices in this file.

### Section 3.2: PSLQ-based parametrization

All code pertaining to Section 3.2 of the thesis is given in the [PSLQ](PSLQ/) directory.
See this directory and the associated [README.md](PSLQ/README.md) for detailed information.
In this section, we link the relevant files directly for quicker access to content.

- A Maple implementation of Algorithm 1 given in Section 3.2.1 is available in [PSLQ_Param_Subroutine.mpl](PSLQ/PSLQ_Param_Subroutine.mpl).
- A Maple implementation of the naive approach for parametrizing using PSLQ is available in [PSLQ_Try_Param_Naive.mpl](PSLQ/PSLQ_Try_Param_Naive.mpl).
- A Maple implementation of the randomized approach for parametrizing using PSLQ is available in [PSLQ_Try_Param_Rand.mpl](PSLQ/PSLQ_Try_Param_Rand.mpl).

### Section 3.3: Grobner basis method for parametrizing

All code associated with the parametrization example described in Section 3.2.2 of the thesis using Grobner bases is given in the [Grobner-Basis](Grobner-Basis/) directory.
For a quick view, the pdf of the Maple worksheet we used for this example is available at [GrobnerBasisMethod.pdf](Grobner-Basis/GrobnerBasisMethod.pdf).

### Section 3.4: Order 12 S-Hadamard parametrized matrices

The set of 3-parameter order 12 S-Hadamard matrices we obtained from the approximations obtained by the optimizer using the PSLQ-based method are available in [pslq_valid_param_mats.m](PSLQ/pslq_valid_param_mats.m).
See the associated [README.md](PSLQ/README.md) for instructions on reading matrices in this file.

## Section 4: Identifying Butson Hadamard matrices

In "Chapter 4: Algebraic Construction of a $BH(12,3)$ Matrix" of the thesis, in "Section 4.1: Identifying $BH(12,3)$ Matrices", we do not explicitly reference any code within this repository since the process is relatively straightforward.
We include the code that we used in this section in [BH-Equivalence](BH-Equivalence/), see this directory for more details.

## Section 5: Analysis of parametrized matrices

This section describes references in "Chapter 5: Analysis of Parametrized Order 12 S-Hadamard Matrices" of the thesis.
