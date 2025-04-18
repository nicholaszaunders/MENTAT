# MENTAT: a *Mathematica*-based human-accessible computer algebra system for Fock state computation

## Overview

MENTAT is a computer algebra system for *Mathematica* 14.1, designed with the aim of closely mimicking the familiar process of doing pen-and-paper Fock-state calculations while leveraging Mathematica's computational efficiency for symbolic algebra calculation. 

MENTAT allows users to freely define fully generalised quantum states, operations, and measurements by simply writing them down, without any prerequisite knowledge of a specialised API or code structure. Since it is possible to build most Fock-space operations out of only a few primitives, the limiting factor in what functionality MENTAT can offer is only the user’s imagination. MENTAT also offers functionality for a variety of common operators and processes found in theoretical quantum optics, such as conjugate transposition, inner and outer vector products, tensor products, etc. as well as shortcuts for common quantum operations and states such as coherent states, beamsplitters, partial tracing of states, calculation of von Neumann entropy, and more.

## Installation and testing

MENTAT is provided as a pre-built *Mathematica* 14.1 package file. To use MENTAT in a Mathematica notebook, simply run the following line at the beginning of the document:

```
<< "path\\MENTAT.wl"
```

where `path` is the directory location of the folder in which the MENTAT.wl file is stored, e.g. `C:\\Users\\user\\...`, etc. The package has been imported successfully when running the above command shows the welcome message. You can also check that the package is working correctly by referencing the included [test notebook and examples](https://github.com/nicholaszaunders/MENTAT/tree/main/examples).

MENTAT requires Mathematica version 14.1.

## Features
Features in MENTAT include:

 - Definition of quantum states via native ket, bra, and density matrix representations
 - General inner and outer product
 - Tensor product
 - Creation and annihilation operator
 - Full and partial traces
 - Conjugate transposition
 - CV state generation (coherent state, two-mode squeezed vacuum state) with cutoff and minimum norm options
 - Conversion between Fock representation and vector representation
 - von Neumann entropy

A basic guide to MENTAT's syntax can be found in [help.md](https://github.com/nicholaszaunders/MENTAT/blob/main/examples/help.md).

## Contribution

If you have a suggestion for functionality, or find a bug or issue, please feel free to raise an issue directly on Github. For general inquiries, feel free to email at **n.zaunders@uq.edu.au**.
If you have used MENTAT for your research, please cite it using the DOI [![DOI](https://zenodo.org/badge/885098771.svg)](https://doi.org/10.5281/zenodo.15171089).
