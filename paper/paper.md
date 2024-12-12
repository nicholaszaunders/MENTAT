---
title: 'MENTAT: a Mathematica-based CAS for accessible and efficient Fock-state computation'
tags:
  - Mathematica
  - quantum
  - optics
  - Fock state
  - computer algebra system
authors:
  - name: Nicholas Zaunders
    orcid: 0009-0008-5612-3903
    equal-contrib: true
    affiliation: "1"
affiliations:
 - name: Centre for Quantum Computation and Communication Technology, School of Mathematics and Physics, University of Queensland, St Lucia, Queensland 4072, Australia.
   index: 1
   ror: 00rqy9422
date: 12 December 2024
bibliography:
 - paper.bib
---

# Summary

In the development of any practical quantum technology, from quantum computers to quantum sensors, the first step is to establish a robust theoretical framework describing the way quantum resources are generated, manipulated, and measured. In quantum optics, an important medium for both computation [PsiQuantum] and communication [Micius], the primary quantum resource comes in the form of discrete packets of light called photons. While a variety of theoretical representations for manipulating photon-based quantum states exist, the most commonly-used one is the Fock state representation, where quantum states are characterised simply by how many photons they contain. However, a drawback of the Fock representation is that Fock-state computations scale exponentially with state size. Computational tools are therefore almost always required for efficient evaluation of research-size problems.

# Statement of Need

 `MENTAT` is a Mathematica-based Computer Algebra System (CAS) designed to provide efficient, but also human-intuitive, symbolic evaluation of quantum problems expressed in the Fock representation. ``MENTAT`` was designed primarily to address a major deficit in the suite of Fock-state-based software packages available to researchers: to our knowledge, no easily-accessible symbolic computational package for performing Fock-state computation exists. Informal surveys of researchers in the field showed that many researchers felt stifled by a lack of computational software, and were forced to write their own bespoke scripts or were forced to perform calculations by hand, an extremely inefficient task for problems greater than simple toy models. While some industry-produced computational packages do exist, such as Xanadu’s Python-based [Strawberry Fields](https://strawberryfields.ai/) and [MrMustard](https://github.com/XanaduAI/MrMustard) [@Killoran2019], these tools are broadly numerical-only and specialised towards a specific purpose (i.e. efficient simulation of quantum-optical circuits), and lack the theoretical tools needed to effectively solve research problems. (For example, StrawberryFields has no functionality for determining the density matrix of a given state, or for even defining custom quantum states, which limit its usefulness in a research context.)

In ``MENTAT``, we develop a heuristic, human-accessible CAS package that mimics the familiar process of doing pen-and-paper Fock-state calculations as closely as possible, which allows for unparalleled generality while also leveraging Mathematica's computational efficiency for symbolic algebra. The real power of ``MENTAT`` lies in allowing users to freely define fully generalised quantum states, operations, and measurements by simply writing them down, without any prerequisite knowledge of a specialised API or code structure; since it is possible to build most Fock-space operations out of only a few primitives, the limiting factor in what functionality `MENTAT` can offer is only the user’s imagination. `MENTAT` also offers functionality for a variety of common operators and processes found in theoretical quantum optics, such as conjugate transposition, inner and outer vector products, tensor products, etc. as well as shortcuts for common quantum operations and states such as coherent states, beamsplitters, partial tracing of states, calculation of von Neumann entropy, and more.

`MENTAT` is designed primarily for use by researchers and PhD students in any application where the Fock state representation is used to analyse a quantum protocol or circuit. Research works enabled specifically by `MENTAT` are currently in preparation by the author, and we anticipate it to become an important resource for future works.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

The Australian Government supported
this research through the Australian Research Coun-
cil’s Linkage Projects funding scheme (Project No.
LP200100601). The views expressed herein are those of
the authors and are not necessarily those of the Aus-
tralian Government or the Australian Research Council.

# References
