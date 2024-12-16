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

In the development of any practical quantum technology, from quantum computers to quantum sensors, the first step is to establish a robust theoretical framework describing the way quantum resources are generated, manipulated, and measured. In quantum optics, an important medium for both computation and communication, the primary quantum resource comes in the form of discrete packets of light called photons. While a variety of theoretical representations for manipulating photon-based quantum states exist, the most commonly-used one is the Fock state representation, where quantum states are characterised simply by how many photons they contain. However, a drawback of the Fock representation is that Fock-state computations scale exponentially with state size. Computational tools are therefore almost always required for efficient evaluation of research-size problems.

# Statement of Need

MENTAT is a Mathematica-based Computer Algebra System (CAS) designed to provide efficient, but also human-intuitive, symbolic evaluation of quantum problems expressed in the Fock representation. MENTAT was designed primarily to address a major deficit in the suite of Fock-state-based software packages available to researchers: to our knowledge, no easily-accessible symbolic computational package for performing Fock-state computation exists. Informal surveys of researchers in the field showed that many researchers felt stifled by a lack of computational software, and tools for Fock state computation were generally limited to *ad hoc* programs or in some cases hand calculation, an extremely inefficient task for problems greater than simple toy models. While some industry-produced computational packages do exist, such as Xanadu’s Python-based [Strawberry Fields](https://strawberryfields.ai/) and [MrMustard](https://github.com/XanaduAI/MrMustard) [@Killoran2019], these tools are broadly numerical-only and specialised towards a specific purpose (i.e. efficient simulation of quantum-optical circuits), and can therefore often lack the broad-scale theoretical tools needed to effectively solve research problems across a variety of disciplines.

![An example *Mathematica* 14.1 notebook using the MENTAT package. The user defines a two-mode squeezed vacuum state $|\psi\rangle$ of parameter $\chi$ up to a cutoff $n = 5$ and the Bell state $|\phi\rangle = (|0,1\rangle + |1,0\rangle)/\sqrt{2}$, finds the product state $|\xi\rangle = |\psi\rangle \otimes |\phi\rangle$, and then mixes the second mode of the two-mode squeezed state with the first Bell state mode on a beamsplitter of transmissivity $\eta$.](notebook_example.png)

In MENTAT, we develop a heuristic, human-accessible CAS package that mimics the familiar process of doing pen-and-paper Fock-state calculations as closely as possible, which affords high levels of generality while also leveraging Mathematica's computational efficiency for symbolic algebra. MENTAT allows users to freely define fully generalised quantum states, operators, and measurements by simply writing them down, without requiring any prerequisite knowledge of a specialised API or code structure; since most Fock-space operations can be constructed out of only a few primitives, MENTAT offers a large range of practical functionality in a compact package. MENTAT also offers functionality for a variety of common operators and processes found in theoretical quantum optics, such as conjugate transposition, inner and outer vector products, tensor products, etc. as well as shortcuts for common quantum operations and states such as coherent states, beamsplitters, partial tracing of states, calculation of von Neumann entropy, and more. Figure 1 demonstrates the simplicity and power of the MENTAT package by evaluating a simple optical circuit which, while conceptually simple, would otherwise be laborious to compute by hand.

MENTAT is designed primarily for use by researchers and PhD students in any application where the Fock state representation is used to analyse a quantum protocol or circuit. Research works enabled specifically by MENTAT are currently in preparation by the author, and we anticipate it to become a extremely valuable resource for future works.

# Acknowledgements

The Australian Government supported
this research through the Australian Research Council’s Linkage Projects funding scheme (Project No. LP200100601). The views expressed herein are those of the authors and are not necessarily those of the Australian Government or the Australian Research Council.

# References
