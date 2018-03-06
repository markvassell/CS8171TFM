# TFM

### Overall Goal:

---
Develop a simple prototype of fragment
assembly template-free modeling system 

###  Plan

---
* Representative protein structure database: [Bioinfo Dataset](http://bioinfo.mni.th-mh.de/pdbselect/)
* Fragment generation: Rosetta (database approach) or FRAGSIOIN (model-based approach)
* Energy function: Rosetta 3, Dfire energy function (executable available), Yang Zhang’s RW potential
(executable available), or something else
* Sampling approach
* Testing: 3 CASP12 TFM targets [CASP11](http://predictioncenter.org/casp12/index.cgi)

### Environment

---
Project interpreter: Python 3.6

Packages:  Biopython


### Helpful Links:

* [Rosetta Fragment Database Files](http://www.robetta.org/downloads/casp/casp12/fragments/)
* [Rosetta File Contents Breakdown](https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/fragment-file)
* [CASP12 Files](http://predictioncenter.org/casp12/targetlist.cgi) 
* [Biopython Help](http://biopython.org/DIST/docs/tutorial/Tutorial.html)
