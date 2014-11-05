# Running mCADRE
James Eddy, Price Lab @ Institute for Systems Biology  
Last updated **Nov 5, 2014**

## Overview

The **mCADRE** method operates on a global reconstruction for the organism of interest and attempts to remove invidual reactions to obtain a representative tissue-specific model. Expression-based and literature-based evidence inputs, along with the connectivity of reactions within the global model, are used to rank reactions. The reactions with lowest supporting evidence for the tissue/organism are thereby given highest priority for removal.

All steps in the method, from ranking reactions to pruning the generic model, are performed within the main `mcadre` function, which can be called as follows:

```
[GM, C, NC, PM, Z, model_C, pruneTime, cRes] = …
    mcadre(model, G, U, confidenceScores, salvageCheck, C_H_genes, method);
```

## System requirements

All code was tested in **Matlab R2014a**, but should be compatible with earlier versions of the software (at least back to 2011-2012). `mcadre` requires the following tools to be installed and *added to the Matlab path*:  

+ COBRA Toolbox: can be cloned from this GitHub [repo](https://github.com/opencobra/cobratoolbox)  
+ fastFVA<sup>(2)</sup>: available for download [here](https://notendur.hi.is/ithiele/software/fastfva.html)
+ fastcc<sup>(2)</sup>: available for download as part of the FASTCORE package [here](http://wwwen.uni.lu/recherche/fstc/life_sciences_research_unit/research_areas/systems_biology/software)

Assuming the above tools are installed and work correctly, the `mcadre` code should work the same, regardless of operating system (tested on Mac OS X 10.7.5).

<span style="font-size:9pt"><sup>2</sup>The fastFVA and fastcc methods perform the same function, so only one is required; the `fastFVA` function uses the free `glpk` solver, which comes with the COBRA Toolbox, while `fastcc` uses the commercial `cplex` solver (available with an academic license).</span>


## Input data

While the global model and confidence scores are included with the package, additional expression-based evidence inputs must be provided by the user. However, example inputs can be found in the `data/testInputs.mat` file.

### Generic/global reconstruction

This package comes with two ready-to-use global reconstruction inputs: `humanModel`, which represents *Recon 1*<sup>(3)</sup>, and `mouseModel`, which represents an updated and modified version of *iMM1415*<sup>(4)</sup>. Each model is saved in a `.mat` file, and can be loaded from the `data/` directory. These `.mat` files each also include a variable called `confidenceScores` representing literature/experimental-based confidence assigned to reactions in the generic model by the original authors.

<span style="font-size:9pt"><sup>3</sup>[Duarte et al. Global reconstruction of the human metabolic network based on genomic and bibliomic data. *Proc Natl Acad Sci USA.* 2007, **104(6)**:1777-82](http://www.ncbi.nlm.nih.gov/pubmed/17267599)  
<sup>4</sup>[Sigurdsson et al. A detailed genome-wide reconstruction of mouse metabolism based on human Recon 1. *BMC Syst Biol.* 2010, **4**:140](http://www.ncbi.nlm.nih.gov/pubmed/20959003)</span>

### Tissue-specific expression evidence

In addition to the global reconstruction input, `mcadre` requires two vectors containing tissue-specific evidence for individual gene expression, derived from one or more omics data sets:

+ `G`: a cellstr-formatted column vector, where each entry ***g*** is the human- or mouse-specific **Entrez ID** for gene ***g***  
+ `U`: a numeric column vector, where each entry represents the calculated "ubiquity score" for the corresponding gene ***g***

Ubiquity scores are calculated from binarized expression data, where ***U(g)*** is equal to the sum of tissue samples in which gene ***g*** is expressed, divided by the total number of samples. For example, if after binarization (**1** = expressed, **0** = not expressed), a gene is determined to be expressed in 4 of 5 tissue samples, the value of ***U*** for this gene would be **0.80**. 

**NOTE: **Given the variety of possible procedures for collecting and collating expression data, converting from probe (or other) ID to Entrez ID, and binarizing expression values, these steps are not included as part of the **mCADRE** package.

### High-confidence tissue-specific genes

An additional, optional input that may be included is `C_H_genes`. If included, this should be a cellstr column vector containing Entrez IDs for any genes with particularly strong evidence for activity in the selected tissue (e.g., detected expression of the corresponding protein product). Within the `rank_reactions` module of `mcadre`, the reactions corresponding to these genes will automatically elevated to the "core" set to avoid removal.

## Options

