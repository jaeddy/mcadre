# mCADRE

This package contains my implementation of the **Matlab** code for the method described in Wang et al., 2012<sup>(1)</sup>. The method is designed to take a generic global metabolic reconstruction (e.g., *Human Recon 1*) plus a collection of tissue-specific gene expression evidence and produce a draft context-specific model for the selected tissue or cell type.

<span style=font-size:10pt><sup>1</sup>[Wang et al. Reconstruction of genome-scale metabolic models for 126 human tissues using mCADRE. *BMC Syst Biol.* 2012, **6**:153](http://www.ncbi.nlm.nih.gov/pubmed/23234303)</span>

## Getting started

All code was tested in **Matlab R2014a**, but should be compatible with earlier versions of the software (at least back to 2011-2012). `mcadre` requires the following tools to be installed and added to the Matlab path:  

+ COBRA Toolbox: can be cloned from this GitHub [repo](https://github.com/opencobra/cobratoolbox)  
+ fastFVA<sup>(2)</sup>: available for download [here](https://notendur.hi.is/ithiele/software/fastfva.html)
+ fastcc<sup>(2)</sup>: available for download as part of the FASTCORE package [here](http://wwwen.uni.lu/recherche/fstc/life_sciences_research_unit/research_areas/systems_biology/software)

To run the **mCADRE** method, simply clone this repo and add the directory to your Matlab path. The script `run_mcadre` provides an example of how to call the main function `mcadre`.

<span style=font-size:10pt><sup>2</sup>The fastFVA and fastcc methods perform the same function, so only one is required; the `fastFVA` function uses the `glpk` solver, which comes with the COBRA Toolbox, while `fastcc` uses the `cplex` solver.</span>

## Inputs

### Generic/global reconstruction

This package comes with two ready-to-use global reconstruction inputs: `humanModel`, which represents *Recon 1*<sup>(3)</sup>, and `mouseModel`, which represents an updated and modified version of *iMM1415*<sup>(4)</sup>. Each model is saved in a `.mat` file, and can be loaded from the `data/` directory. These `.mat` files each also include a variable called `confidenceScores` representing literature/experimental-based confidence assigned to reactions in the generic model by the original authors.

### Tissue-specific expression evidence