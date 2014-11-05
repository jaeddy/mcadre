# mCADRE

This package contains my implementation of the **Matlab** code for the method described in Wang et al., 2012<sup>(1)</sup>. The method is designed to take a generic global metabolic reconstruction (e.g., *Human Recon 1*) plus a collection of tissue-specific gene expression evidence and produce a draft context-specific model for the selected tissue or cell type.

<sup>1) [Wang et al. Reconstruction of genome-scale metabolic models for 126 human tissues using mCADRE. *BMC Syst Biol.* 2012, **6**:153](http://www.ncbi.nlm.nih.gov/pubmed/23234303)</sup>


## Getting started

All code was tested in **Matlab R2014a**, but should be compatible with earlier versions of the software (at least back to 2011-2012). `mcadre` requires the following tools to be installed and *added to the Matlab path*:  

+ COBRA Toolbox: can be cloned from this GitHub [repo](https://github.com/opencobra/cobratoolbox)  
+ fastFVA<sup>(2)</sup>: available for download [here](https://notendur.hi.is/ithiele/software/fastfva.html)
+ fastcc<sup>(2)</sup>: available for download as part of the FASTCORE package [here](http://wwwen.uni.lu/recherche/fstc/life_sciences_research_unit/research_areas/systems_biology/software)

To run the **mCADRE** method, simply clone this repo and add the directory to your Matlab path. The script `run_mcadre` provides an example of how to call the main function `mcadre`. The command to run `mcadre` is structured as follows:

```
[PM, GM, C, NC, Z, model_C, pruneTime, cRes] = …
    mcadre(model, G, U, confidenceScores, salvageCheck, C_H_genes, method);
```

<sup>2) The fastFVA and fastcc methods perform the same function, so only one is required; the `fastFVA` function uses the free `glpk` solver, which comes with the COBRA Toolbox, while `fastcc` uses the commercial `cplex` solver (available with an academic license).</sup>


## Documentation

For details on running **mCADRE** and on various input data and options, please refer to the `Manual.md` file under `docs`. Further information about the algorithm and published results can be found in the original paper referenced above.


## Contact

Any errors or bugs with the code in this repository should be reported as issues on GitHub. Additional questions about using or further developing the method can be directed to the [Nathan Price lab](mailto:nprice@systemsbiology.org) at the Institute for Systems Biology. 

**NOTE:** no future releases or new features are planned for this package, so any users looking to modify the method are encouraged to fork this repo and do so on their own. 