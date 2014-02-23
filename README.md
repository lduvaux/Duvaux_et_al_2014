AphidCNV
========

Copy number variation and speciation in Pea Aphid. A collaboration with [Ludovic Duvaux](l.duvaux@sheffield.ac.uk).


Description
------

This is the source code we have used for our paper (TODO: LINK).


Usage
------

You need `R(>=3.0)`.
And the following packages:

```R
install.packages("ape");
install.packages("epicalc");
install.packages("foreign");
install.packages("ggplot2");
install.packages("lme4");
install.packages("MASS");
install.packages("methods");
install.packages("microbenchmark");
install.packages("multicore");
install.packages("MuMIn");
install.packages("nnet");
install.packages("optimalCaptureSegmentation");
install.packages("parallel");
install.packages("randomForest");
```

You also need the unofficial package `optimalCaptureSegmentation`
you can download it at http://bioinformatics.nki.nl/ocs/downloads/optimalCaptureSegmentation_0.9-4.tar.gz
And to install it :

```R
install.packages("/PATH/TO/THE/PKG/optimalCaptureSegmentation_0.9-4.tar.gz", repos =NULL, type="source")
```


Code organisation
-------

The pipeline is entirely implemented into the `src` directory.
The succession of processing steps is orchestrated by a master `src/Makefile`.
Each processing steps has its own folder and local `Makefile`.
The master `Makefile` essentially calls `sub-Makefiles` in the good order.
After if finishes, each processing steps generate an transitory `*.Rdata` file that can be loaded by another step.
