AphidCNV
========

Copy number variation and speciation in Pea Aphid. A collaboration with [Ludovic Duvaux](l.duvaux@sheffield.ac.uk).


Description
------

This is the source code we have used for our paper (TODO: LINK).


Usage
------

You need `R(>=3.0)`.
And the folowing packages:

* `TODO`
* `randomForest`
* `TODO`
* optimalCapt (TODO : clone!!)



Code organisation
-------

The pipeline is entirely implemented into the `src` directory.
The succession of processing steps is orchestrated by a master `src/Makefile`.
Each processing steps has its own folder and local `Makefile`.
The master `Makefile` essentially calls `sub-Makefiles` in the good order.
After if finishes, each processing steps generate an transitory `*.Rdata` file that can be loaded by another step.
