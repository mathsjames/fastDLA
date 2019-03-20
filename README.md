# fastDLA

This repository provides a program for simulating off-grid diffusion limited aggregation (DLA) due to Peter Meakin [TODO].

The code runs in O(n polylog(n)) time (within a polylog(n) factor of optimal) and is explained and analysed in a forthcoming paper. The header file alone can be included in a C++ project or the cpp file can be compiled to generate clusters and save them to file.

## Usage

Downloading the whole repository and compiling fastDLA.cpp will provide a binary file which can then be executed with the command

fastDLA n filename seed noiseReductionFactor

(Where "fastDLA" should be replaced with the name of the binary), generates a cluster of n particles and saves their x and y co-ordinates (in binary format) into the file "filename". The arguments seed and noiseReducitonFactor are optional. seed sets the seed of the random number generator allowing for reproducible output and is set using the time if no seed is provided. noiseReductionFactor can be ignored unless the user wishes to simulate the noise reduced diffusion limited aggregation of [TODO], in which case this sets the noise reduction factor to use.

Alternatively the header file fastDLA.hpp can be included in a C++ project, it provides a single class FastCluster whose constructor takes arguments n, seed and noiseReductionFactor (as before the latter two are optional) and whose variable points is an array of std::complex<double> which after construction contains the points of the cluster.

## Licence

TODO

## Future Alterations

This code may be updated over time to improve performance but we will endeavour to keep the usage the same (of course old versions will remain available in any case). Feature suggestions are welcome, though no promise is made to implement any requests if there is a research need and it is not too complicated to implement I would be happy to dive back into the code to stop someone else from having to do so.

## Performance

The most recent performance tests were done on a TODO and the runtimes were as follows