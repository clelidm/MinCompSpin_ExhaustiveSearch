# Exhaustive search for the best Minimally Complex Spin Models

This repository contains a code developed for the paper on *Statistical Inference of Minimally Complex Models* available in [arXiv:2008.00520]. The code performs an exhaustive search for the best Minimally Complex Spin Model (MC-Spin Model) on a given basis. 

The code go through all possible MC-spin models of a given rank `r`, where an MC-spin model is defined by a given partition of the `r` basis operators provided. The comparison between models is based on their evidence (posterior probability that the model produces the data, integrated over the parameter values — see paper). The selected model is the one with the largest evidence.

One big advantage of the family of models (the MC-spin models) is that the computation of the evidence doesn’t require fitting the parameters of the model, which allows to significantly accelerate the comparison between models. The limiting factor of this code is the exhaustive search for the best model in a model space that is exponentially increasing with `r`.

For a given number of `r` of basis elements, the number of possible partitions the `r` basis operators is equal to the Bell number of `r`, which grows exponentially with `r`. Besides, the running time of the program also grows linearly with the number of different states observed in the system (so approximatively linearly with the number of datapoints in the dataset). For these reasons, this code is good for use on small systems, typically with r <~15 variables. Note that for r~15, the code may take several days to perform the exhaustive search.

Contains a code that performs an exhaustive search for the best MC-Spin Model on a given basis

This repository contains codes developped for the paper on *Statistical Inference of Minimally Complex Models* described in [arXiv:2008.00520].

## Requirements

The code uses the C++11 version of C++.

## Usage

### Set global variables in the file `data.h`:

Before compiling specify the following global variables:
 - `const unsigned int n`, with the number of variables of the dataset; This number can also be the number of basis operators in the provided basis; it must not be smaller than the number of basis operators.
 - `const string OUTPUT_directory`, with the name of the output directory; All the generated files will be placed in this folder.
 - `const string datafilename`, with the location and name of the binary datafile.
 - (Optional) `const string basis_IntegerRepresentation_filename`, with the location and name of the input file containing the basis element written in the integer representation.
 - (Optional) `const string basis_BinaryRepresentation_filename`,  with the location and name of the input file containing the basis element written in the binary representation.

### Specify the spin basis:

The basis can be specified by hand directly at the beginning of the `int main()` function in `uint32_t Basis_Choice[]`, or by using an input file.

#### Structure of the basis:

#### Reading the basis from an input file:

### Examples:

All useful functions that can be called from `int main()` are declared at the beginning of the `main.cpp` file and described in the section "Available functions" below. For hands-on and simple tests of the program, check the examples in the function `int main()` of the `main.cpp` file.

### Compile:

`g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp`

## Input and Output files:

## Available functions:



