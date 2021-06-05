# Exhaustive search for the best Minimally Complex Spin Models
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



