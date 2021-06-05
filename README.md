# Exhaustive search for the best Minimally Complex Spin Models

This repository contains a code developed for the paper on *Statistical Inference of Minimally Complex Models* available in [arXiv:2008.00520](https://arxiv.org/abs/2008.00520). The code performs an exhaustive search for the best Minimally Complex Spin Model (MC-Spin Model) on a given basis. 

The code go through all possible MC-spin models of a given rank `r`, where an MC-spin model is defined by a given partition of the `r` basis operators provided (see paper). The comparison between models is based on their evidence (posterior probability that the model produces the data, integrated over the parameter values). The selected model is the one with the largest evidence.

One big advantage of this family of models (the MC-spin models) is that the computation of the evidence doesn’t require fitting the parameters of the model, which allows to significantly accelerate the comparison between models. The limiting factor of this code is the exhaustive search for the best model as the space of models is exponentially increasing with `r`.

For a given number of `r` of basis elements, the number of possible partitions the `r` basis operators is equal to the Bell number of `r`, which grows exponentially with `r`. The running time of the program also grows linearly with the number of different states observed in the system (so approximatively linearly with the number of datapoints in the dataset). For these reasons, this code is good for use on small systems, typically with `r <~15` variables. Note that for `r~15`, the code will take several days to perform the exhaustive search.

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

### Specify the spin basis (see functions in `Basis_Choice.cpp`):

The element of the basis for building the Minimally Complex Model (MCM) has to be specified by the user before compiling the code.

The basis can be written “by hand” directly at the beginning of the `int main()` function in `uint32_t Basis_Choice[]`, or by using an input file (see subsections below).

If you don’t know which basis to use, you can run the minimally complex model algorithm on the original basis in which the dataset is written. This can be done by using the function `list<uint32_t> Original_Basis()` to define the basis.

In general, we advise you to use the basis in which the dataset is the closest to be generated from an independent model (see discussion in the paper). Finding this basis can be done using the heuristic/exhaustive search algorithm available separately *here*.

In the code, a basis is stored as a list of integers `list<uint32_t> Basis`, where each integer defines a basis operator (see explanation below).

#### Structure of the basis:

Basis elements are spin operators that are all independent from each others (see paper). You can use the function (to come) to check if the elements you have specified in `list<uint32_t> Basis` actually form a basis, i.e. if the set is only composed of independent operators.

Each element of the basis must be a spin operator. A spin operator is the product of a subset of spin variables (see paper). For instance, Op = s1 * s2 is a spin operator (physically, it is associated to a pairwise interactions between s1 and s2); Op = s1*s2*s3 is also a spin operator (this time associated to a three-body interaction between s1, s2 and s3).

In the code, spin operators are encoded on a binary number of `n` bits, where `n` is the number of spin variables in the system (which you must define in the file `data.h`). The binary representation of a given operator has a bit `1` for each spin included in the operator, and `0` for all the other spins. Importantly, spin variables are numbered from the right to the left. 
For instance, take the operator Op = s1 s2, this operator would be represented in the code by the binary number `Op = 000000011` (for `n=9`).

>      -->  Op = s1 s2           Spin operator
>      -->  Op = 000000011       Binary representation
>      -->  Op = 3               Integer representation   ( 000000011 = 3 )

Finally, to simplify the definition of a spin operator, you can directly use the integer corresponding to the binary number. For instance, to defined the operator Op = s1 s2, you can use the binary representation Op = 000000011 or the integer representation Op = 3.

The number of elements in the basis must be at most equal to the number `n` of variables in the system. Note that the rank `r` of the basis can also be smaller than `n`. In this case, the code will automatically truncate the data to reduce it to the sub-space defined by the `r` specified basis elements.

#### Reading the basis from an input file:

### Examples:

All useful functions that can be called from `int main()` are declared at the beginning of the `main.cpp` file and described in the section "Available functions" below. For hands-on and simple tests of the program, check the examples in the function `int main()` of the `main.cpp` file.

### Compile:

`g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp`

## Input and Output files:

## Available functions:



