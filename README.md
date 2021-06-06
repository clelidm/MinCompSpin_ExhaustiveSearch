# Exhaustive search for the best Minimally Complex Spin Models

This repository contains a code developed for the paper on *Statistical Inference of Minimally Complex Models* available in [arXiv:2008.00520](https://arxiv.org/abs/2008.00520). The code performs an exhaustive search for the best Minimally Complex Spin Model (MC-Spin Model) on a given basis. 

The code go through all possible MC-spin models of a given rank `r`, where an MC-spin model is defined by a given partition of the `r` basis operators provided (see paper). The comparison between models is based on their evidence (posterior probability that the model produces the data, integrated over the parameter values). The selected model is the one with the largest evidence.

One big advantage of this family of models (the MC-spin models) is that the computation of the evidence doesn’t require fitting the parameters of the model, which allows to significantly accelerate the comparison between models. The limiting factor of this code is the exhaustive search for the best model as the space of models is exponentially increasing with `r`.

For a given number of `r` of basis elements, the number of possible partitions the `r` basis operators is equal to the Bell number of `r`, which grows exponentially with `r`. The running time of the program also grows linearly with the number of different states observed in the system (so approximatively linearly with the number of datapoints in the dataset). For these reasons, this code is good for use on small systems, typically with `r <~15` variables. Note that for `r~15`, the code will take several days to perform the exhaustive search.

## Requirements

The code uses the C++11 version of C++.

### Compile:

`g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp`

## Examples:

All useful functions that can be called from `int main()` are declared at the beginning of the `main.cpp` file and described in the sections below. For hands-on and simple tests of the program, check the examples in the function `int main()` of the `main.cpp` file. (Add details on the dataset taken as an example).

----

# Usage

## Set global variables in the file `data.h`:

Before compiling specify the following global variables:
 - `const unsigned int n`, with the number of variables of the dataset; This number can also be the number of basis operators in the provided basis; it must not be smaller than the number of basis operators.
 - `const string OUTPUT_directory`, with the name of the output directory; All the generated files will be placed in this folder.
 - `const string datafilename`, with the location and name of the binary datafile.
 - (Optional) `const string basis_IntegerRepresentation_filename`, with the location and name of the input file containing the basis element written in the integer representation.
 - (Optional) `const string basis_BinaryRepresentation_filename`,  with the location and name of the input file containing the basis element written in the binary representation.


## Specify the spin basis (see functions in `Basis_Choice.cpp`):

The element of the basis for building the Minimally Complex Model (MCM) has to be specified by the user before compiling the code.

The basis can be written “by hand” directly at the beginning of the `int main()` function in `uint32_t Basis_Choice[]`, or by using an input file (see subsections below).

If you don’t know which basis to use, you can run the minimally complex model algorithm on the original basis in which the dataset is written. This can be done by using the function `list<uint32_t> Original_Basis()` to define the basis.

In general, we advise you to use the basis in which the dataset is the closest to be generated from an independent model (see discussion in the paper). Finding this basis can be done using the heuristic/exhaustive search algorithm available separately *here*.

In the code, a basis is stored as a list of integers `list<uint32_t> Basis`, where each integer defines a basis operator (see explanation below).

### Structure of the basis:

Basis elements are spin operators that are all independent from each others (see paper). You can use the function (*to come*) to check if the elements you have specified in `list<uint32_t> Basis` actually form a basis, i.e. if the set is only composed of independent operators.

Each element of the basis must be a spin operator. A spin operator is the product of a subset of spin variables (see paper). For instance, `Op = s1 * s2` is a spin operator (physically, it is associated to a pairwise interactions between `s1` and `s2`); `Op = s1*s2*s3` is also a spin operator (this time associated to a three-body interaction between `s1`, `s2` and `s3`).

In the code, spin operators are encoded on a binary number of `n` bits, where `n` is the number of spin variables in the system (which you must define in the file `data.h`). The binary representation of a given operator has a bit `1` for each spin included in the operator, and `0` for all the other spins. Importantly, spin variables are numbered from the right to the left. 
For instance, take the operator `Op = s1 s2`, this operator would be represented in the code by the binary number `Op = 000000011` (for `n=9`). Finally, to simplify the definition of a spin operator in the code, you can directly use the integer corresponding to this binary number. For instance, to defined the operator `Op = s1 s2`, you can use the binary representation `Op = 000000011` or the integer representation `Op = 3`.

> Example: the three representations of a spin operator: 
>      -->  Op = s1 s2           Spin operator
>      -->  Op = 000000011       Binary representation
>      -->  Op = 3               Integer representation   ( 000000011 = 3 )

Finally, the number of elements in the basis must be at most equal to the number `n` of variables in the system. Note that the rank `r` of the basis can also be smaller than `n`. In this case, the code will automatically truncate the data to reduce it to the sub-space defined by the `r` specified basis elements.

Note, in the example above, that the convention taken for writing the spin operators is to label the spin variables from the right to the left in the binary representation. This is just a convention and doesn't change the ordering of the spin variables from their order in the dataset, i.e., the first variable on the left in the binary representation of the operators is the same as the first variable on the left in the dataset provided as input file.

However, an important point to be able to interpret the results of the program, is that we adopted the same convention for the new basis: the first operator provided in the basis will correspond to the variable `sigma1`, which will be the most on the right in the transformed dataset.

### Defining the basis manually (at beginning of the `int main()` function):
The basis can be specified by hand directly at the beginning of the `int main()` function in `uint32_t Basis_Choice[]`. In this case, you we advise you to use the integer representation of the basis operators. In the example provided in the `int main()` function: `uint32_t Basis_Choice[] = {36, 10, 3, 272, 260, 320, 130, 65, 4}` defines a basis with `9` independent operators. Here are the different representations for these spin operators (first, the integer representation; second, the binary representation; finally the corresponding spin operators)
>   36     000100100     s3 s6
>   10 	   000001010     s2 s4
>   3 	    000000011     s1 s2
>   272 	  100010000     s5 s9
>   260 	  100000100     s3 s9
>   320 	  101000000     s7 s9
>   130 	  010000010     s2 s8
>   65 	   001000001     s1 s7
>   4 	    000000100     s3

### Reading the basis from an input file (see `Basis_Choice.cpp`):

The following functions allow you to define a basis from an input file.
 - `list<uint32_t> Read_BasisOp_BinaryRepresentation()`, if operators are written using the binary representation (see example file in the `INPUT` folder); the location of the file must be specified in `data.h` in the variable `basis_BinaryRepresentation_filename`.
 - `list<uint32_t> Read_BasisOp_IntegerRepresentation()`, if operators are written using the integer representation (see example file in the `INPUT` folder); the location of the file must be specified in `data.h` in the variable `basis_IntegerRepresentation_filename`.

For any of these two functions, operators should be written in one single column in the file. 

### Print the basis in the terminal (see `Basis_Choice.cpp`):
To check the information about a basis, you can print it in the terminal using the function `void PrintTerm_Basis(list<uint32_t> Basis_li)`.

## Read the input dataset:

The function `map<uint32_t, unsigned int> read_datafile(unsigned int *N)` reads the dataset available at the location specified in the variable `const string datafilename` in `data.h`. The dataset is then stored in the a structure `map<uint32_t, unsigned int> Nset` that map each observed states to the number of times they occur in the dataset.

## Transform the input dataset in the new basis:

The function `map<uint32_t, unsigned int> build_Kset(map<uint32_t, unsigned int> Nset, list<uint32_t> Basis, bool print_bool=false)` changes the basis of the dataset from its original basis (or the one in which `Nset` provided as an argument is written) to the one provided as the argument `Basis`.

## Find the Best MC-Spin Model:

Two functions are available to perform a search among MC-Spin Models:
 - `map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int r, unsigned int N, double *LogE_best)` compares all the MCM of rank r, based on the `r` first elements of the new basis `Basis_li` (the one used to build Kset);
 - `map<uint32_t, uint32_t> MCM_allRank(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best)` compares all the MCM based on the `r` first elements of the new basis `Basis_li` for all `r=1` to the size of `Basis_li`;

## Print information about your model


## Input and Output files:

### Input files:
Input files must be stored in the `INPUT` folder, you must provide the following input files:
— a binary datafile. The name of the datafile must be specified in `data.h` in the variable `datafilename`. Datapoints must be written as binary strings of 0’s and 1’s encoded on at least `n` bits (with no spaces between the bits), where `n` is the number of spin variables specified in `data.h`. The file must contain one datapoint per line — see example in the INPUT folder.
— (optional) a basis input file (see section "Reading the basis from an input file” above).

### Output files:
All the output files will be stored in the output folder whose name is specified in `data.h`.

## Available functions:



