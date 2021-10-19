# Exhaustive search
# for the best Minimally Complex spin Models (MCM)

This program allows to **uncover AND model community structures in binary data**, while **taking into account possible high order patterns** of data in the detection of the communities (i.e., possible high-order correlations between the variables). The idea of the algorithm is based on performing statistical inference with a family of spin models (maximum entropy models for binary data) that have minimal information theoretic complexity. These models are called Minimally Complex Models (MCM). The selected model can be used as a generative model for data. Details can be found in Ref. [1].

----

## Description

This repository contains a code developed for the paper Ref. [1] on *Statistical Inference of Minimally Complex Models* available in [arXiv:2008.00520](https://arxiv.org/abs/2008.00520). The code performs an exhaustive search for the best Minimally Complex spin Model (MCM) in a chosen basis for the data. The basis is chosen by the user.

To do so, the code go through all possible MCMs of a given rank `r`, where an MCM is defined by a partition of the `r` basis operators provided (see Ref. [1]). The comparison between models is based on their evidence (posterior probability that the model generates the data, integrated over its parameter values). The selected model is the one with the largest evidence.

One big advantage of this family of models (the MCMs) is that the computation of the evidence doesn’t require fitting the parameters of the model, which allows to significantly accelerate the comparison between models. The limiting factor of this code is the exhaustive search for the best model, as the space of models is exponentially increasing with `r`.

Provided `r` of basis elements, the number of possible partitions of this set of `r` elements is equal to the Bell number of `r`, which grows roughly exponentially with `r`. Besides, the running time of the program also grows linearly with the number of different states observed in the dataset (so approximatively linearly with the number of datapoints in the dataset). For these reasons, this code is good for use on small systems, typically with `r <~15` variables. Note that for `r~15`, the code will take several days to perform the exhaustive search.

To efficiently generate all possible set partitions of these `r` operators, we used the algorithm described in Ref. [2] and in Ref. [3] (Algorithm E) that generates the partitions in Gray-code order. Thus, the code goes through all possible partitions by changing only single bit for each new partition.

For larger systems (`r>~15`), please check the greedy algorithm available [here](https://github.com/clelidm/MinCompSpin_Greedy). 

[1]  C. de Mulatier, P. P. Mazza, M. Marsili, *Statistical Inference of Minimally Complex Models*, [arXiv:2008.00520](https://arxiv.org/abs/2008.00520)

[2]  Ehrlich, Gideon. *Loopless algorithms for generating permutations, combinations, and other combinatorial configurations.* Journal of the ACM (JACM) 20.3 (1973): 500-513.

[3]  D.E. Knuth, *The Art of Computer Programming*, Volume 4, Combinatorial Algorithms: Part 1.693 (Addison-Wesley Professional), (2011).

## Requirements

The code uses the C++11 version of C++.

**To compile:** `g++ -std=c++11 -O3 main.cpp Data_Manipulation.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp MCM_info.cpp`

**To execute:** `./a.out`

## Examples

All the functions useful for the user are declared in the file `library.h` and can be called from the function `int main()` in the file `main.cpp`. They are described in the sections below.

For hands-on and simple tests of the program, check the examples in the function `int main()` of the `main.cpp` file. In the `INPUT` folder, we provided the binary dataset `Dataset_Shapes_n9_N1e5.dat`, which is the dataset used in the section "Boolean description of a binary dataset" of Ref. [1].

Note that, for illustration purposes, the example in `int main()` run three different versions of the search algorithm (see section "Find the Best MCM"). To run with more variables, please choose only one of these versions (we advise you to use the function `MCM_GivenRank_r()`, see `Version 1` in the examples). 

## License

This code is an open source project under the GNU GPLv3.

----

# Usage

## Input and Output files

**Input files:** Input files can be placed in the `INPUT` folder; you must provide the following input files:
 - a binary datafile. The name of the datafile must be specified in `data.h` in the variable `datafilename`. Datapoints must be written as binary strings of 0’s and 1’s encoded on at least `n` bits (with no spaces between the bits), where `n` is the number of spin variables specified in `data.h`. The file must contain one datapoint per line — see example file `Dataset_Shapes_n9_N1e5.dat` in the `INPUT` folder.
 - (Optional) a file containing the basis elements (see section "Reading the basis from an input file” below).

**Output files:**
All the output files will be stored in the output folder, whose name is specified in `data.h` in the variable `OUTPUT_directory`.

## Set the global variables, in the file `data.h`

Before compiling specify the following global variables:
 - `const unsigned int`**`n`**, with the number of variables of the dataset. This number must be smaller or equal to the number of columns in the input dataset. If it is smaller, the program will only read the `n` first columns of the dataset (from the left). This number must be larger or equal to the number `m` of basis elements provided in the `main()` function.
 - `const string`**`datafilename`**, with the location and name of the input binary datafile.
 - (Optional) `const string`**`basis_IntegerRepresentation_filename`**, with the location and name of the input file containing the basis element written in the integer representation (see section "Reading the basis from an input file” below).
 - (Optional) `const string`**`basis_BinaryRepresentation_filename`**,  with the location and name of the input file containing the basis element written in the binary representation (see section "Reading the basis from an input file” below).
 - `const string`**`OUTPUT_directory`**, with the name of the output directory. All the generated files will be placed in this folder.

## Specify the spin basis

The following functions are defined in `Basis_Choice.cpp`.

The elements of the basis on which to build the Minimally Complex Model (MCM) have to be specified by the user before running the program.

In the code, a basis is stored as a list of 32-bit integers `list<uint32_t>`**`Basis`**, where each integer defines a spin operator (see explanation below, in the section `Structure of the basis`). There are three ways to specify the basis:

 - The basis can be written “by hand” directly at the beginning of the `int main()` function, in the variable `uint32_t Basis_Choice[]` (see section `Defining the basis manually` below).
 
 - The basis can be provided in an input file (see section `Reading the basis from an input file` below).
 
 - The basis can simply be the original basis in which the data is already written. If you don’t know which basis to use, you can run the minimally complex model algorithm on the "original basis" of the data, which is the basis in which the dataset is written. This can be done by using the function `list<uint32_t>`**`Original_Basis`**`()` to define the basis.

In general, we advise you to use the basis in which the dataset is the closest to be generated from an independent model (see discussion in Ref. [1]). Finding this basis can be done using the greedy search algorithm available separately [here](https://github.com/clelidm/BinaryData_HighOrderInteractions_GreedyAlgo).

### Structure of the basis:

**Basis:** The basis elements are spin operators that are all independent from each others (see Ref. [1]). You can use the function (*to come*) to check if the elements you have specified in `list<uint32_t> Basis` actually form a basis, i.e. if the set is only composed of independent operators.

The number of elements in the basis can be at most equal to the number `n` of variables in the system. Note that if you provide more than `n` elements in your basis, then the elements in the set you provided are not independent, and, at most, there is `n` of them that are independent. 

The rank `r` of the basis you provide can be smaller than `n`. In this case, the code will automatically truncate the data to reduce it to the sub-space defined by the `r` first basis elements specified.

**Spin operators:** Each element of the basis must be a spin operator, where a spin operator is defined as the product of a subset of spin variables. For instance, `Op = s1 * s2` is a spin operator (graphically, it is associated to a pairwise interactions between `s1` and `s2`); `Op = s1*s2*s3` is also a spin operator (associated to a three-body interaction between `s1`, `s2` and `s3`).

In the code, spin operators are encoded by a binary integer on `n` bits, where `n` is the number of spin variables in the system (which you must define in the file `data.h`). The binary representation of a given operator has a bit `1` for each spin included in the operator, and `0` for all the other spins. Importantly, spin variables are numbered from the right to the left. 
For instance, take a system of `n=9` spins and the operator `Op = s1 s2`, this operator would be represented in the code by the binary number `Op = 000000011`. Finally, to simplify the definition of a spin operator in the code, you can directly use the integer corresponding to this binary number. For instance, to defined the operator `Op = s1 s2`, you can use the binary representation `Op = 000000011` or the integer representation `Op = 3`.
>      Example: the three representations of a spin operator: 
>      -->  Op = s1 s2           Spin operator
>      -->  Op = 000000011       Binary representation
>      -->  Op = 3               Integer representation   ( 000000011 = 3 )

The convention taken for writing the spin operators is to label the spin variables from the right to the left in the binary representation. Note that this doesn't change the ordering of the spin variables from their order in the dataset, i.e., the first variable on the left in the binary representation of the operators is the same as the first variable on the left in the dataset provided as input file.

However, an important point, to be able to interpret the results of the program, is that we adopted the same convention for the new basis: the first operator provided in the basis will correspond to the variable `sigma1`, which will be the variable the most on the right in the transformed dataset (see the example in the next section).

### Defining the basis manually (see example at the beginning of the `int main()` function):

The basis can be specified by hand directly at the beginning of the `int main()` function in `uint32_t Basis_Choice[]`. In this case, we advise you to use the integer representation of the basis operators. For example, `uint32_t Basis_Choice[] = {36, 10, 3, 272, 260, 320, 130, 65, 4}` defines a basis with `9` independent spin operators. Below we give the different representations for these operators: first, the integer representation in the original basis; second, the corresponding binary representation; third, the corresponding spin operator; and finally representations of these operators in the new basis:
>      36     000100100     s3 s6     -->>  new basis :     000000001     1       sigma1
>      10     000001010     s2 s4     -->>  new basis :     000000010     2       sigma2
>      3      000000011     s1 s2     -->>  new basis :     000000100     4       sigma3
>      272    100010000     s5 s9     -->>  new basis :     000001000     8       sigma4
>      260    100000100     s3 s9     -->>  new basis :     000010000     16      sigma5
>      320    101000000     s7 s9     -->>  new basis :     000100000     32      sigma6
>      130    010000010     s2 s8     -->>  new basis :     001000000     64      sigma7
>      65     001000001     s1 s7     -->>  new basis :     010000000     128     sigma8
>      4      000000100     s3        -->>  new basis :     100000000     256     sigma9

### Reading the basis from an input file:

The following functions allow you to define the basis from an input file:

 - `list<uint32_t>`**`Read_BasisOp_BinaryRepresentation`**`(string Basis_binary_filename = basis_BinaryRepresentation_filename)`, for operators written with the binary representation (see the example file `Dataset_Shapes_n9_Basis_Binary.dat` in the `INPUT` folder). The location of the file must be given as an argument of the function. By default the function will read the file specified in `data.h` in the variable `basis_BinaryRepresentation_filename`.

 - `list<uint32_t>`**`Read_BasisOp_IntegerRepresentation`**`(string Basis_integer_filename = basis_IntegerRepresentation_filename)`, for operators written with the integer representation (see example file `Dataset_Shapes_n9_Basis_Integer.dat` in the `INPUT` folder). The location of the file must be given as an argument of the function. By default the function will read the file specified in `data.h` in the variable `basis_IntegerRepresentation_filename`.

For both functions, operators must be written in the file in one single column. The operator at the top of the column will correspond to the variable `sigma1` in the new basis (i.e. to the bit the most to the right), the second to `sigma2`, etc.

### Printing the basis in the terminal:
To print information about a basis in the terminal, use the function `void`**`PrintTerm_Basis`**`(list<uint32_t> Basis_li)`.

## Read and Transform the Input Data:

The following functions are defined in `Data_Manipulation.cpp`.

### Read the input dataset

The function `map<uint32_t, unsigned int>`**`read_datafile`**`(unsigned int *N, string filename = datafilename)` reads the dataset with location and name provided in argument (`string filename`). By default the function will read the file specified in the variable `const string datafilename` in `data.h`. The dataset is then stored in the structure `map<uint32_t, unsigned int>`**`Nset`** that maps each observed state to the number of times they occur in the dataset. Note that each state of the system is encoded as an `n`-bit integer.

### Re-write the dataset in the new basis

The function `map<uint32_t, unsigned int>`**`build_Kset`**`(map<uint32_t, unsigned int> Nset, list<uint32_t> Basis, bool print_bool=false)` changes the basis of the dataset from its original basis (or the one in which `Nset`, provided as an argument, is written) to the basis provided as an argument in `Basis`. It is possible to print this new map (i.e., the frequency of occurrence of each state in the new basis) by changing the default value of `print_bool` to `true`.

## Find the Best MCM:

The following functions are defined in `Best_MCM.cpp`.

### Encoding MCMs:

The code only compares MCMs that correspond to partitions of the basis operators provided.
Once the dataset is converted in the new basis, these MCMs are simply encoded on `r` digits, where `r` is the number of basis operators on which one want to perform the search for the best MCM (see next section).

In the code, partitions of the `r` digits are encoded in two different ways:

 - **Version a:**  a vector of `r` integers, which is used for enumerating partitions in Gray-code order. In this representation, spin variables with the same digit belong to the same part. This representation is also used for printing MCMs in files. See example below.

 - **Version b:**  a list of 32-bits integers, which is used for computing the log-evidence and the log-likelihood of MCMs (see section below), and which is the representation in which the best MCM is printed to the terminal. In this representation, each integer encodes a part of the partition: in the binary representation of each integer, the variables with a bit equal to `1` belong to the part. See example below.

>      Example, for a system with 9 binary variables:
>      
>      Version a: The partition 001000112 has 3 parts:
>                            - one with spins s4, s5, s6, s8, s9, corresponding to the digit `0`;
>                            - one with spins s2, s3, s7, corresponding to the digit `1`;
>                            - one with the spin s1 alone, corresponding to the digit `2`.
>                            
>      Version b: The same partition can be encoded as the list of three integers, 
>                  440, 70, 1, with binary representations:
>                            - 440 = 110111000, this part contains the spins s4, s5, s6, s8, s9;
>                            -  70 = 001000110, this part contains the spins s2, s3, s7;
>                            -   1 = 000000001, this part contains the spin s1 alone.

### Exhaustive search for the best MCM:

Three functions are available to perform an exhaustive search for the best MCM (i.e., the MCM with the largest log-evidence `logE`):

 - **Function 1:** The function **`MCM_GivenRank_r`** compares all the MCMs of rank `r`, based on the `r` first elements of the new basis (i.e., the basis used to build Kset). The total number of these models is given by the Bell number of `r`, denoted `Bell(r)`. See declaration:
 >        map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false)

 - **Function 2:** The function **`MCM_AllRank_SmallerThan_r_Ordered`** compares all the MCMs based on the `k` first elements of the new basis for all `k=1 to r`. The total number of these model is given by the sum for `k=1` to `r` of `Bell(k)`. See declaration:
>        map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_Ordered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false)

 - **Function 3:** The function **`MCM_AllRank_SmallerThan_r_nonOrdered`** compares all the MCMs based on **any** `k` elements of the new basis for all `k=1 to r`. The total number of these model is given by the sum for `k=1` to `r` of `[n choose k] x Bell(k)`. See declaration:
>        map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_nonOrdered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false)

These three functions enumerate all possible partitions of a set using variantes of the algorithm described in Ref. [2] and [3]. The algorithm efficiently generates all set partitions in Gray-code order.

For all three functions: 
 - the default value of `r` is the number `n` of spin variables;
 - it is possible to print in a file the values of the log-Evidence of **all the tested models**. To do so, change the value of the input variable `print_bool` to true (the defaut value is `false`).

**Recommendations:**  These three functions were created for you to test what happens if you compare the models for the three different cases. However, for general use, we recommend that, once you have defined which are the `r` operators on which you want to search for the best MCM, you directly run the search among MCMs exactly based on these `r` operators (i.e., using the function 1). We suggest to reduce beforehand the selection to the `r` most relevant basis operators (similarly to a dimensionality reduction step).

### Print information about your model

The function `void`**`PrintTerminal_MCM_Info`**`(map<uint32_t, unsigned int >Kset, unsigned int N, map<uint32_t, uint32_t> MCM_Partition)` prints in the terminal information about the MCM given as an argument in `MCM_Partition`:

 - **For the whole MCM,** the function prints the number of Independent Components of the MCM (i.e., the number of partitions -- or communities -- of the MCM, see ref. [1]), as well as the log-likelihood (`LogL`), the parameter complexity (`C_param`), the geometric complexity (`C_geom`), the total complexity (`C_tot`), the Minimum Description Length (`MDL`) and the log-evidence (`LogE`). 
See Ref. [Entropy 2018, 20(10), 739](https://www.mdpi.com/1099-4300/20/10/739) for the definition of the complexity of spin models (in connection with the Minimum Description Length Principle).

 - **For each independent part in the MCM,** the function prints the integer and binary representation of the part (see section "Encoding MCMs" above), as well as  the following information for the MCM with a single community based on that part: the log-likelihood (`LogL`), the parametric complexity (`C_param`), the geometric complexity (`C_geom`), the total complexity (`C_tot`), and the log-evidence (`LogE`).

### Define your own MCM

Any MCM can also be specified by the user:

  - **directly in the program**, with the function `map<uint32_t, uint32_t>`**`Create_MCM`**`(uint32_t MCM_table[], int k)` (see the example in the `int main()` function). In this case, we advise you to use the integer representation of the parts. For example, `uint32_t MCM_Choice[] = {384, 64, 32, 16, 8, 4, 2, 1}` defines an MCM with `8` independent parts, based on `n=9` spins. This is the best MCM obtained in the example of section "Boolean description of a binary dataset" in Ref. [1]. 

  - **using an input file**, with the function `map<uint32_t, uint32_t>`**`Read_MCMParts_BinaryRepresentation`**`(string MCM_binary_filename)`. Parts must be written with the binary representation, one part per line (see the example file `Dataset_Shapes_n9_MCM_Binary.dat` in the INPUT folder). The location of the file must be given as an argument of the function (see the example in `int main()`).

Below we give the different representations for the parts defined by `uint32_t MCM_Choice[] = {384, 64, 32, 16, 8, 4, 2, 1}`: in the first column, the integer representation in the original basis; in the second column, the corresponding binary representation: 
>      384    110000000
>      64     001000000
>      32     000100000
>      16     000010000 
>      8      000001000  
>      4      000000100  
>      2      000000010  
>      1      000000001 

You can check that the model (i.e., list of parts) that you have provided properly defines an MCM by calling the function `bool`**`check_partition`**`(map<uint32_t, uint32_t> Partition)`. This function will return `false` if there is an overlap between the parts.

## Likelihood, Complexity and Evidence:

The following functions are defined in `LogL.cpp`, `Complexity.cpp` and `LogE.cpp`.

Users can also get **specific information about an MCM** with the following functions:
- `double`**`LogL_MCM`**`(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false)` returns the log-likelihood of the MCM defined by `Partition`;
- `double`**`LogE_MCM`**`(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false)` returns the log-evidence of the MCM defined by `Partition`;
- `double`**`Complexity_MCM`**`(map<uint32_t, uint32_t> Partition, unsigned int N, double *C_param, double *C_geom)` place the parameter complexity and the geometric complexity of the MCM model defined in `Partition` respectively at the addresses `*C_param` and `*C_geom`. Finally, the function returns the total complexity of the model.

Users can also get **specific information about any subcomplete part (SC-part)** of an MCM with the functions:
 - `double`**`LogL_SubCM`**`(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, bool print_bool = false)` returns the log-likelihood of the SC-part, where `Kset` is the dataset written in the new basis, and where `Ai` is the binary representation of the SC-part (see section "Encoding MCMs").
 - `double`**`LogE_SubCM`**`(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, bool print_bool = false)` returns the log-evidence of the SC-part, where `Kset` is the dataset written in the new basis, and where `Ai` is the binary representation of the SC-part (see section "Encoding MCMs").
 - `double`**`ParamComplexity_SubCM`**`(unsigned int m, unsigned int N)` returns the model complexity of the SC-part due to the number of parameters in the part ("parameter complexity"); this is the first order complexity term in the Minimum Description Length principle (this term is of the order of `O(log N)` where `N` is the number of datapoints -- see Ref. [1]).
 - `double`**`GeomComplexity_SubCM`**`(unsigned int m)` returns the geometric complexity of the SC-part; it is the second order complexity term in the Minimum Description Length principle (which is of the order of `O(1)` -- see Ref. [1]).

