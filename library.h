#include <sstream>
#include <list>
#include <map>


/******************************************************************************/
/******************************************************************************/
/*********************     SPECIFY the BASIS    *******************************/
/******************************************************************************/
/******************************************************************************/

/*** Original Basis:    ***********************************************/
/******************************************************************************/
list<uint32_t> Original_Basis();   // return the original basis, i.e., {s1, s2, ..., sn}

/*** READ BASIS from a FILE:    ***********************************************/
/******************************************************************************/
list<uint32_t> Read_BasisOp_BinaryRepresentation(string Basis_binary_filename = basis_BinaryRepresentation_filename);   // default filename to specify in data.h
list<uint32_t> Read_BasisOp_IntegerRepresentation(string Basis_integer_filename = basis_IntegerRepresentation_filename); 

/*** Print Basis Info in the Terminal:    *************************************/
/******************************************************************************/
void PrintTerm_Basis(list<uint32_t> Basis_li);




/******************************************************************************/
/******************************************************************************/
/******************     READ and TRANSFORM DATA    ****************************/
/******************************************************************************/
/******************************************************************************/

/*** READ DATA and STORE data in Nset:    *************************************/
/******************************************************************************/
map<uint32_t, unsigned int> read_datafile(unsigned int *N, string filename = datafilename); // filename to specify in data.h

/*** DATA CHANGE of BASIS:    *************************************************/
/******************************************************************************/
// *** Build Kset with the following definitions:
// *** mu_m = states of the systems written in the basis specified in `list<uint32_t> Basis`
// *** Kset[sig_m] = Number of times the state mu_m appears in the transformed dataset
//
// *** Rem: the new basis can have a lower dimension then the original dataset; 
// *** in which case the function will reduce the dataset to the subspace defined by the specified basis.
map<uint32_t, unsigned int> build_Kset(map<uint32_t, unsigned int> Nset, list<uint32_t> Basis, bool print_bool=false);




/******************************************************************************/
/******************************************************************************/
/***************** Log-LIKELIHOOD (LogL), Log-EVIDENCE (LogE) *****************/
/***************************  and COMPLEXITY   ********************************/
/******************************************************************************/
/******************************************************************************/

/****************   for a sub-Complete Model (SubCM)   ************************/
/******************************************************************************/
// *** the SubCM is the one specified in Ai;
// *** Ai must be an integer encoded on at least n bits, where each 1 indicates the basis elements included in the part:
// *** For ex. Ai = 01001 is encoded on n=5 basis elements, and element Op1 and Op4 belong to the part;
// *** Rem: Basis elements are ordered from the right to the left.

double LogL_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, bool print_bool = false);
double LogE_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, bool print_bool = false);

// *** Complexity of a SC model based on m basis Operators: m >= 1. Rem: C_geom(m=1) = log(pi):
double GeomComplexity_SubCM(unsigned int m);                  // Geometric complexity
double ParamComplexity_SubCM(unsigned int m, unsigned int N); // Complexity due to the number of parameters

/******************   for a Complete Model (CM)   *****************************/
/******************************************************************************/
double LogL_CM(map<uint32_t, unsigned int > Kset, unsigned int N);

/****************************    for a MCM     ********************************/
/******************************************************************************/
double LogL_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false);
double LogE_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false);
double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, double *C_param, double *C_geom);




/******************************************************************************/
/******************************************************************************/
/********************   DEFINE MCMs and PRINT INFO   **************************/
/******************************************************************************/
/******************************************************************************/

// *** Define an MCM by hand:
map<uint32_t, uint32_t> Create_MCM(uint32_t MCM_table[], int k);

// *** Define an MCM from a file; Each part must be encoded in a binary number over n spins:
map<uint32_t, uint32_t> Read_MCMParts_BinaryRepresentation(string MCM_binary_filename);

// *** Check that the provided model corresponds to a partition of the basis variables (i.e. properly defines an MCM):
bool check_partition(map<uint32_t, uint32_t> Partition);

// *** Print information about the MCM specified in `MCM_Partition`:
void PrintTerminal_MCM_Info(map<uint32_t, unsigned int > Kset, unsigned int N, map<uint32_t, uint32_t> MCM_Partition);

// *** Create successive independent models defined on the new basis, and print the corresponding information:
void PrintInfo_All_Indep_Models(map<uint32_t, unsigned int> Kset, unsigned int N);

// *** Create successive Sub-complete models defined on the new basis, and print the corresponding information:
void PrintInfo_All_SubComplete_Models(map<uint32_t, unsigned int> Kset, unsigned int N);




/******************************************************************************/
/******************************************************************************/
/***************************   Find Best MCM   ********************************/
/******************************************************************************/
/******************************************************************************/
// *** Compute all partitions of a set using Algorithm H:

/******************************************************************************/
// *** Version 1: Compare all the MCMs of rank r, 
// ***            based on the r first elements of the basis used to build Kset:
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false);

/******************************************************************************/
// *** Version 2:  
// ***            Compare all the MCMs 
// ***            based on the k first elements of the basis used to build Kset
// ***            for all k=1 to r, where r <= basis.size() 
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_Ordered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false);

/******************************************************************************/
// *** Version 3:  
// ***            Compare all the MCMs based on any subset of k elements 
// ***            of the of r first elements of the basis used to build Kset
// ***            for all k=1 to r, where r <= basis.size() 
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_nonOrdered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false);

