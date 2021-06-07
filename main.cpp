//g++ -std=c++11 -O3 main.cpp Operations_OnData.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp
//time ./a.out
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <bitset>
#include <map>
#include <cmath>       /* tgamma */

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/******************     READ and TRANSFORM DATA    ****************************/
/******************************************************************************/

/*** READ DATA and STORE data in Nset:    *************************************/
/******************************************************************************/
map<uint32_t, unsigned int> read_datafile(unsigned int *N); // filename to specify in data.h

/*** READ BASIS from a FILE:    ***********************************************/
/******************************************************************************/
list<uint32_t> Read_BasisOp_BinaryRepresentation();   // filename to specify in data.h
list<uint32_t> Read_BasisOp_IntegerRepresentation(); 

/*** Original Basis:    ***********************************************/
/******************************************************************************/
list<uint32_t> Original_Basis();   // return the original basis, i.e., {s1, s2, ..., sn}

/*** Print Basis Info in the Terminal:    *************************************/
/******************************************************************************/
void PrintTerm_Basis(list<uint32_t> Basis_li);

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
/***************** Log-LIKELIHOOD (LogL), Log-EVIDENCE (LogE) *****************/
/***************************  and COMPLEXITY   ********************************/
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
/******************************   Find Best MCM   *****************************/
/******************************************************************************/
// *** Compute all partitions of a set using Algorithm H:

// *** Version 1: Compare all the MCM of rank r, 
// ***            based on the r first elements of the basis used to build Kset:
map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int r, unsigned int N, double *LogE_best);

// *** Version 2: Compare all the MCM 
// ***            based on the r first elements of the basis used to build Kset
// ***            for all r=1 to basis.size()  
map<uint32_t, uint32_t> MCM_allRank(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best);

// *** Print information about the MCM specified in `MCM_Partition`:
void PrintTerminal_MCM_Info(map<uint32_t, unsigned int > Kset, unsigned int N, map<uint32_t, uint32_t> MCM_Partition);

/******************************************************************************/
/*******************************   main function   ****************************/
/******************************************************************************/
int main()
{  
  cout << "--->> Create OUTPUT Folder: (if needed) ";
  system( ("mkdir -p " + OUTPUT_directory).c_str() );
  cout << endl;

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "******************************  Choice of the basis:  *************************************";
  cout << endl << "*******************************************************************************************" << endl;
// *** Choice of the basis for building the Minimally Complex Model (MCM):

// *** Basis elements are written using the integer representation of the operator
// *** For instance, a basis element on the last two spin variable would be written: 
// ***      -->  Op = s1 s2           Spin operator
// ***      -->  Op = 000000011       Binary representation
// ***      -->  Op = 3               Integer representation   ( 000000011 = 3 )

  // *** The basis can be specified by hand here:
  uint32_t Basis_Choice[] = {36, 10, 3, 272, 260, 320, 130, 65, 4};    // Ex. This is the best basis for the USSC dataset

  unsigned int m = sizeof(Basis_Choice) / sizeof(uint32_t);
  list<uint32_t> Basis_li;  Basis_li.assign (Basis_Choice, Basis_Choice + m); 

  // *** Basis can also be read from a file:
  // list<uint32_t> Basis_li = Read_BasisOp_IntegerRepresentation();

  // *** Or one can use the original basis:
  // list<uint32_t> Basis_li = Original_Basis();

  // *** Print info about the Basis:
  PrintTerm_Basis(Basis_li);

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "***********************************  Read the data:  **************************************";
  cout << endl << "*******************************************************************************************" << endl;
  unsigned int N=0; // will contain the number of datapoints in the dataset
  map<uint32_t, unsigned int> Nset = read_datafile(&N);

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "*********************************  Change the data basis   ********************************"; 
  cout << endl << "**************************************  Build Kset:  **************************************";
  cout << endl << "*******************************************************************************************" << endl;
  // *** Transform the data in the specified in Basis_SCModel[];
  map<uint32_t, unsigned int> Kset = build_Kset(Nset, Basis_li, false);

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "************************************  All Indep Models:  **********************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << "Independent models in the new basis:" << endl;
  
  map<uint32_t, uint32_t> Partition_Indep;  uint32_t Op = 1;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_Indep[i] = Op;
    cout << "Added Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_Indep, N) << " \t LogL = " << LogL_MCM(Kset, Partition_Indep, N) << endl;    
    Op = Op << 1;
  }

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "**************************  All Successive Sub-Complete Models:  **************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << "Sub-Complete models in the new basis:" << endl;

  map<uint32_t, uint32_t> Partition_SC;  Op = 1;
  Partition_SC[0] = 0;
  for (uint32_t i = 0 ; i<n; i++)
  {
    Partition_SC[0] += Op;
    cout << "Added Op = " << Op << " \t LogE = " << LogE_MCM(Kset, Partition_SC, N) << " \t LogL = " << LogL_MCM(Kset, Partition_SC, N) << endl;    
    Op = Op << 1;
  }

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "********************************  Compare all MCModels:  **********************************";
  cout << endl << "*******************************************************************************************" << endl;

// *** Search among all MCM based on the r first basis operators (models of rank exactly equal to `r`):
  int r=9;
  double LogE_BestMCM = 0;
  map<uint32_t, uint32_t> MCM_Partition = MCM_GivenRank_r(Kset, r, N, &LogE_BestMCM);
//cout << "\t Best LogE = " << LogE_BestMCM << endl;

// *** Search among all MCM based on the r first basis operators,
// ***  where   1 <= r <= m, where m is the size of Basis_li:
  //map<uint32_t, uint32_t> MCM_Partition = MCM_allRank(Kset, N, &LogE_BestMCM);

// *** 
  PrintTerminal_MCM_Info(Kset, N, MCM_Partition);

  return 0;
}


