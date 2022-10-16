// To compile: g++ -std=c++11 -O3 main.cpp Data_Manipulation.cpp LogE.cpp LogL.cpp Complexity.cpp Best_MCM.cpp Basis_Choice.cpp MCM_info.cpp
// To run: time ./a.out
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <map>
#include <cmath>       /* tgamma */

using namespace std;

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"
#include "library.h"

/******************************************************************************/
/*******************************   main function   ****************************/
/******************************************************************************/
int main()
{  
  cout << "--->> Create OUTPUT Folder: (if needed) ";
  system( ("mkdir -p " + OUTPUT_directory).c_str() );
  cout << endl;

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "***********************************  Read the data:  **************************************";
  cout << endl << "*******************************************************************************************" << endl;
  unsigned int N=0; // will contain the number of datapoints in the dataset
  map<uint32_t, unsigned int> Nset = read_datafile(&N);


  cout << endl << "*******************************************************************************************";  
  cout << endl << "********************************** !! IMPORTANT !! ****************************************";
  cout << endl << "*******************************************************************************************";  
  cout << endl << "******************************  CHOICE OF THE BASIS:  *************************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << endl << "Choice of the basis for building the Minimally Complex Model (MCM):" << endl;

// *** Basis elements are written using the integer representation of the operator
// *** For instance, a basis element on the last two spin variable would be written: 
// ***      -->  Op = s1 s2           Spin operator
// ***      -->  Op = 000000011       Binary representation
// ***      -->  Op = 3               Integer representation   ( 000000011 = 3 )

  // *** One can simply use the original basis of the data:   // This is the most natural choice a priori
//   list<uint32_t> Basis_li = Original_Basis();

  // *** The basis can be specified by hand here:
  uint32_t Basis_Choice[] =  {3, 5, 9, 48, 65, 129, 272, 81, 1};    // Ex. This is the best basis for the "Shapes" dataset

  unsigned int m = sizeof(Basis_Choice) / sizeof(uint32_t);
  list<uint32_t> Basis_li;  Basis_li.assign (Basis_Choice, Basis_Choice + m); 

  // *** The basis can also be read from a file:
//   list<uint32_t> Basis_li = Read_BasisOp_IntegerRepresentation();
//   list<uint32_t> Basis_li = Read_BasisOp_BinaryRepresentation();

  // *** Print info about the Basis:
  PrintTerm_Basis(Basis_li);

  cout << "Number of spin variables, n=" << n << endl;
  cout << "Number of basis elements, m=" << Basis_li.size() << endl;
  if (Basis_li.size() > n) { cout << " -->  Error: the number 'm' of basis elements is larger than the size 'n' of the system." << endl;  }
  else { 
    cout << " -->  m <= n :  Everything seems fine." << endl;
    cout << "Make sure that the set of basis elements provided are orthogonal to each other." << endl;
  }

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "******************************  CHANGE THE BASIS OF THE DATA  *****************************"; 
  cout << endl << "**************************************  Build Kset:  **************************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << endl << "INFORMATION:\n\t To look for the best MCM on a chosen basis, the data needs first to be re-written in that basis.";
  cout << endl << "\t The following function re-write the histogram of the data in the basis specified above.";
  cout << endl << "\t If the specified basis is the original basis, then this transformation is not needed." << endl;

  cout << endl << "/!\\ IMPORTANT: If m<n:";
  cout << endl << "\t If the size 'm' of the basis is strictly smaller than the number 'n' of variables, ";
  cout << endl << "\t then the data will be truncated to the 'm' first basis elements." << endl;

  cout << endl << "Transform the data in the specified basis:" << endl;  
  map<uint32_t, unsigned int> Kset = build_Kset(Nset, Basis_li, false);

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "********************************  All Independent Models:  ********************************";
  cout << endl << "*******************************************************************************************" << endl << endl;

  cout << "Independent models in the new basis:" << endl;
  
  PrintInfo_All_Indep_Models(Kset, N);

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "**************************  All Successive Sub-Complete Models:  **************************";
  cout << endl << "*******************************************************************************************" << endl << endl;

  cout << "Sub-Complete models in the new basis:" << endl;

  PrintInfo_All_SubComplete_Models(Kset, N);


  cout << endl << "*******************************************************************************************"; 
  cout << endl << "*********************************  Define your own MCM:  **********************************";
  cout << endl << "*******************************************************************************************" << endl << endl;

  // *** The MCM can be specified by hand here:
  uint32_t MCM_Choice[] =  {384, 64, 32, 16, 8, 4, 2, 1};

  unsigned int k = sizeof(MCM_Choice) / sizeof(uint32_t);  // Number of parts
  map<uint32_t, uint32_t> MCM_Partition0 = Create_MCM(MCM_Choice, k);

  // *** The MCM can also be read from a file:
//  map<uint32_t, uint32_t> MCM_Partition0 = Read_MCMParts_BinaryRepresentation("./INPUT/Shapes_n9_MCM_Binary.dat");

  if(check_partition(MCM_Partition0))
  {
    PrintTerminal_MCM_Info(Kset, N, MCM_Partition0);
  }
  else { cout << "The set of 'parts' provided does not form a partition of the basis elements." << endl;  }


  cout << endl << "*******************************************************************************************";  
  cout << endl << "**************************************  EXAMPLES:  ****************************************";
  cout << endl << "*****************  The three following functions search for the BEST MCM ******************";
  cout << endl << "*************************  based on the BASIS specified above *****************************";
  cout << endl << "*******************************************************************************************" << endl; 

  cout << endl << "/!\\ IMPORTANT: prior to using the three following functions:"; 
  cout << endl << "\tThe choice of the basis must have been provided in the variable \'Basis_li\',";
  cout << endl << "\tand the data (initially stored in \'Nset\') must have been re-written in that basis using the function \'build_Kset()\'.";
  cout << endl << "\tThe variable \'Kset\' then contains the transformed data." << endl;

  cout << endl << "/!\\ SPECIAL CASE:  ANALYSIS IN THE ORIGINAL BASIS:"; 
  cout << endl << "\tThe data may also be analyzed untransformed, i.e. in its original basis.";
  cout << endl << "\tIn this case:\n\t\t-- there is no need to specify \'Basis_li\', nor to transformed the data with \'build_Kset()\',";  
  cout << endl << "\t\t-- the three following functions can be used directly with the variable \'Nset\' instead of \'Kset\'." << endl;

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "*******************************  Find the Best MCM:  **************************************";
  cout << endl << "***********************************  VERSION 1  *******************************************"; 
  cout << endl << "*******************************************************************************************";  
  cout << endl << "**********************  Compare all MCMs of a given rank 'r' ******************************";
  cout << endl << "*********************  based on the 'r' first basis Operators:  ***************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << endl << "/!\\ INFORMATION:"; 
  cout << endl << "\tThe following function searches for the best MCM among all the MCMs based on the 'r' first basis operators";
  cout << endl << "\ti.e., among only among MCMs of rank exactly equal to 'r'" << endl;

  cout << endl << "/!\\ IMPORTANT: Condition on the value of 'r':  r <= m <= n ";
  cout << endl << "\t'r' must be smaller or equal to the number 'm' of basis element provided, 'm=Basis_li.size()',";
  cout << endl << "\twhich must be smaller or equal to the number 'n' of spin variables." << endl << endl;

  int r1 = 9;
  double LogE_BestMCM1 = 0;

  if (r1 <= Basis_li.size())
  {
    map<uint32_t, uint32_t> MCM_Partition1 = MCM_GivenRank_r(Kset, N, &LogE_BestMCM1, r1, false);
    //cout << "\t Best LogE = " << LogE_BestMCM1 << endl;
    PrintTerminal_MCM_Info(Kset, N, MCM_Partition1);
  }
  else { cout << "The condition on the value of 'r' is not respected" << endl;  }

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "*******************************  Find the Best MCM:  **************************************";
  cout << endl << "***********************************  VERSION 2  *******************************************"; 
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "****************  Compare all MCMs of rank 'k', with 1 <= k <=r' **************************";
  cout << endl << "******************  based on the 'k' first basis Operators:  ******************************";
  cout << endl << "*******************************************************************************************" << endl;

  cout << endl << "/!\\ INFORMATION:"; 
  cout << endl << "\tThe following function searches among all the MCMs based on the 'k' FIRST basis operators provided, for all k in [1; r]" << endl;

  cout << endl << "/!\\ IMPORTANT: Condition on the value of 'r':  r <= m <= n ";
  cout << endl << "\t'r' must be smaller or equal to the number 'm' of basis element provided, 'm=Basis_li.size()',";
  cout << endl << "\twhich must be smaller or equal to the number 'n' of spin variables." << endl;

  cout << endl << "\tPlease check the function declaration for the default arguments." << endl << endl;

/******************************************************************************/
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
/******************************************************************************/

  int r2 = 9;
  double LogE_BestMCM2 = 0;

  if (r2 <= Basis_li.size())
  {
    map<uint32_t, uint32_t> MCM_Partition2 = MCM_AllRank_SmallerThan_r_Ordered(Kset, N, &LogE_BestMCM2, r2, false);
    //cout << "\t Best LogE = " << LogE_BestMCM2 << endl;
    PrintTerminal_MCM_Info(Kset, N, MCM_Partition2);
  }
  else { cout << "The condition on the value of 'r' is not respected" << endl;  }

  cout << endl << "*******************************************************************************************"; 
  cout << endl << "*******************************  Find the Best MCM:  **************************************";
  cout << endl << "***********************************  VERSION 3  *******************************************";   
  cout << endl << "*******************************************************************************************"; 
  cout << endl << "****************  Compare all MCMs of rank 'k', with 1 <= k <=r' **************************";
  cout << endl << "***********  based on any 'k' subset of the basis Operators provided  *********************";
  cout << endl << "*******************************************************************************************" << endl;
  
  cout << endl << "/!\\ INFORMATION:"; 
  cout << endl << "The following function searches among all MCMs based on ANY SUBSET of 'k' operators of the basis provided, for all k=1 to r" << endl;

  cout << endl << "/!\\ IMPORTANT: Conditions on the value of 'r':  r <= m <= n ";
  cout << endl << "\t'r' must be smaller or equal to the number 'm' of basis element provided, 'm=Basis_li.size()',";
  cout << endl << "\twhich must be smaller or equal to the number 'n' of spin variables." << endl;

  cout << endl << "\tPlease check the function declaration for the default arguments." << endl << endl;

/******************************************************************************/
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
/******************************************************************************/

  int r3 = 9;
  double LogE_BestMCM3 = 0;

  if (r3 <= Basis_li.size())
  {
    map<uint32_t, uint32_t> MCM_Partition3 = MCM_AllRank_SmallerThan_r_nonOrdered(Kset, N, &LogE_BestMCM3, r3, false);
    //cout << "\t Best LogE = " << LogE_BestMCM3 << endl;
    PrintTerminal_MCM_Info(Kset, N, MCM_Partition3);
  }
  else { cout << "The condition on the value of 'r' is not respected" << endl;  }


  return 0;
}


