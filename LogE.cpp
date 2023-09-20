#include <bitset>
#include <cmath>       /* tgamma */
#include <map>
#include <vector>

using namespace std;

#include "data.h"


/******************************************************************************/
/**************************   MODEL COMPLEXITY   ******************************/
/******************************************************************************/
double GeomComplexity_SubCM(unsigned int m);

/******************************************************************************/
/**************** Log-Evidence (LogE) of a sub-complete model  ****************/
/******************************************************************************/
// Compute the log-evidence of a sub-complete model based on m basis elements
// ! Kset must have been previously reduced to these m basis elements !
// This function is mainly used for call by `LogE_PartMCM`,
// but can also be used to compute the log-likelihood of a complete model
//
double LogE_SubC_forMCM(map<uint32_t, unsigned int > Kset, uint32_t m, unsigned int N)
{
  double LogE = 0;

  map<uint32_t, unsigned int >::iterator it;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;

  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) {cout << "problem Ks = 0 for mu_m = " << (it->first) << endl; }
    LogE += lgamma(Ks + 0.5);
  }  
  if (Ncontrol != N) { cout << "Error Likelihood function: Ncontrol != N" << endl;  }

  //LogE +=  ((1UL << m) - Kset.size()) * lgamma(0.5); // for all the states that are not observed
  //return LogE - GeomComplexity_SubCM(m) - lgamma( (double)( N + (1UL << (m-1)) ) );
  return LogE + lgamma((double)( 1UL << (m-1) )) - (Kset.size()/2.) * log(M_PI) - lgamma( (double)( N + (1UL << (m-1)) ) ); 
}

/******************************************************************************/
/*********  Log-Evidence (LogE) of a sub-complete part of a MCM   *************/
/******************************************************************************/
// Compute the log-evidence of the sub-complete part (of an MCM) defined by Ai.
// This function could be also used directly by the user
// to compute the log-evidence of a sub-complete model

double LogE_SubCM_Vect(vector<pair<uint32_t, unsigned int>> Kset_Vect, uint32_t Ai, unsigned int N)
{
  map<uint32_t, unsigned int > Kset_new;

  uint32_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset

//Build Kset_new:
  for (auto const& it : Kset_Vect)
  {
    s = ((it).first) & Ai;          // troncated state: take only the bits of s (=it.first) indicated by Ai
    Kset_new[s] += ((it).second);   // # of times s appears in the data set
  }

  return LogE_SubC_forMCM(Kset_new, bitset<n>(Ai).count(), N);
}

/******************************************************************************/
/****************************   LogE of a MCM   *******************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
bool check_partition(map<uint32_t, uint32_t> Partition);

double LogE_MCM_Vect(vector<pair<uint32_t, unsigned int>> Kset_Vect, map<uint32_t, uint32_t> Partition, unsigned int N)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogE = 0; 
    unsigned int rank = 0;
    map<uint32_t, uint32_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogE += LogE_SubCM_Vect(Kset_Vect, (*Part).second, N);
      rank += bitset<n>((*Part).second).count();
    }  
    return LogE - ((double) (N * (n-rank))) * log(2.);
  //}
}




