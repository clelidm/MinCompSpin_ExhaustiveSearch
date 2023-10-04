#include <bitset>
#include <cmath>       /* tgamma */
#include <map>
#include <vector>

using namespace std;

#include "data.h"

/******************************************************************************/
/************************ Build Kset for a single ICC  ************************/
/******************************************************************************/
map<uint32_t, unsigned int > build_Kset_ICC(vector<pair<uint32_t, unsigned int>> Kset, uint32_t Ai)
{
    map<uint32_t, unsigned int > Kset_ICC;

    uint32_t s;        // state
    unsigned int ks=0; // number of time state s appears in the dataset

  //Build Kset_ICC:
    for (auto const& it : Kset)
    {
      s = ((it).first) & Ai;          // troncated state: take only the bits of s (=it.first) indicated by Ai
      Kset_ICC[s] += ((it).second);   // # of times s appears in the data set
    }

    return Kset_ICC;
}

/***********************************************************************************************************************/
/***********************************************************************************************************************/
/**************************************************   LOG-E   **********************************************************/
/***********************************************************************************************************************/
/***********************************************************************************************************************/
/*
double LogE_CM(vector<pair<uint32_t, unsigned int>> Kset, uint32_t m, unsigned int N)
{
  double LogE = 0;

  for (auto& it : Kset)
  {
    LogE += lgamma((it.second) + 0.5);
  }  

  return LogE + lgamma((double)( 1UL << (m-1) )) - (Kset.size()/2.) * log(M_PI) - lgamma( (double)( N + (1UL << (m-1)) ) ); 
}*/

/******************************************************************************/
/*********  Log-Evidence (LogE) of an ICC part of a MCM   *********************/
/******************************************************************************/
// Compute the LogE of the ICC defined by Ai;
//    i.e. of the sub-part of an MCM identififed by Ai;
// this function doesn't account of the contribution to LogE due to the non-modeled spins (i.e. N*log(2) per spin)

double LogE_ICC(vector<pair<uint32_t, unsigned int>> Kset, uint32_t Ai, unsigned int N)
{
  map<uint32_t, unsigned int > Kset_ICC = build_Kset_ICC(Kset, Ai);  /// Question: make an exception for the CM? i.e. if Ai = 1111...11 ??

  uint32_t m = bitset<n>(Ai).count();

  double LogE = 0;

  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;

  map<uint32_t, unsigned int >::iterator it;

  for (it = Kset_ICC.begin(); it!=Kset_ICC.end(); ++it)
  {
    Ks = (it->second);  
    Ncontrol += Ks;
    if (Ks == 0) {cout << "problem Ks = 0 for mu_m = " << (it->first) << endl; }
    LogE += lgamma(Ks + 0.5);
  }  
  if (Ncontrol != N) { cout << "Error Likelihood function: Ncontrol != N" << endl;  }

  //LogE +=  ((1UL << m) - Kset.size()) * lgamma(0.5); // for all the states that are not observed
  //return LogE - GeomComplexity_ICC(m) - lgamma( (double)( N + (1UL << (m-1)) ) );
  return LogE + lgamma((double)( 1UL << (m-1) )) - (Kset_ICC.size()/2.) * log(M_PI) - lgamma( (double)( N + (1UL << (m-1)) ) ); 
}

/******************************************************************************/
/****************************   LogE of a MCM   *******************************/
/******************************************************************************/
// check if *Partition* is an actual partition of the basis elements, 
//   i.e., that no basis element appears in more than 1 part of the partition.
//   i.e., that each basis element only appears in a single part of the partition.
//bool check_partition(map<uint32_t, uint32_t> Partition);

double LogE_MCM(vector<pair<uint32_t, unsigned int>> Kset, map<uint32_t, uint32_t> Partition, unsigned int N)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogE = 0; 
    unsigned int rank = 0;
    map<uint32_t, uint32_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogE += LogE_ICC(Kset, (*Part).second, N);
      rank += bitset<n>((*Part).second).count();
    }  
    return LogE - ((double) (N * (n-rank))) * log(2.);
  //}
}


/***********************************************************************************************************************/
/***********************************************************************************************************************/
/**************************************************   LOG-L   **********************************************************/
/***********************************************************************************************************************/
/***********************************************************************************************************************/

/******************************************************************************/
/**************** Log-likelihood (LogL) of a Complete Model  ******************/
/******************************************************************************/
// Compute the log-likelihood of a Complete Model on Kset:

double LogL_CM(vector<pair<uint32_t, unsigned int>> Kset, unsigned int N)
{
  double LogL = 0;
  double Nd = N;

  for (auto& it : Kset)
  {
    LogL += ((it.second) * log((double) (it.second) / Nd) ); 
  } 

  return LogL;
}

/******************************************************************************/
/*********  Log-Likelihood (LogL) of an ICC part of a MCM   *******************/
/******************************************************************************/
// Compute the LogL of the ICC defined by Ai;
//    i.e. of the sub-part of an MCM identififed by Ai;
// this function doesn't account of the contribution to LogL due to the non-modeled spins (i.e. N*log(2) per spin)

double LogL_ICC(vector<pair<uint32_t, unsigned int>> Kset, uint32_t Ai, unsigned int N)
{
  map<uint32_t, unsigned int > Kset_ICC = build_Kset_ICC(Kset, Ai);

  double LogL = 0;

  map<uint32_t, unsigned int >::iterator it;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;
  double Nd = N;

  for (it = Kset_ICC.begin(); it!=Kset_ICC.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) {cout << "problem Ks = 0 for mu_m = " << (it->first) << endl; }
    LogL += (Ks * log((double) Ks / Nd) );
  }
  if (Ncontrol != N) { cout << "Error in function 'LogLikelihood_SCforMCM': Ncontrol != N" << endl;  }

  return LogL;
}

/******************************************************************************/
/******************** Log-likelihood (LogL) of a MCM  *************************/
/******************************************************************************/

double LogL_MCM(vector<pair<uint32_t, unsigned int>> Kset, map<uint32_t, uint32_t> Partition, unsigned int N)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogL = 0; 
    unsigned int rank = 0;
    map<uint32_t, uint32_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogL += LogL_ICC(Kset, (*Part).second, N);
      rank += bitset<n>((*Part).second).count();
    }  
    return LogL - ((double) (N * (n-rank))) * log(2.);
  //}
}

