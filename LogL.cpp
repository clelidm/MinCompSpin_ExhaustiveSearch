#include <bitset>
#include <cmath>       /* tgamma */
#include <map>

using namespace std;

#include "data.h"

/******************************************************************************/
/**************** Log-likelihood (LogL) of a complete model  ******************/
/******************************************************************************/
// Compute the log-likelihood of a complete model on Kset:
// This function is mainly used for call by `LogL_SC_PartMCM`,
// but can also be used to compute the log-likelihood of a complete model
//
double LogL_CM(map<uint32_t, unsigned int > Kset, unsigned int N)
{
  double LogL = 0;

  map<uint32_t, unsigned int >::iterator it;
  unsigned int Ncontrol = 0; // for control
  unsigned int Ks = 0;
  double Nd = N;

  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    Ks = (it->second);  Ncontrol += Ks;
    if (Ks == 0) {cout << "problem Ks = 0 for mu_m = " << (it->first) << endl; }
    LogL += (Ks * log((double) Ks / Nd) );
  }
  if (Ncontrol != N) { cout << "Error in function 'LogLikelihood_SCforMCM': Ncontrol != N" << endl;  }

  return LogL;
}

/******************************************************************************/
/***************************** Log-likelihood (LogL) **************************/
/***********************   of a sub-complete part of a MCM   ******************/
/******************************************************************************/
// Compute the log-likelihood of the sub-complete part (of an MCM) defined by Ai.
// This function could be also used directly by the user
// to compute the log-likelihood of a sub-complete model

double LogL_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, bool print_bool = false)
{
  map<uint32_t, unsigned int>::iterator it;
  map<uint32_t, unsigned int > Kset_new;

  uint32_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset

  if (print_bool)  { 
  cout << endl << "--->> Build Kset for SC Model based on "  << Ai << " = " << bitset<n>(Ai) << " for MCM.." << endl;
  }
//Build Kset_new:
  for (it = Kset.begin(); it!=Kset.end(); ++it)
  {
    s = it->first;      // initial state s 
    ks = it->second;    // # of times s appears in the data set
    if (print_bool)  {  cout << s << ": \t" << bitset<n>(s) << " \t" ;  }

    s &= Ai;   // troncated state: take only the bits indicated by Ai
//    sig_m = bitset<m>(bitset<m>(mu).to_string()).to_ulong(); //bitset<m>(mu).to_ulong(); // mu|m
    if (print_bool)  {  cout << s << ": \t" << bitset<n>(s) << endl; }

    Kset_new[s] += ks;
    //Kset[mu_m].second.push_back(make_pair(mu, N_mu));
  }
  if (print_bool)  {  cout << endl;  }

  return LogL_CM(Kset_new, N);
}

/******************************************************************************/
/******************** Log-likelihood (LogL) of a MCM  *************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
bool check_partition(map<uint32_t, uint32_t> Partition);

double LogL_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false)
{
  //if (!check_partition(Partition)) {cout << "Error, the argument is not a partition." << endl; return 0;  }

  //else
  //{
    double LogL = 0; 
    unsigned int rank = 0;
    map<uint32_t, uint32_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      LogL += LogL_SubCM(Kset, (*Part).second, N);
      rank += bitset<n>((*Part).second).count();
    }  
    return LogL - ((double) (N * (n-rank))) * log(2.);
  //}
}


