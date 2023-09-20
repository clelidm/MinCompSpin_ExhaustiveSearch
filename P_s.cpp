#include <bitset>
#include <cmath>       /* tgamma */
#include <map>
#include <list>
#include <fstream>
#include <vector>

using namespace std;

#include "data.h"


/********************************************************************/
/**************************    Structure    *************************/
/********************************************************************/
struct Proba {
    uint32_t s;         // state in the original basis
    uint32_t sig;       // state in the new basis

    double P_D_s = 1.;    // empirical probability of s  
    double P_D_sig = 1.;  // empirical probability of sig --> this should be the same then s if r=n, i.e. if the MCM models all the n spins
    double P_MCM = 1.;  // model probability of s 
};

/******************************************************************************/
/****************   Return Kset over a chosen sub-basis b_a    ****************/
/******************************************************************************/
map<uint32_t, unsigned int> Build_Kset_ba_Vect(vector<pair<uint32_t, unsigned int>> Kset_Vect, uint32_t Ai)
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

  return Kset_new;
}

/******************************************************************************/
/****************************      Check partition     ************************/
/******************************************************************************/
//check if *Partition* is an actual partition of the basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.
pair<bool, uint32_t> check_partition(map<uint32_t, uint32_t> Partition);

/******************************************************************************/
/*****************   Compute the contribution to P_MCM(s)   *******************/
/******************  due to the sub-CM defined by Kset_ba   *******************/
/******************************************************************************/

void update_proba_MCM(map<uint32_t, Proba> &all_P, map<uint32_t, unsigned int> Kset_ba, uint32_t Ai, unsigned int N)
{
  map<uint32_t, Proba>::iterator it_P;

  uint32_t s, sig;        // states
  unsigned int ks=0;      // number of time state s appear in the dataset
  double Nd = (double) N;

//Contribution of the basis "ba" of spins to the model probability P_MCM(s):
  for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)
  {
    s = it_P->first;      // initial state s 
    sig = s & Ai;         // troncated state: take only the bits indicated by Ai
    all_P[s].P_MCM *= Kset_ba[sig]/Nd;
  }
}

/******************************************************************************/
/*******************      Compute the model probability     *******************/
/*********************  for the MCM constructed from Kset  ********************/
/************************   with the given partition   ************************/
/******************************************************************************/
// This function can be used directly on the original basis, by replacing Kset by Nset:

map<uint32_t, Proba> P_sig_Vect(vector<pair<uint32_t, unsigned int>> Kset_Vect, map<uint32_t, uint32_t> Partition, unsigned int N) // Probabilities in the sigma basis
{
  // Fill in the data probability:
  map<uint32_t, Proba> all_P;
  double Nd = (double) N;

  uint32_t s;        // state
  unsigned int ks=0; // number of time state s appear in the dataset

  // Check partition:
  pair<bool, uint32_t> Is_partition = check_partition(Partition);
  uint32_t rank = Is_partition.second;

  if (!Is_partition.first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    double pre_factor = 1./((double) (1 << (n-rank))); 


//Build Kset_new:
    for (auto const& it : Kset_Vect)
    {
      s = (it).first;          // initial state s
      ks = (it).second;    // # of times s appears in the data set

      all_P[s].P_D_s = ((double) ks)/Nd;
      all_P[s].P_MCM = pre_factor;
    }

    // Compute the Kset over each part: Kset_ba:
    map<uint32_t, unsigned int> Kset_ba;
    map<uint32_t, uint32_t>::iterator Part;

    for (Part = Partition.begin(); Part != Partition.end(); Part++)
    {
      Kset_ba = Build_Kset_ba_Vect(Kset_Vect, (*Part).second);         // (*Part).second) = Ai = integer indicated the spin elements included in b_a
      update_proba_MCM(all_P, Kset_ba, (*Part).second, N);
      Kset_ba.clear();
    }  
  }

  return all_P;
}


void PrintFile_StateProbabilites_NewBasis_Vect(vector<pair<uint32_t, unsigned int>> Kset_Vect, map<uint32_t, uint32_t> MCM_Partition, unsigned int N, string filename = "Result")
{
  // Probabilities in the sigma basis:
  map<uint32_t, Proba> P_all = P_sig_Vect(Kset_Vect, MCM_Partition, N);
  map<uint32_t, Proba>::iterator it_P;

  string Psig_filename = filename + "_DataVSMCM_Psig.dat";

  fstream file_P_sig((OUTPUT_directory + Psig_filename), ios::out);
  file_P_sig << "## 1:sig \t 2:P_D(sig) \t 3:P_MCM(sig)" << endl;

  for (it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    file_P_sig << bitset<n>(it_P->first) << "\t " << (it_P->second).P_D_s << "\t " << (it_P->second).P_MCM << endl;
  }

  file_P_sig.close();
}

/******************************************************************************/
/*******************      Compute the model probability     *******************/
/**************************  for the MCM constructed   ************************/
/***************   in a given basis, with a given partition   *****************/
/******************************************************************************/
uint32_t transform_mu_basis(uint32_t mu, list<uint32_t> basis);

map<uint32_t, Proba> P_s_Vect(vector<pair<uint32_t, unsigned int>> Nset_Vect, list<uint32_t> Basis, map<uint32_t, uint32_t> Partition, unsigned int N) // Probabilities in the sigma basis
{
  double Nd = (double) N;

  // Fill in the data probability:
  map<uint32_t, Proba> all_P;

  if (!check_partition(Partition).first) {cout << "Error, the argument is not a partition: the function returned an empty map for P[s]." << endl; }
  else
  { 
    // Build Kset:
    map<uint32_t, unsigned int > Kset;
    uint32_t s;        // initial state
    uint32_t sig_m;    // transformed state and to the m first spins
    unsigned int ks=0; // number of time state s appear in the dataset

  //Build Kset and fill in P[s] from the data:
    cout << "--->> Build Kset and fill in P[s] from the data..." << endl;
    for (auto const& it : Nset_Vect)
    {
      s = (it).first;          // original state s
      sig_m = transform_mu_basis(s, Basis);   // new state

      // Fill in Kset:
      Kset[sig_m] += (it).second;    // += ks = number of time state s appear in the dataset

      // Fill in P[s] empirical and value of transformed state:
      all_P[s].P_D_s = ((double) ((it).second))/Nd; // = ks/Nd
      all_P[s].sig = sig_m;
    }

    // convert map to a vector
    vector<pair<uint32_t, unsigned int>> Kset_Vect(Kset.size());
    int i=0;
    for (auto& my_pair : Kset)
    {
      Kset_Vect[i]=my_pair;
      i++;
    }

  // Compute the model probability of the MCM based on Kset using "Partition":
    cout << "--->> Compute P[s] for the chosen MCM..." << endl << endl;
    map<uint32_t, Proba> all_P_sig = P_sig_Vect(Kset_Vect, Partition, N);   // Compute P[s] for the MCM in the new basis

  // Report the values of P_MCM in the original all_P:
    map<uint32_t, Proba>::iterator it_P;
    for (it_P = all_P.begin(); it_P!=all_P.end(); ++it_P)           // Report P[s] for the MCM in the original basis
    {
      sig_m = (it_P->second).sig;
      (it_P->second).P_MCM = all_P_sig[sig_m].P_MCM;
      (it_P->second).P_D_sig = all_P_sig[sig_m].P_D_s;  // not necessary a probability distribution anymore
    }
  
    all_P_sig.clear();
    Kset.clear();
  }

  return all_P;  
}

/******************************************************************************/
/*****************      PRINT FILE: INFO about an MCM     *********************/
/******************************************************************************/

void PrintFile_MCM_Info(list<uint32_t> Basis, map<uint32_t, uint32_t> MCM_Partition, string filename = "Result")
{
  //***** PRINT BASIS: 
  fstream file_MCM_info((OUTPUT_directory + filename + "_MCM_info.dat"), ios::out);

  file_MCM_info << "## sig_vec = states in the chosen new basis (ideally the best basis), defined by the basis operators:" << endl;
  int i = 1;
  for (list<uint32_t>::iterator it = Basis.begin(); it != Basis.end(); it++)
  {
    file_MCM_info << "##\t sig_" << i << " = " << bitset<n>(*it) << " = " << (*it) << endl; i++;
  } file_MCM_info << "##" << endl;

  // Print info about the model -- Print MCM:
  file_MCM_info << "## The MCM Partition is defined on the sig_vec basis by the following Parts:" << endl;

  //***** PRINT MCM: 
  i = 1;
  for (map<uint32_t, uint32_t>::iterator it = MCM_Partition.begin(); it != MCM_Partition.end(); it++)
  {    
    uint32_t Part = (*it).second;
    file_MCM_info << "##\t MCM_Part_" << i << " = " << bitset<n>(Part) << " = " << Part << endl; i++;
  }
  file_MCM_info << "##" << endl;

  file_MCM_info.close();
}

/*  
void PrintFile_MCM_Info(map<uint32_t, unsigned int > Kset, unsigned int N, map<uint32_t, uint32_t> MCM_Partition)
{
  uint32_t Part = 0, m=0;
  double C_param=0, C_geom=0;
  Complexity_MCM(MCM_Partition, N, &C_param, &C_geom);
  double LogL = LogL_MCM(Kset, MCM_Partition, N);

  cout << "********** General Information about the MCM: **********" << endl; 
  cout << "The chosen MCM has " << MCM_Partition.size() << " partitions and the following properties:" << endl;
  cout << "\t LogL = " << LogL << endl;
  cout << " \t C_param = " << C_param << " \t \t C_geom = " << C_geom << endl;
  cout << " \t Total complexity = " << C_param + C_geom << endl;
  cout << " \t MDL = " << LogL - C_param - C_geom << endl;
  cout << "  \t LogE = " << LogE_MCM(Kset, MCM_Partition, N) << endl;

  cout << endl << "********** Information about each part of the MCM: **********";
  cout << endl << "\t (the total LogE of the model is the sum of the values for each part, which are reported in the table below)";
  cout << endl << "\t !! The first operator of the specified basis corresponds to the rightmost bit !!";
  cout << endl << "\t !! The last operator corresponds to the leftmost bit !!" << endl << endl;;
  cout << "## 1:Part_int \t 2:Part_binary \t 3:LogL \t 4:C_param \t 5:C_geom \t 6:C_tot \t 7:LogE" << endl;

  for (map<uint32_t, uint32_t>::iterator i = MCM_Partition.begin(); i != MCM_Partition.end(); i++)
  {    
    Part = (*i).second;
    m = bitset<n>(Part).count();  // rank of the part (i.e. rank of the SCM)
    C_param = ParamComplexity_SubCM(m, N);
    C_geom = GeomComplexity_SubCM(m);

    cout << " \t " << Part << " \t " << bitset<n>(Part) << " \t";
    cout << LogL_SubCM(Kset, Part, N) << " \t";
    cout << C_param << " \t " << C_geom << " \t" << C_param + C_geom << " \t";
    cout << LogE_SubCM(Kset, Part, N) << endl;
  }
  cout << endl;
}*/

/******************************************************************************/
/*************      Print the model probabilities in a file     ***************/
/******************************************************************************/
void PrintFile_StateProbabilites_OriginalBasis_Vect(vector<pair<uint32_t, unsigned int>> Nset_Vect, list<uint32_t> Basis, map<uint32_t, uint32_t> MCM_Partition, unsigned int N, string filename = "Result")
{
  // Compute all the state probabilities:
  map<uint32_t, Proba> P_all = P_s_Vect(Nset_Vect, Basis, MCM_Partition, N);

  double *Pk_D = (double *)malloc((n+1)*sizeof(double)); 
  double *Pk_MCM = (double *)malloc((n+1)*sizeof(double)); 
  unsigned int k = 0;

  for(k=0; k<=n; k++)
  {
    Pk_D[k] = 0;
    Pk_MCM[k] = 0;
  }

  string Ps_filename = filename + "_DataVSMCM_Ps.dat";
  string Pk_filename = filename + "_DataVSMCM_Pk.dat";

  cout << "--->> Print information about the MCM in the file: \'" << filename << "_MCM_info.dat\'" << endl;
  cout << "--->> Print the state probabilities P(s) in the file: \'" << Ps_filename << "\'" << endl;
  cout << "--->> Print the probability of a state with k \'+1\' bits: \'" << Pk_filename << "\'" << endl << endl;

  //***** Print info about the model -- Print Basis and MCM:  **************/
  PrintFile_MCM_Info(Basis, MCM_Partition, filename);

  //***** Print P(s):  *****************************************************/
  uint32_t s;
  fstream file_Ps((OUTPUT_directory + Ps_filename), ios::out);

  file_Ps << "## s = states in the original basis" << endl;
  file_Ps << "## sig = states in the chosen new basis (ideally the best basis)" << endl;
  file_Ps << "## The chosen \'sig\'-basis and the chosen MCM are printed in the file " << filename << "_MCM_info.dat" << endl;

  // Print P(s): 
  file_Ps << "## " << endl;
  file_Ps << "## 1:s \t 2:P_D(s) \t 3:P_MCM(s) \t 4:sig" << endl;

  for (map<uint32_t, Proba>::iterator it_P = P_all.begin(); it_P!=P_all.end(); ++it_P)
  {   
    s = it_P->first;
    file_Ps << bitset<n>(s) << "\t" <<  (it_P->second).P_D_s << "\t" << (it_P->second).P_MCM << "\t" << bitset<n>((it_P->second).sig) << endl;

    k = bitset<n>(s).count();
    Pk_D[k] += (it_P->second).P_D_s;      // P[k] in the data
    Pk_MCM[k] += (it_P->second).P_MCM;    // P[k] from the MCM
  }
  file_Ps.close();

  //***** Print P(k):   ***************************************************/
  fstream file_Pk((OUTPUT_directory + Pk_filename), ios::out);

  file_Pk << "## 1:k \t 2:P_D(k) \t 3:P_MCM(k)" << endl;

  for(k=0; k<=n; k++)
  {
    file_Pk << k << "\t" << Pk_D[k] << "\t" << Pk_MCM[k] << endl;
  }
  file_Pk.close();
}

