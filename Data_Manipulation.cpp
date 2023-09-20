#include <iostream>
#include <fstream>
#include <sstream>
#include <list>
#include <bitset>
#include <map>
#include <vector>

/********************************************************************/
/**************************    CONSTANTS    *************************/
/********************************************************************/
#include "data.h"

/******************************************************************************/
/**************************     READ FILE    **********************************/
/******************************************************************************/
/**************    READ DATA and STORE them in Nset    ************************/
vector<pair<uint32_t, unsigned int>> read_datafile_Vect(unsigned int *N, string filename = datafilename)    // O(N)  where N = data set size
{
  string line, line2;     uint32_t state = 0;
  (*N) = 0;            // N = dataset size
  cout << endl << "--->> Read \"" << filename << "\",\t Build Nset...";

// ***** data are store in Nset:  ********************************
  map<uint32_t, unsigned int> Nset; // Nset[mu] = #of time state mu appears in the data set
  
  ifstream myfile (filename.c_str());
  if (myfile.is_open())
  {
    while ( getline (myfile,line))
    {
      line2 = line.substr (0,n);          //take the n first characters of line
      state = bitset<n>(line2).to_ulong();   //convert string line2 into a binary integer
      Nset[state] += 1;
      //cout << line << endl;   //cout << state << " :  " << bitset<n>(state) << endl;
      (*N)++;
    }
    myfile.close();
  }
  else cout << "Unable to open file"; 

  cout << "\t\t data size N = " << (*N) << endl;

  // convert map to a vector
  vector<pair<uint32_t, unsigned int>> Nset_Vect(Nset.size());    //Nset_Vect.resize(Nset.size()); // reserve(n)

  int i=0;
  for (auto& my_pair : Nset)
  {
    Nset_Vect[i]=my_pair;
    i++;
  }

  return Nset_Vect;
}

/*
void get_data_vect(string fname, Partition_vect &p_struct) {

  // https://github.com/clelidm/MinCompSpin_Greedy
  map<__uint128_t, unsigned int> data;

  __uint128_t state;
  string fpath = "../input/data/" + fname + ".dat";
  string line, subline;
  ifstream myfile(fpath);

  while (getline(myfile, line)) {

    subline = line.substr(0, p_struct.n);
    state = string_to_int(subline, p_struct.n);
    data[state]++;
    p_struct.N++;

  }

  myfile.close();

  cout << "- loaded: " << fpath;
  cout << " (" << p_struct.N << " samples)\n" << endl;
  //cout << "test test N = " << data.size() << endl;

  (p_struct.data).resize(data.size()); // reserve(n)
    int i=0;
    for (auto& my_pair : data)
    {
      p_struct.data[i]=my_pair;
      i++;
    }
}*/

/******************************************************************************/
/*********************     CHANGE of BASIS: one datapoint  ********************/
/******************************************************************************/
// Given a choice of a model (defined by the m basis vector) --> return the new m-state (state in the new m-basis)
// Rem: must have m <= n 
uint32_t transform_mu_basis(uint32_t mu, list<uint32_t> basis)
{
  uint32_t bit_i = 1;
  uint32_t final_mu = 0;

  list<uint32_t>::iterator phi_i;

  for(phi_i = basis.begin(); phi_i != basis.end(); ++phi_i)
  {
    if ( (bitset<n>( (*phi_i) & mu ).count() % 2) == 1) // odd number of 1, i.e. sig_i = 1
      {
        final_mu += bit_i;
      }
    bit_i = (bit_i << 1);
  }

  return final_mu;
}

/******************************************************************************/
/************************** K_SET *********************************************/
/******************************************************************************/
// Build Kset for the states written in the basis of the m-chosen independent 
// operator on which the SC model is based:

vector<pair<uint32_t, unsigned int>> build_Kset_Vect(vector<pair<uint32_t, unsigned int>> Nset_Vect, list<uint32_t> Basis, bool print_bool=false)
// sig_m = sig in the new basis and cut on the m first spins 
// Kset[sig_m] = #of time state mu_m appears in the data set
{
  map<uint32_t, unsigned int> Kset;
  uint32_t sig_m;    // transformed state and to the m first spins

  cout << endl << "--->> Build Kset..." << endl;

//Build Kset:
  for (auto const& it : Nset_Vect)
  {
    sig_m = transform_mu_basis((it).first, Basis); // transform the initial state s=(it).first into the new basis
    Kset[sig_m] += ((it).second); // number of time state s appear in the dataset

    if (print_bool)  {  cout << ((it).first) << ": \t" << bitset<n>((it).first) << " \t" << sig_m << ": \t" << bitset<n>(sig_m) << endl; }
  }
  cout << endl;

  // convert map to a vector
  vector<pair<uint32_t, unsigned int>> Kset_Vect(Kset.size());

  int i=0;
  for (auto& my_pair : Kset)
  {
    Kset_Vect[i]=my_pair;
    i++;
  }

  return Kset_Vect;
}




