#include <map>
#include <fstream>
#include <bitset>

#include "data.h"

/******************************************************************************/
/**************** Log-likelihood (LogL), Geometric Complexity *****************/
/*************************  and Log-evidence (LogE) ***************************/
/******************************************************************************/
double LogE_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false);
double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, double *C_param, double *C_geom);

/********************************************************************/
/*************    CHECK if "Partition" IS A PARTITION   *************/
/********************************************************************/
//check if *Partition* is an actual partition of r basis elements, 
// i.e., that no basis element appears in more than 1 part of the partition.
// i.e., that each basis element only appears in a single part of the partition.

bool check_partition(map<uint32_t, uint32_t> Partition);

/********************************************************************/
/**********************    PRINT PARTITION   ************************/
/********************************************************************/
void Print_Partition(uint32_t *a)
{
  for (int i=0; i<n; i++)
  {    cout << a[i];  }
}

void Print_Partition_Converted(map<uint32_t, uint32_t>  partition)
{
  for (map<uint32_t, uint32_t>::iterator i = partition.begin(); i != partition.end(); i++)
  {    cout << (*i).second << " = " << bitset<n>((*i).second) << "\n";  }
  cout << endl;
}

/********************************************************************/
/****************    CONVERSION  of a partition    ******************/
/**********************   SPECIFIC TO MCM    ************************/
/********************************************************************/
// *** map<uint32_t, uint32_t>   --> .first = i = index    --> .second = a[i] = number of element in the part
map<uint32_t, uint32_t> Convert_Partition_forMCM(uint32_t *a, unsigned int r=n)
{
  map<uint32_t, uint32_t> Partition;
  uint32_t element = 1;

  for (int i=r-1; i>=0; i--)  // read element from last to first
    {    
      Partition[(a[i])] += element;
      element = element << 1;      //cout << a[i] << "\t :";   Print_Partition_Converted(Partition);
    }

//  cout << "Convert, " << Partition.size() << " parts: \t " ;
//  bool ok = check_partition(Partition);
//  cout << " \t --> Partition Ok? " << ok << endl << endl;

  return Partition;
}

// *** map<uint32_t, uint32_t> --> .first = i = index of the partition    --> .second = a[i] = number of element in the part
map<uint32_t, uint32_t> Convert_Partition_forMCM_withSubPart(uint32_t *a, bool *keep_SubPartition, unsigned int r=n)
{
  map<uint32_t, uint32_t> Partition;

  uint32_t element = 1;
  bool switch_ = false;
  *keep_SubPartition = true;

  for (int i=r-1; i>=0; i--)  // read element from last to first
    {    
      Partition[(a[i])] += element;  // cout << a[i] << "\t ";
      element = element << 1;
      if(switch_ == true && a[i] != 0)  { *keep_SubPartition = false;  }
      else if(a[i] == 0) { switch_ = true;  }
    }

//  cout << "Convert, " << Partition.size() << " parts: \t " ;
//  bool ok = check_partition(Partition);
//  cout << " \t --> Partition Ok? " << ok << endl << endl;

  return Partition;
}

/******************************************************************************/
/*********************  Compute all Partitions of a set   *********************/
/***************************   with Algorithm H   *****************************/
/******************************************************************************/
// *** find the first index j (from the right) such that a[j] != b[j]
int find_j(uint32_t *a, uint32_t *b, unsigned int r)
{
  int j = r-2;
  while (a[j] == b[j])  {   j--;  }
  return j;
}
/******************************************************************************/
// *** Version 1: 
// ***            Compare all the MCM of rank r, 
// ***            based on the r first elements of the basis used to build Kset:
/******************************************************************************/
map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false)
{
  cout << "--->> Search for the best MCM.." << endl << endl;

  int counter = 0, i = 0;
  string xx_st = "";
  for(int i=0; i<n-r; i++)
    {  xx_st += "_";  }

  // *** Print in file Best MCMs:
  fstream file_BestMCM((OUTPUT_directory + "BestMCM_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;

  // *** Print in file all MCMs:
  fstream file_MCM_Rank_r((OUTPUT_directory + "AllMCMs_Rank_r" + to_string(r) + ".dat").c_str(), ios::out);
  if(print_bool)
  {
    cout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
    cout << (OUTPUT_directory + "AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;
    file_MCM_Rank_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_MCM_Rank_r << "To activate the prints for all the MCMs of rank r="<< r << ","<< endl;
    file_MCM_Rank_r << " specify `print_bool=true` in the last argument of the function MCM_GivenRank_r();"; 
  }

  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;

  // *** LogE and Complexity
  double LogE = 0;
  double C_param = 0, C_geom = 0;
  map<uint32_t, uint32_t> Partition;

  // *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(n*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }

  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM(a, r), N);

  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_MCM_Rank_r << counter << ": \t";
    Partition = Convert_Partition_forMCM(a, r);
    LogE = LogE_MCM(Kset, Partition, N);     //LogE

    // *** Print in file:
    if(print_bool)
    {
      file_MCM_Rank_r << xx_st;
      for (i=0; i<r; i++)   {    file_MCM_Rank_r << a[i];  }     //Print_Partition(a);

      Complexity_MCM(Partition, N, &C_param, &C_geom);    //Complexity
      file_MCM_Rank_r << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter << endl;
    }

    // *** Best MCM LogE:
    if ( LogE > (*LogE_best)) 
    { 
      *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];  aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New \t " << counter << endl;  
    }
    else if ( LogE == (*LogE_best) )
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];  aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem \t " << counter << endl;    
    }

    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[r-1] up to reaching b[r-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[r-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }

  file_BestMCM.close();
  file_MCM_Rank_r.close();

  cout << "--> Number of MCModels (of rank r=" << r << ") that were compared: " << counter << endl;

  cout << endl << "********** Best MCM: **********";
  cout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  cout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;

  cout << "\t >> Best Model = ";
  cout << xx_st;
  for(int i=0; i<r; i++) {  cout << aBest[i];  }
  cout << "\t \t LogE = " << (*LogE_best) << endl << endl;

  Partition = Convert_Partition_forMCM(aBest, r);
  free(a); free(b); free(aBest);

  return Partition;
}
/******************************************************************************/
// *** Version 2:  
// ***            Compare all the MCM 
// ***            based on the k first elements of the basis used to build Kset
// ***            for all k=1 to r, where r <= basis.size()
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
/******************************************************************************/
map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_Ordered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false)
{
  int counter = 0, i = 0;
  int counter_subMCM = 0;

  string xx_st = "";
  for(int i=0; i<n-r; i++)
    {  xx_st += "_";  }

  // *** Print in file Best MCMs:
  fstream file_BestMCM(OUTPUT_directory + "BestMCM_Rank_r<=" + to_string(r) + "_Ordered.dat", ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;

  // *** Print in file all models:
  fstream file_allMCM_r((OUTPUT_directory +"AllMCMs_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  fstream file_allSubMCM((OUTPUT_directory +"AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat").c_str(), ios::out);

  if(print_bool)
  {
    cout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
    cout << (OUTPUT_directory +"AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;

    cout << "--> Print the LogE-value of all the MCM of rank k<" << r << " in the file '";
    cout << (OUTPUT_directory +"AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat") << "'" << endl << endl;

    file_allMCM_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
    file_allSubMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_allMCM_r << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allMCM_r << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
    file_allSubMCM << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allSubMCM << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
  }

  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;

  // *** LogE and Complexity
  double LogE = 0;
  double C_param = 0, C_geom = 0;
  map<uint32_t, uint32_t> Partition;

  // *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(r*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }

  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM(a, r), N);


  // *** SubPartitions (rank < n):
  bool keep_SubPartition = false;

  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_allMCM_r << counter << ": \t";

    // *** Original Partition:
    Partition = Convert_Partition_forMCM_withSubPart(a, &keep_SubPartition, r);     //Print_Partition_Converted(Partition); 
    LogE = LogE_MCM(Kset, Partition, N);     //LogE

    // *** Print in file:
    if(print_bool)
    {
      file_allMCM_r << xx_st;
      for (i=0; i<r; i++)   {    file_allMCM_r << a[i];  }     //Print_Partition(a);

      Complexity_MCM(Partition, N, &C_param, &C_geom);    //Complexity
      file_allMCM_r << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter << endl;
    }

    // *** Best MCM LogE:
    if ( LogE > (*LogE_best)) 
    { 
      *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New \t " << counter << endl;  
    }
    else if ( LogE == (*LogE_best) )
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem \t " << counter << endl;    
    }

    // *** Sub-Partition:
    if (keep_SubPartition)
    {
      counter_subMCM++;

      Partition.erase(0); //Print_Partition_Converted(Partition); 
      LogE = LogE_MCM(Kset, Partition, N);     //LogE

      // *** Print in file:
      if(print_bool)
      { 
        file_allSubMCM << xx_st;
        for (i=0; i<r; i++)     //Print_Partition(a);
        {
          if (a[i] == 0 )  {  file_allSubMCM << "x";  } 
          else {  file_allSubMCM << (a[i]-1);  } 
        }

        Complexity_MCM(Partition, N, &C_param, &C_geom);    //Complexity
        file_allSubMCM << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter_subMCM << endl;
      }

      // *** Best MCM LogE:
      if ( LogE > (*LogE_best) )
      { 
        *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a); 
        file_BestMCM << xx_st; 
        for (i=0; i<r; i++)   
          {    
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
          }
        file_BestMCM << "\t " << LogE << " \t New" << endl;  
      }
      else if ( LogE == (*LogE_best) )
      {  
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)  
        {
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);   aBest[i] = (a[i]-1);  } 
          else {  file_BestMCM << "x";  aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
      }
    }

    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[n-1] up to reaching b[n-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[n-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }

  file_BestMCM.close();
  file_allMCM_r.close();
  file_allSubMCM.close();

  cout << "--> Number of MCM models (of rank <=" << r << ") that have been compared: " << counter + counter_subMCM << endl << endl;
 
  cout << endl << "********** Best MCM: **********";
  cout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  cout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  
  cout << "\t >> Best Model = ";
  cout << xx_st;
  for(int i=0; i<r; i++) {  if(aBest[i] != -1)  {cout << aBest[i];}   else {cout << "x";}   }
  cout << "\t \t LogE = " << (*LogE_best) << endl << endl;

  Partition = Convert_Partition_forMCM(aBest, r); 
  free(a); free(b); free(aBest);

  return Partition;
}

/******************************************************************************/
// *** Version 3:  
// ***            Compare all the MCMs based on any subset of k elements 
// ***            of the of r first elements of the basis used to build Kset
// ***            for all k=1 to r, where r <= basis.size() 
// *** By default: - r=n
// ***             - the function doesn't print the logE-values for all the tested MCMs. To activate --> print_bool = true 
/******************************************************************************/
map<uint32_t, uint32_t> MCM_AllRank_SmallerThan_r_nonOrdered(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best, unsigned int r=n, bool print_bool=false)
{
  cout << "All MCM based on all subsets of r operators among n chosen independent operators, r<=n: " << endl;

  int counter = 0, i = 0;
  int counter_subMCM = 0;

  string xx_st = "";
  for(int i=0; i<n-r; i++)
    {  xx_st += "_";  }

  // *** Print in file Best MCMs:
  fstream file_BestMCM(OUTPUT_directory + "BestMCM_Rank_r<=" + to_string(r) + "_NonOrdered.dat", ios::out);
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;

  // *** Print in file:
  fstream file_allMCM_r((OUTPUT_directory +"AllMCMs_Rank_r=" + to_string(r) + ".dat").c_str(), ios::out);
  fstream file_allSubMCM((OUTPUT_directory +"AllMCMs_Rank_r<" + to_string(r) + "_NonOrdered.dat").c_str(), ios::out);

  if(print_bool)
  {
    cout << "--> Print the LogE-value of all the MCM of rank r=" << r << " in the file '";
    cout << (OUTPUT_directory +"AllMCMs_Rank_r=" + to_string(r) + ".dat") << "'" << endl << endl;

    cout << "--> Print the LogE-value of all the MCM of rank k<" << r << " in the file '";
    cout << (OUTPUT_directory +"AllMCMs_Rank_r<" + to_string(r) + "_Ordered.dat") << "'" << endl << endl;

    file_allMCM_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
    file_allSubMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;
  }
  else 
  { 
    file_allMCM_r << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allMCM_r << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
    file_allSubMCM << "To activate the prints for all the MCMs of rank r<="<< r << ","<< endl;
    file_allSubMCM << " specify `print_bool=true` in the last argument of the function MCM_AllRank_SmallerThan_r_Ordered();"; 
  }

  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(r*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(r*sizeof(uint32_t));
  for (int i=0; i<r; i++)
  {    a[i]=0; b[i]=1;  }
  int j = r-1;

  // *** LogE and Complexity
  double LogE = 0;
  double C_param = 0, C_geom = 0;
  map<uint32_t, uint32_t> Partition, Partition_buffer;

  //  *** Save Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(r*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }
  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM(a, r), N);

  // *** for SubModels:
  uint32_t amax = 0, atest = 0;

  //SubPartitions (rank < n):

  //ALGO H:
  while(j != 0) // && counter < 200)
  {
    // *** H2: Visit:   ******
    counter++;

    // *** Partition:
    Partition = Convert_Partition_forMCM(a, r); 
    LogE = LogE_MCM(Kset, Partition, N);     //LogE
    Complexity_MCM(Partition, N, &C_param, &C_geom);    //Complexity

    // *** Print in file:
    if(print_bool)
    {
      file_allMCM_r << xx_st;
      for (i=0; i<r; i++)   {    file_allMCM_r << a[i];  }
      file_allMCM_r << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter << endl;
    }

    // *** Best MCM LogE:
    if ( LogE > (*LogE_best)) 
    { 
      *LogE_best = LogE; 
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New" << endl;  
    }
    else if ( LogE == (*LogE_best)) 
    {  
      file_BestMCM << xx_st;
      for (i=0; i<r; i++)   {    file_BestMCM << a[i];     aBest[i]=a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
    }

    // *** Find max value in a[]:  //amax=0;   for(i=0; i<n; i++)  {  if (a[i] > amax) { amax = a[i]; } }
    if ( a[r-1] == b[r-1] ) { amax = b[r-1]; } else { amax = b[r-1]-1; } 

    // *** Sub-Partition: ***************************** //
    for(atest=0; atest<=amax; atest++)
    {
      counter_subMCM++;

      // *** Partition:
      Partition_buffer = Partition;
      Partition_buffer.erase(atest);

      LogE = LogE_MCM(Kset, Partition_buffer, N);     //LogE
      Complexity_MCM(Partition_buffer, N, &C_param, &C_geom);    //Complexity

      // *** Print in file:
      if(print_bool)
      {
        file_allSubMCM << xx_st;
        for (i=0; i<r; i++) 
        {
          if (a[i] == atest )  {  file_allSubMCM << "x";  } 
          else {  file_allSubMCM << a[i];  } 
        }
        file_allSubMCM << " \t" << LogE << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom) << " \t" << counter_subMCM << endl;
      }

      // *** Best MCM LogE:
      if ( LogE > (*LogE_best)) 
      { 
        *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  

        file_BestMCM << xx_st;
        for (i=0; i<r; i++)   
        {    
          if (a[i] > atest )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else if (a[i] < atest )  {  file_BestMCM << a[i];    aBest[i] = a[i];  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t New" << endl;  
      }
      else if ( LogE == (*LogE_best)) 
      {  
        file_BestMCM << xx_st;
        for (i=0; i<r; i++)  
        {
          if (a[i] > atest )  {  file_BestMCM << (a[i]-1);    aBest[i] = (a[i]-1);  } 
          else if (a[i] < atest )  {  file_BestMCM << a[i];    aBest[i] = a[i];  } 
          else {  file_BestMCM << "x";   aBest[i] = -1;  } 
        }
        file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
      }
    }

    if(a[r-1] != b[r-1])  {  a[r-1] += 1;  }   // H3: increase a[n-1] up to reaching b[n-1]
    else
    {  
      j = find_j(a,b,r);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[r-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[n-1]
        while ( j < (r-1) )
        {
          a[j] = 0;
          b[j] = b[r-1]; // = m
          j++; 
        }
        a[r-1] = 0;
      }
    }
  }

  file_BestMCM.close();
  file_allMCM_r.close();
  file_allSubMCM.close();

  cout << "--> Number of MCM models (of rank <=" << r << ") that have been compared: " << counter + counter_subMCM << endl;
 
  cout << endl << "********** Best MCM: **********";
  cout << endl << "\t !! The first operator of the basis provided corresponds to the bit the most on the right !!";
  cout << endl << "\t !! The last operator corresponds to the bit the most on the left !!" << endl << endl;;
  

  cout << "\t >> Best Model = ";
  cout << xx_st;
  for(int i=0; i<n; i++) {  if(aBest[i] != -1)  {cout << aBest[i];}   else {cout << "x";}   }
  cout << "\t \t LogE = " << (*LogE_best) << endl << endl;

  Partition = Convert_Partition_forMCM(aBest, r); 
  free(a); free(b); free(aBest);

  return Partition;
}



