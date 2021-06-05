#include <map>
#include <fstream>

#include "data.h"

/******************************************************************************/
/**************** Log-likelihood (LogL), Geometric Complexity *****************/
/*************************  and Log-evidence (LogE) ***************************/
/******************************************************************************/
double LogL_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false);
double LogE_MCM(map<uint32_t, unsigned int > Kset, map<uint32_t, uint32_t> Partition, unsigned int N, bool print_bool = false);
void Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, double *C_param, double *C_geom);

double LogE_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, bool print_bool = false);
double LogL_SubCM(map<uint32_t, unsigned int > Kset, uint32_t Ai, unsigned int N, bool print_bool = false);
double GeomComplexity_SubCM(unsigned int m);
double ParamComplexity_SubCM(unsigned int m, unsigned int N);

/********************************************************************/
/*******    PRINT INFO on each PART of an MCM (= a partition)   *****/
/********************************************************************/
void PrintTerminal_MCM_Info(map<uint32_t, unsigned int > Kset, unsigned int N, map<uint32_t, uint32_t> MCM_Partition)
{
  uint32_t Part = 0, m=0;
  double C_param=0, C_geom=0;
  Complexity_MCM(MCM_Partition, N, &C_param, &C_geom);
  double LogL = LogL_MCM(Kset, MCM_Partition, N);

  cout << "General Model information:" << endl; 
  cout << "\t LogL = " << LogL << endl;
  cout << " \t C_param = " << C_param << " \t C_geom = " << C_geom << endl;
  cout << " \t Total complexity = " << C_param + C_geom << endl;
  cout << " \t MDL = " << LogL - C_param - C_geom << endl;
  cout << "  \t LogE = " << LogE_MCM(Kset, MCM_Partition, N) << endl;

  cout << endl << "Information about each part of the selected MCM:" << endl;
  cout << "## 1:Part_int \t 2:Part_binary \t 3:LogL \t 4:C_param \t 5:C_geom \t 6:C_tot \t 7:LogE" << endl;

  for (map<uint32_t, uint32_t>::iterator i = MCM_Partition.begin(); i != MCM_Partition.end(); i++)
  {    
    Part = (*i).second;
    m = bitset<n>(Part).count();  // rank of the part (i.e. rank of the SCM)
    C_param = ParamComplexity_SubCM(m, N);
    C_geom = GeomComplexity_SubCM(m);

    cout << Part << " \t " << bitset<n>(Part) << " \t";
    cout << LogL_SubCM(Kset, Part, N) << " \t";
    cout << C_param << " \t " << C_geom << " \t" << C_param + C_geom << " \t";
    cout << LogE_SubCM(Kset, Part, N) << endl;
  }
  cout << endl;
}

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
map<uint32_t, uint32_t> Convert_Partition_forMCM(uint32_t *a)
{
  map<uint32_t, uint32_t> Partition;
  uint32_t element = 1;

  for (int i=n-1; i>=0; i--)  // read element from last to first
    {    
      Partition[(a[i])] += element;
      element = element << 1;      //cout << a[i] << "\t :";   Print_Partition_Converted(Partition);
    }

//  cout << "Convert, " << Partition.size() << " parts: \t " ;
//  bool ok = check_partition(Partition);
//  cout << " \t --> Partition Ok? " << ok << endl << endl;

  return Partition;
}

map<uint32_t, uint32_t> Convert_Partition_forMCM_rank(uint32_t *a, unsigned int r)
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
map<uint32_t, uint32_t> Convert_Partition_forMCM_withSubPart(uint32_t *a, bool *keep_SubPartition)
{
  map<uint32_t, uint32_t> Partition;

  uint32_t element = 1;
  bool switch_ = false;
  *keep_SubPartition = true;

  for (int i=n-1; i>=0; i--)  // read element from last to first
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
map<uint32_t, uint32_t> MCM_GivenRank_r(map<uint32_t, unsigned int > Kset, unsigned int r, unsigned int N, double *LogE_best)
{
  int counter = 0, i = 0;
  string xx_st = "";
  for(int i=0; i<n-r; i++)
    {  xx_st += "x";  }

  //fstream file_MCM_Rank_r(("MCMs_Rank_r" + to_string(r) + ".dat").c_str(), ios::out);
  //file_MCM_Rank_r << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;

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

  // *** Best MCMs:
  uint32_t *aBest = (uint32_t *)malloc(n*sizeof(uint32_t));
  for(int i=0; i<r; i++) {  aBest[i]=a[i];  }

  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM_rank(a, r), N);
  fstream file_BestMCM((OUTPUT_directory + "MCMsBest_Rank_r" + to_string(r) + ".dat").c_str(), ios::out);
    //list<uint32_t *> Best_MCM;
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;

  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_MCM_Rank_r << counter << ": \t";
    //file_MCM_Rank_r << xx_st;
    //for (i=0; i<r; i++)   {    file_MCM_Rank_r << a[i];  }     //Print_Partition(a);

    Partition = Convert_Partition_forMCM_rank(a, r);

    LogE = LogE_MCM(Kset, Partition, N);     //LogE
    //file_MCM_Rank_r << " \t" << LogE;

    //Complexity_MCM(Partition, N, &C_param, &C_geom);    //Complexity
    //file_MCM_Rank_r << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom);
    //file_MCM_Rank_r << " \t" << counter << endl;

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
  //file_MCM_Rank_r.close();

  cout << "Number of MCM models (of rank r) that was compared: " << counter << endl << endl;
  cout << "\t >> Best Model = ";
  cout << xx_st;
  for(int i=0; i<r; i++) {  cout << aBest[i];  }
  cout << "\t \t LogE = " << (*LogE_best) << endl << endl;

  Partition = Convert_Partition_forMCM_rank(aBest, r);

  return Partition;
}
/******************************************************************************/
// *** Version 2:  
// ***            Compare all the MCM 
// ***            based on the r first elements of the basis used to build Kset
// ***            for all r=1 to basis.size()  
/******************************************************************************/
map<uint32_t, uint32_t> MCM_allRank(map<uint32_t, unsigned int > Kset, unsigned int N, double *LogE_best)
{
  int counter = 0, i = 0;
  int counter_subMCM = 0;
  fstream file_allMCM("All_MCMs_Rank_n.dat", ios::out);
  file_allMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;

  fstream file_allSubMCM("All_MCMs_AllRank.dat", ios::out);
  file_allSubMCM << "# 1:Partition \t 2:LogE \t 3:C_K \t 4:C_geom \t 5:C_tot \t 6:counter" << endl;

  // *** H1: Initialisation:
  uint32_t *a = (uint32_t *)malloc(n*sizeof(uint32_t));
  uint32_t *b = (uint32_t *)malloc(n*sizeof(uint32_t));
  for (int i=0; i<n; i++)
  {    a[i]=0; b[i]=1;  }
  int j = n-1;

  // *** LogE and Complexity
  double LogE = 0;
  double C_param = 0, C_geom = 0;
  map<uint32_t, uint32_t> Partition;

  // *** Best MCMs:
  *LogE_best = LogE_MCM(Kset, Convert_Partition_forMCM(a), N);
  fstream file_BestMCM(OUTPUT_directory + "Best_MCMs.dat", ios::out);
    //list<uint32_t *> Best_MCM;
  file_BestMCM << "# 1:Partition \t 2:LogE " << endl;

  // *** SubPartitions (rank < n):
  bool keep_SubPartition = false;

  // *** ALGO H:
  while(j != 0)
  {
    // *** H2: Visit:
    counter++;  //file_allMCM << counter << ": \t";
    for (i=0; i<n; i++)   {    file_allMCM << a[i];  }     //Print_Partition(a);

    Partition = Convert_Partition_forMCM_withSubPart(a, &keep_SubPartition);     //Print_Partition_Converted(Partition); 

    // *** Original Partition:
    LogE = LogE_MCM(Kset, Partition, N);     //LogE
    file_allMCM << " \t" << LogE;

    Complexity_MCM(Partition, N, &C_param, &C_geom);    //Complexity
    file_allMCM << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom);
    file_allMCM << " \t" << counter << endl;

    // *** Best MCM LogE:
    if ( LogE > (*LogE_best)) 
    { 
      *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
      for (i=0; i<n; i++)   {    file_BestMCM << a[i];  } 
      file_BestMCM << "\t " << LogE << " \t New" << endl;  
    }
    else if ( LogE == (*LogE_best) )
    {  
      for (i=0; i<n; i++)   {    file_BestMCM << a[i];  }
      file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
    }

    // *** Sub-Partition:
    if (keep_SubPartition)
    {
      for (i=0; i<n; i++) 
      {
        if (a[i] == 0 )  {  file_allSubMCM << "x";  } 
        else {  file_allSubMCM << (a[i]-1);  } 
      }

      counter_subMCM++;
      Partition.erase(0); //Print_Partition_Converted(Partition); 

      LogE = LogE_MCM(Kset, Partition, N);     //LogE
      file_allSubMCM << " \t" << LogE;

      Complexity_MCM(Partition, N, &C_param, &C_geom);    //Complexity
      file_allSubMCM << " \t" << C_param << " \t" << C_geom << " \t" << (C_param + C_geom);
      file_allSubMCM << " \t" << counter_subMCM << endl;

      // *** Best MCM LogE:
      if ( LogE > (*LogE_best) )
      { 
        *LogE_best = LogE;  //Best_MCM.clear(); Best_MCM.push_back(a);  
        for (i=0; i<n; i++)   
          {    
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);  } 
          else {  file_BestMCM << "x";  } 
          }
        file_BestMCM << "\t " << LogE << " \t New" << endl;  
      }
      else if ( LogE == (*LogE_best) )
      {  
        for (i=0; i<n; i++)  
        {
          if (a[i] != 0 )  {  file_BestMCM << (a[i]-1);  } 
          else {  file_BestMCM << "x";  } 
        }
        file_BestMCM << "\t " << LogE << " \t Idem" << endl;    
      }
    }

    if(a[n-1] != b[n-1])  {  a[n-1] += 1;  }   // H3: increase a[n-1] up to reaching b[n-1]
    else
    {  
      j = find_j(a,b,n);  //H4: find first index j (from the right) such that a[j] != b[j]
      if (j==0) { break;  }   //H5: Increase a[j] unless j=0 [Terminate]
      else 
      {
        a[j] += 1;
        b[n-1] = b[j] + ((a[j]==b[j])?1:0);  // m
        j++;      //H6: zero out a[j+1], ..., a[n-1]
        while ( j < (n-1) )
        {
          a[j] = 0;
          b[j] = b[n-1]; // = m
          j++; 
        }
        a[n-1] = 0;
      }
    }
  }

  file_BestMCM.close();
  file_allMCM.close();

  cout << "Number of MCM models (of rank <=n) that was compared: " << counter << endl << endl;
  cout << "\t >> Best Model = ";
  //cout << xx_st;
  for(int i=0; i<n; i++) {  cout << a[i];  }
  cout << "\t \t LogE = " << (*LogE_best) << endl << endl;

  return Partition;
}



