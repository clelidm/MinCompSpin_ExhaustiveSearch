#include <bitset>
#include <cmath>       /* tgamma */
#include <map>

using namespace std;

#include "data.h"


/******************************************************************************/
/**************************   MODEL COMPLEXITY   ******************************/
/******************************************************************************/

/********* of a Sub-Complete Model based on m basis operators     *************/
/******************************************************************************/
// for 1 <= m <= n . Rem C(m=1) = log(pi)
double GeomComplexity_SubCM(unsigned int m)     // Geometric complexity
{
  double pow = (double) ( 1UL << (m-1) );
  return (log(M_PI) * pow - lgamma(pow)  );   // lgamma(x) = log(gamma(x))
}

double ParamComplexity_SubCM(unsigned int m, unsigned int N)  // Parameter Complexity
{
  uint32_t K = (1UL << m) - 1;  // number of interactions
  return K * log(((double) N)/2./M_PI) / 2.;
}

/********* of the Minimally Complex Model (MCM) defined by "Partition"   ******/
/******************************************************************************/
// Compute separately: -- the first order complexity    --> stored in C_param
//                     -- and the geometric complexity  --> stored in C_geom

double Complexity_MCM(map<uint32_t, uint32_t> Partition, unsigned int N, double *C_param, double *C_geom)
{
  *C_param = 0;   *C_geom = 0;
  uint32_t m_i = 0;  // number of elements in Ai

  for (map<uint32_t, uint32_t>::iterator Part = Partition.begin(); Part != Partition.end(); Part++)
  {
    m_i = bitset<n>((*Part).second).count();
    (*C_param) += ParamComplexity_SubCM(m_i, N);
    (*C_geom) += GeomComplexity_SubCM(m_i);
  }  

  return (*C_param) + (*C_geom);
}

