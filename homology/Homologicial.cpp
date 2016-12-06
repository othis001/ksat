#include "findHomology.h"

/* The constructor for our Homologicial complex. */
Homologicial::Homologicial(char *thefilename) : Complex::Complex(thefilename)
{
  if(sizeA > 32) // Each simplex is stored as a 32 bit integer. Bigger integers are impractical with this algorithm.
  {
    cerr << "Error: This version won't work on As with more than 32 hyperedges." << endl;
    exit(1);
  }

  homology = new int[sizeA];
  for(int i=0; i<sizeA; i++)
    homology[i] = 0;
}

/* Destructor for Homologicial. */
Homologicial::~Homologicial()
{
  delete [] homology;
}

/* Calculates the rational homology of our complex. */
void Homologicial::calculateHomology() 
{
  int i;
  unsigned int *Ci, *Cip1;

  int brank[sizeA]; // Array of boundary matrix ranks.
  for(int i=0; i<sizeA; i++)
    brank[i] = 0;

  int j, m, alt, pos; unsigned int catemp, face;
  PID_integer ZZ;
  SparseMatrix<PID_integer> *BMatrix;
  for(i=0; i< sizeA-1; i++)
  {
    Ci   = new unsigned int[ CiLength[i] ];
    Cip1 = new unsigned int[ CiLength[i+1] ];
    makeCi(Ci, Cip1, i, 1); 

    BMatrix = new SparseMatrix<PID_integer>(ZZ, CiLength[i+1], CiLength[i]);   

    for(j=0; j< CiLength[i+1]; j++) // This creates the boundary matrix for di+1:Ci+1->Ci.
    {
      alt = -1;
      catemp = Cip1[j];
      for(m=31; m>= 0; m--)
        if((catemp & (1 << m)) != 0)
        {
          alt = alt * -1;
          face = catemp & (~(1 << m));
          pos = findinArray(face, Ci, CiLength[i]);
          BMatrix->setEntry(j+1, pos+1, alt); // Linbox's matrices are indexed starting at 1 and not 0.
        }
    }
    delete [] Ci;
    delete [] Cip1;
     
    rank(brank[i+1], *BMatrix); // Here we use Linbox to find its rank.

    delete BMatrix;
  }

  for(i=0; i<sizeA-1; i++) // See my master's thesis (http://trace.tennessee.edu/utk_gradthes/228/) for an explanation of this formula.
    homology[i] = CiLength[i] - brank[i] - brank[i+1];
}

/* Prints the rational homology. */
void Homologicial::printHomology()
{
  for(int i=0; i < sizeA-1; i++)
    cout << "H" << i << " is Q^" << homology[i] << endl;
}

/* Uses binary search to find a simplex. */
int Homologicial::findinArray(unsigned int face, unsigned int *Ci, int length)
{
  int low, high, mid;
  
  low = 0;
  high = length-1;
  
  while(low <= high)
  {
    mid = (low+high)/2;
      
    if(Ci[mid] > face) low=mid+1;
    else if(Ci[mid] < face) high=mid-1;
    else return mid;
  }
  return -1;
}

