#include "findEuler.h"

/* The constructor for our complex. */
Complex::Complex(char *filename)
{
  int i, j, k, m, temp;
  ifstream inFile;
  
  inFile.open(filename);
  if (!inFile) {
    cerr << "Error: Unable to open file." << endl;
    exit(1); // Terminate with error.
  }

  /* From here on, we read in an A (our input file). */ 

  /* Here we find the number of subsets in A and the size of F. */ 
  sizeA = 0; sizeF = 0;
  while(inFile >> temp)
  {
    if(temp > sizeF) sizeF = temp;
    if(temp == 0) sizeA++; 
  }
  inFile.close();  

  /* Finally we put the subsets into B, an array of bits that stores A. */
  /* For example {1 4} would be 1 0 0 1.                                */
  inFile.open(filename);

  width = sizeF / 32;
  if(sizeF % 32 > 0) width++; // width is the number of integers dedicated to storing each subset.
  
  B = new long unsigned int[width * sizeA];
  for(i=0; i< width * sizeA; i++) B[i] = 0;

  k = 0;
  for(i=0; i < sizeA; i++)
  {
    m = 0;
    for(j=0; j< sizeF+1; j++)
    {
      inFile >> temp;

      if(temp == 0) break;
      if(temp > 32 * (m+1) ) m++;
      
      B[k+m] = B[k+m] | (1 << (temp - ( 32 * m ) - 1)); 
    }
    k= k + width;
  }
  inFile.close();

  CiLength = new int[sizeA];
  for(i=0;i<sizeA; i++)
    CiLength[i] = 0;

  makeCi(NULL, NULL, sizeA, 0); // Counts the simplices.
}

/* This function computes the i and i+1 simplices. It would be more efficient to compute all simplices, but
   due to memory constraints, we only handle two dimensions at a time.
   Ci   - the list of i dimensional simplices.
   Cip1 - the list of i+1 dimensional simplices.
   thei - specifies what i is.
   mode - if we are in mode 0 then we only compute and store dimensions of chain groups.
          if we are in mode 1 then we store the i and i+1 dimensional simplices.
*/
void Complex::makeCi(unsigned int *Ci, unsigned int *Cip1, int thei, int mode)
{
  /* This is the initalization stage. */
  int i; bool intersPrev;
   
  long unsigned int **a = new long unsigned int*[sizeA];
  for(i=0; i<sizeA; i++)
  {
    a[i] = NULL;
  }
  
  long unsigned int **aPointer = new long unsigned int*[sizeA];
  for(i=0; i<sizeA; i++)
  {
    aPointer[i] = NULL;
  }
  
  long unsigned int **inters = new long unsigned int*[sizeA];
  for(i=0; i<sizeA; i++)
  {
    inters[i] = NULL;
  }
   
  char catemp[thei+1];
  int catempcounter;
  int epos = 0;
  int counter, k, n, m, pos, thepos; char temp[thei];	
  
  int Cip = 0;
  int Cip1p = 0;
  
  int level;
  
  level = 0;
  aPointer[0] = B;
  
  /* This while loop will calculate chain groups of A. */
  while(1)
  {
  
    /* First we grab an element from aPointer. */
    a[level] = aPointer[level];
    /* Now increment aPointer. */
    aPointer[level] = aPointer[level] + width;
    
    intersPrev = false;  
    
    /* Here we calculate the intersection of a with all the previous as. */
    if(level > 0)
    {
      if(inters[level] != NULL) delete [] (inters[level]);
           
      inters[level] = intersection(inters[level-1],a[level]);
    }
    else
      inters[0] = a[0];  
    
    for(i=0; i<width; i++)
      if (inters[level][i])
      {
        intersPrev = true;
        break;
      }
    
    /* If a does intersect all the previous a's: */  
    if(intersPrev)
    { 
      if(mode == 1) // If mode=1, we are storing simplices.
      {
        if(Cip1 != NULL)
        {
          if(level == thei)
          {
	        for(i=0; i<=level; i++)
	        {
	          Cip1[Cip1p] |= (1 << (31 - ( (a[i] - B)/ width  ) ) );
	        }  
  
	        Cip1p++;
          }
          else if(level == thei-1)
          {
	        for(i=0; i<=level; i++)
	        {
	          Ci[Cip] |= (1 << (31 - ( (a[i] - B)/ width  ) ) );
	        }  
  
   	        Cip++;
          }
        }
        else if(level == thei)
        {
	      for(i=0; i<=level; i++)
	      {
	        Ci[Cip] |= (1 << ( 31 - ( (a[i] - B)/ width  ) ) );
	      }  
  
	      Cip++;
        } 
      }
      else if(mode == 0) CiLength[level]++; // If mode=0, we are just counting simplices.
 
      if((aPointer[level] -B < width*sizeA) && level < thei)
      {
        aPointer[level+1] = aPointer[level];
        level++; 
      }
    }
    
    /* If we have reached the end of A we go back a level. */
    if(aPointer[level]-B >= width*sizeA) 
      level--;    
      
    /* Finally we break if level goes below 0. */
    if(level < 0) break;
  }

  /* Finally we do some garbage collection. */
  for(i=1; i<sizeA; i++)
  {
    if(inters[i] != NULL) delete [] inters[i];
  }
  delete [] a;
  delete [] aPointer;
  delete [] inters;
}

/* Returns the Euler characteristic of our Complex. */
int Complex::getEuler()
{
  int i, X =0;

  for(i=0; i<sizeA; i++)
    if(i % 2 == 0) X += CiLength[i];
    else X -= CiLength[i];

  return X;
}

void Complex::print()
{
  int i, j, k;
  
  cout << "Got A = ";
  for(i=0; i<sizeA; i++)
  {

    for(j=0; j<width; j++)
    {
      cout << "{ ";
      for(k=0; k<32; k++)
      {
        if(B[i*width +j] & (1 << k))
          cout << k+1+j*32 << " ";
      }
    }
    cout << "}";
    if(i != sizeA-1) cout << " , ";
  }
  cout << endl << "size F = " << sizeF << " and size A = " << sizeA << endl; 

}
Complex::~Complex()
{
  delete B;
}

/* This function calculates the intersection of cells a and b. */
long unsigned int *Complex::intersection(long unsigned int *a, long unsigned int *b)
{
  int i;
  long unsigned int *inter;
  
  inter = new long unsigned int[width];
  
  for(i=0; i<width; i++)
    inter[i] = a[i] & b[i];  

  return(inter);
}
