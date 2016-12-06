#include "smoosher.h"

/* The constructor for our smooshable complex. */
Smooshable::Smooshable(char *thefilename) : Complex::Complex(thefilename)
{
  if(sizeof(unsigned int) != 4)
  {
    cerr << "Error: Please use 32 bit integers." << endl;
    exit(1);
  }
  if(sizeA > 32) // Each simplex is stored as a 32 bit integer. Bigger integers are impractical with this algorithm.
  {
    cerr << "Error: This version won't work on As with more than 32 hyperedges." << endl;
    exit(1);
  }

  filename = thefilename; // It is necessary to store the file name so we can name our matrix files.

  trashtable = new long unsigned int*[sizeA]; // This array stores arrays of flags to store if we smooshed a simplex.
  for(int i=0; i < sizeA; i++)
    trashtable[i] = NULL; 

  down_pos = sizeA; // This is the highest dimension with potentially nontrivial homology.
}

/* Destructor for Smooshable. */
Smooshable::~Smooshable()
{
  for(int i=0; i<sizeA; i++)
    if(trashtable[i] != NULL) delete [] trashtable[i];
  delete [] trashtable;

  delete [] filename;
}

/* This function does the smooshing. up specifies to start from the bottom and mem specifies to run in memory efficient mode. */
void Smooshable::smoosh(bool up, bool mem)
{
  int i, j, counter, top;
  unsigned int *Ci, *Cip1 = NULL;
  bool up_smooshed, down_smooshed, first = true, first2 = true;
  
  /* Here we find the highest nontrivial Ci. */
  for(i=sizeA-1; i>=0; i--)
    if(CiLength[i] != 0)
    {
      down_pos = i;
      break;
    }

  up_pos   = 1;
  
  trashtable[0] = maketrashrow(sizeA);
  setTrash(0, trashtable[0]); // trashes { 1 } in C0.

  top = down_pos;

  while(1)
  {
  up_smooshed = false; down_smooshed = false;


  if(!(first && up)){
  
    /* Now start the from above smooshing. */ 
    cerr << "\nStarting from above smooshing ... \n";

    /* This loop stops the smooshing if we fail to completely smoosh a Ci. */
   
    while(1)
    {
      if(down_pos < up_pos)
      {
        down_pos++;
        break;
      }
  
      if(down_pos != top && trashtable[down_pos+1] != NULL)
      {    
        counter = 0;
        for(j=0; j<CiLength[down_pos+1]; j++)
          if(getTrash(j, trashtable[down_pos+1]) == false) counter++;
   
        if(counter != 0) 
        {
          down_pos++; 
          break;
        }
        else
        {
          delete [] trashtable[down_pos+1];
          trashtable[down_pos+1] = NULL;
        }
      }

      cerr << "Doing C" << down_pos << " ... ";
    
      Ci = new unsigned int[ CiLength[down_pos-1] ];
      for(i=0; i< CiLength[down_pos-1]; i++)
        Ci[i] = 0;

      if(!mem)
      {
        Cip1 = new unsigned int[ CiLength[down_pos] ];
        for(i=0; i< CiLength[down_pos]; i++)
          Cip1[i] = 0;

        makeCi(Ci, Cip1, down_pos, 1); 
      }
      else
        makeCi(Ci, Cip1, down_pos-1, 1); 
    
      if(trashtable[down_pos  ] == NULL) trashtable[down_pos  ] = maketrashrow(CiLength[down_pos  ]);
      if(trashtable[down_pos-1] == NULL) trashtable[down_pos-1] = maketrashrow(CiLength[down_pos-1]); 
    
      smooshDown(Ci, Cip1, down_pos, &up_smooshed);  

      delete [] Ci;
      if(!mem) delete [] Cip1;    

      cerr << "done\n"; 
    
      down_pos--;        
    }
  } 
  /* Now start the bottom up smooshing. */ 
  cerr << "\nStarting from below smooshing ...\n";

  while(1)
  {
    if(up_pos > down_pos)
    {
      up_pos--;
      break;
    }
    
    if(up_pos !=1 && trashtable[up_pos-2] != NULL)
    {  
      counter = 0;
      for(j=0; j<CiLength[up_pos-2]; j++)
        if(getTrash(j, trashtable[up_pos-2]) == false) counter++;
   
      if(counter == 0)
      {
        delete [] trashtable[up_pos-2];
        trashtable[up_pos-2] = NULL;
      }
      else
      {
        up_pos--; 
        break;
      }
    } 
         
    cerr << "Doing C" << up_pos-1 << " ... ";
    
    Ci = new unsigned int[ CiLength[up_pos-1] ];
    for(i=0; i< CiLength[up_pos-1]; i++)
      Ci[i] = 0;

    if(!mem)
    {
      Cip1 = new unsigned int[ CiLength[up_pos] ];
      for(i=0; i< CiLength[up_pos]; i++)
        Cip1[i] = 0;

      makeCi(Ci, Cip1, up_pos, 1);
    }
    else
      makeCi(Ci, Cip1, up_pos-1, 1);
       
    if(trashtable[up_pos] == NULL) trashtable[up_pos] = maketrashrow(CiLength[up_pos]);
    
    smooshUp(Ci, Cip1, up_pos, &down_smooshed);
   
    delete [] Ci;
    if(!mem) delete [] Cip1;
    
    cerr << "done\n";
    
    up_pos++;
  }
  if(up_smooshed == false && down_smooshed == false) break;
  first = false; 
  }
  cerr << endl;
}

/* These manage the lists of simplex tags. Stored in bit form to save memory. */

long unsigned int *Smooshable::maketrashrow(int thesize)
{
  int size = (thesize / 32) + 1;
  long unsigned int *trashrow = new long unsigned int[size];
  
  for(int i=0; i<size; i++)
    trashrow[i] = 0;

  return trashrow;
}
bool Smooshable::getTrash(int pos, long unsigned int *trashtable)
{ 
  if(trashtable[pos / 32] & (1 << (pos % 32) ) )
    return true;
  else
    return false;
}
void Smooshable::setTrash(int pos, long unsigned int *trashtable)
{
  trashtable[pos / 32] = trashtable[pos / 32] | (1 << (pos % 32) );
}

/* 
This performs the top down smooshing on Ci and Ci-1. The algorithm works as follows:
First a pass is made through Ci to count how many times each top dimensional face occurs. 
Now the program goes through Ci again and sees if any top dimensional faces occur only once. 
If it finds one that does, the program smooshes that face and it's containing simplex in Ci
and the count of the number of simplices is adjusted accordingly. This process continues 
until there are no remaining top dimensional faces in Ci that occur only once.
*/
  
void Smooshable::smooshDown(unsigned int *Ci, unsigned int *Cip1, int i, bool *smooshed)
{
    int j, numone;  

    char *facenums = new char[CiLength[i-1]];
    for(j=0; j < CiLength[i-1]; j++)
      facenums[j] = 0;

    qsmoosh(Ci, Cip1, i, 3, NULL, facenums, 0, NULL, NULL);

    while(1)
    {
      numone = 0;
      for(j=0; j < CiLength[i-1]; j++)
        if(facenums[j] == 1) numone++;
   
      if(numone == 0) break;

      qsmoosh(Ci, Cip1, i, 4, NULL, facenums, numone, NULL, NULL);
    }
    delete [] facenums; 
}

/* 
This performs the bottom up smooshing on Ci and Ci-1. This is a lot simplier than the above algorithm.
The program goes through each simplex in Ci and checks to see if all its top dimensional faces except
for one have been smooshed. If they have, the program eliminates that simplex and its one remaining
top dimensional face. The process keeps iterating until no more simplices may be cancelled. 
*/

void Smooshable::smooshUp(unsigned int *Ci, unsigned int *Cip1, int i, bool *smooshed)
{
   bool trouble;
   long unsigned int *smooshlist;
   smooshlist = maketrashrow(CiLength[i]);
 
    do
    {
      trouble = false;

      qsmoosh(Ci, Cip1, i, 2, &trouble, NULL, 0, smooshlist, NULL);
      
      if(trouble == true) *smooshed = true;
      
    }while(trouble == true);
   
    delete [] smooshlist; 
}

/* 
This function just produces Ci when given Ci-1 to speed up the process of smooshing. Ci and Ci-1 are
not stored at once to conserve memory (unless in fast mode). Operations are performed as Ci is generated.
*/

void Smooshable::qsmoosh(unsigned int *Ci, unsigned int *Cip1, int thei, int mode, bool *trouble, char *facenums, int numone, long unsigned int *smooshlist, unsigned int *facelist)
{
  int j, k, m, n, w, last;

  unsigned int catemp1;
  unsigned int face;
  int thepos, epos = 0, pos, pos1, counter, max;

  long unsigned int caint[width];
  long unsigned int caint2[width];
  bool intersPrev;

  if(Cip1) max = CiLength[thei];
  else max = CiLength[thei-1];


  for(j=0; j< max; j++)
  {
    if(Cip1) /* Fast mode */
    {
      /* This skips some simplices if we can. */
      if(~(trashtable[thei][j / 32]) == 0)
      {
        j = (j / 32) * 32 + 32 - 1;
        epos = (j / 32) * 32 + 32;
        continue;
      } 
      
      goto skip; // I know goto statements are taboo, but since I was just writing this code for myself, I don't care.
    }

    /* The following code is just to generate Ci from Ci-1. Note it is skipped if we are 
       in fast mode as seen by the goto statement above. */
    for(w=0; w<width; w++)
      caint[w] = ~0;
 
    for(k=0; k<32; k++)
      if((Ci[j] & (1 << (31 - k))) != 0)
      {
        last = k;
        for(w=0; w<width; w++)
          caint[w] = caint[w] & B[k*width + w];
      }
    
    for(k= last+1; k<sizeA; k++)
    {
      intersPrev = false;

      for(w=0; w<width; w++)
        caint2[w] = caint[w] & B[k*width + w]; 

      for(w=0; w<width; w++)
        if(caint2[w] != 0)
        {
          intersPrev = true;
          break;
        }

      if(intersPrev)
      {
        catemp1 = Ci[j];
        catemp1 = catemp1 | (1 << (31 -k));
        skip:
        /* This tests to see if a simplex has been smooshed already or had all its top
           dimensional faces smooshed. */
        if(!getTrash(epos, trashtable[thei]))
         {
          if(Cip1) catemp1 = Cip1[j];

          if(mode==2  && !getTrash(epos, smooshlist)) /* This mode is for the bottom up smooshing. */
          {
            counter = 0;
      
            for(m=31; m>=0; m--)
              if((catemp1 & (1 << m)) != 0)
              {
                face = catemp1 & (~(1 << m));

                pos = findinArray(face, Ci, CiLength[thei-1]); 
     
                if(getTrash(pos, trashtable[thei-1]) == false)
                {
                  counter++;
                  if(counter > 1) break;
                  thepos = pos;
                }
              }
   
            if(counter == 1) /* All but one of the top dimensional faces are smooshed. */
            {
              setTrash(thepos, trashtable[thei-1]);
              setTrash(epos, trashtable[thei]);
              *trouble = true;
            }
            else if(counter == 0)
              setTrash(epos, smooshlist);
          }
          else if(mode == 3) /* Here we count the top dimensional faces. */
          { 
            for(m=31; m>= 0; m--)
              if((catemp1 & (1 << m)) != 0)
              {
                face = catemp1 & (~(1 << m));

                pos = findinArray(face, Ci, CiLength[thei-1]);

                if(getTrash(pos, trashtable[thei-1]) == false)
                  facenums[pos]++;
              }
          }
          else if(mode == 4) /* Top down smooshing. */
          {
            for(m=31; m>= 0; m--)
              if((catemp1 & (1 << m)) != 0)
              {
                face = catemp1 & (~(1 << m));
                pos = findinArray(face, Ci, CiLength[thei-1]);

                if(facenums[pos] == 1)
                {
                  setTrash(pos, trashtable[thei-1]);
                  setTrash(epos, trashtable[thei]);
               
                  for(n=31; n>=0; n--) /* Decrements the face counters. */
                    if((catemp1 & (1 << n)) != 0)
                    {
                      face = catemp1 & (~(1 << n));
                      pos1 = findinArray(face, Ci, CiLength[thei-1]); 
                      facenums[pos1]--;
                    }
                  break;
                }
              }
          }
        }
        epos++; /* This is the counter for Ci. */
        if(Cip1) break;
      }
    }
  }
}

/* Uses binary search to find a simplex. */
int Smooshable::findinArray(unsigned int face, unsigned int *Ci, int length)
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

/* Here we find the dimensions of the chain groups after the smooshing. */
void Smooshable::printLengths()
{
  int i, j; int counter;
  
  for(i=0; i <= down_pos; i++)
  {
    if(i>up_pos && i<down_pos-1) // This is for simplices that may be out of reach of the smoosher. 
    {                            // Each of these Ci have no trashtable initialized but they are still all there.
      counter =  CiLength[i];
    }
    else
    {
      counter = 0;
      if(trashtable[i] != NULL)
        for(j=0; j< CiLength[i]; j++)
          if(getTrash(j, trashtable[i]) == false) counter++;
   
      if(i==0) counter++;
    }

    cerr << "C" << i << " dimension is " << counter << endl;
  }
}

/* This function saves the boundary matrices of any dimensions the smoosher tried and failed to completely smoosh. */
void Smooshable::printMatrix()
{ 
  int i, j, k, m, counter1, counter2, pos, alt, counter3, test, counter4, counter5;
  unsigned int *Ci, *Cip1, *nosmooshed, *supersmooshed, catemp1, face;

  char *thefilename;

  ofstream myfile; 


  long unsigned int *smooshlist;
  
  for(i=0; i< down_pos; i++)
  {
     if(trashtable[i] == NULL || trashtable[i+1] == NULL) continue;
     if(i>up_pos && i<down_pos-1) continue; // Skips generating matrices for dimensions the smoosher can't reach.

     counter1 = 0; counter2 = 0;
     for(j=0; j<CiLength[i]; j++)
       if(getTrash(j, trashtable[i]) == false) counter1++;
     for(j=0; j<CiLength[i+1]; j++)
       if(getTrash(j, trashtable[i+1]) == false) counter2++;
     if(counter1 == 0 || counter2 == 0) continue;

     counter3 = 0;

     Ci = new unsigned int[ CiLength[i] ];
     for(j=0; j< CiLength[i]; j++)
       Ci[j] = 0;
     makeCi(Ci,   NULL, i  , 1); 
     for(counter4=0, j=0; j< CiLength[i]; j++)
       if(getTrash(j, trashtable[i]) == false) counter4++;
     nosmooshed = new unsigned int[ counter4 ];
     for(k=0, j=0; j< CiLength[i]; j++)
       if(getTrash(j, trashtable[i]) == false) nosmooshed[k++] = Ci[j];
     delete [] Ci;
       
     Cip1 = new unsigned int[ CiLength[i+1] ];
     for(j=0; j< CiLength[i+1]; j++)
       Cip1[j] = 0;
     makeCi(Cip1, NULL, i+1, 1);  

     /* I know this is really ugly but the version of C++ I have installed won't allow most string manipulations. */  
     thefilename = new char[100]; k = 0;
     for(int l=0; filename[l] != '\0'; l++)
       thefilename[k++] = filename[l];

     thefilename[k++] = '_';

     thefilename[k++] = 'C';
     if(i>8)
     {
       thefilename[k++] = ((i+1) / 10) + 48;
       thefilename[k++] = ((i+1) % 10) + 48;
     }
     else thefilename[k++] = (i+1) + 48;
     thefilename[k++] = 'C';
     if(i>9)
     {
       thefilename[k++] = (i / 10) + 48;
       thefilename[k++] = (i % 10) + 48;
     }
     else thefilename[k++] = i + 48;
     thefilename[k++] = '.';
     thefilename[k++] = 't';
     thefilename[k++] = 'x';
     thefilename[k++] = 't';
     thefilename[k++] = '\0';

     myfile.open(thefilename);
     if(!myfile.is_open())
     {
       cerr << "\nFile I/O Error!\n";
       exit(1); 
     }

     smooshlist = maketrashrow(counter4);
     counter3 = 0;
     for(j=0; j< CiLength[i+1]; j++)
       if(getTrash(j, trashtable[i+1]) == false)
       {
         test = 0;
         catemp1 = Cip1[j];
         for(m=31; m>= 0; m--)
           if((catemp1 & (1 << m)) != 0)
           {
             face = catemp1 & (~(1 << m));
             pos = findinArray(face, nosmooshed, counter4);

             if(pos != -1)
             {
               test = 1;
               setTrash(pos, smooshlist);
             }
           }
         if(test == 1) counter3++;
       }
     
     for(counter5=0, j=0; j< counter4; j++)
       if(getTrash(j, smooshlist) == true) counter5++;
     supersmooshed = new unsigned int[ counter5 ];
     for(k=0, j=0; j< counter4; j++)
       if(getTrash(j, smooshlist) == true) supersmooshed[k++] = nosmooshed[j];
     delete [] nosmooshed;
     delete [] smooshlist;
      

     myfile << counter3 << " " << counter5 << " M" << endl;

     counter3 = 0;
     for(j=0; j< CiLength[i+1]; j++)
     {
       test = 0; alt = -1;
       if(getTrash(j, trashtable[i+1]) == true) continue;

       catemp1 = Cip1[j];
       for(m=31; m>= 0; m--)
         if((catemp1 & (1 << m)) != 0)
         {
           alt = alt * -1;
           face = catemp1 & (~(1 << m));
           pos = findinArray(face, supersmooshed, counter5);

           if(pos == -1) continue;

           myfile << counter3+1 << " " << pos+1 << " " << alt << endl;
           test = 1;
         }
       if(test == 1) counter3++;
     }
     myfile << "0 0 0";
    
     delete [] supersmooshed;
     delete [] Cip1;

     myfile.close();
  
     cerr << "Saved C" << i+1 << "C" << i << "matrix to file " << thefilename << endl;
     delete [] thefilename; 
  }


} 

/* Prints the Simplices. */
void Smooshable::printSimps()
{
  int counter, i, j, k, l;
  char *p;
  unsigned int *Ci;
  
  if(trashtable[0] == NULL)
  {
    trashtable[0] = maketrashrow(CiLength[0]);
    for(i=1; i<CiLength[0]; i++)
      setTrash(i, trashtable[0]);
  }
  else
    trashtable[0][0] = trashtable[0][0] & ~1;

  for(i=0; i<= down_pos; i++)
  {
    counter = 0;
    
    if(trashtable[i] != NULL)  
      for(j=0; j<CiLength[i]; j++)
        if(getTrash(j, trashtable[i]) == false) counter++; 

    cout << "C" << i << " is: ";  
   
    if(counter != 0 )  
    {
      Ci = new unsigned int[ CiLength[i] ];
      for(j=0; j< CiLength[i]; j++)
      Ci[j] = 0;
      
      makeCi(Ci, NULL, i, 1); 

      for(j=0; j<CiLength[i]; j++)
      {
        if(getTrash(j, trashtable[i]) == false)
        {
          cout << "{";
          l = 0;
          for(k=31; k>=0; k--)
            if((Ci[j] & (1 << k)) != 0)
            {
              cout << (31 - k) + 1;
              if(l != i) cout << " ";
              l++; 
            }
          cout << "} ";
        }    
      }
      delete [] Ci;
    }
    cout << endl;
  } 
}
