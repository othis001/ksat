#include "findHomology.h"

int main(int argc, char *argv[])
{
  if(argc != 2)
  {
    printf("Usage is: findHomology filename\nwhere filename is a file containing A.\n");
    exit(1);
  
  }
  char *filename = argv[1];

  Homologicial *theComplex = new Homologicial(filename);
  theComplex->calculateHomology();
  theComplex->printHomology();

  delete theComplex;
  return(0);  
}
