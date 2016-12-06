#include "findEuler.h"

int main(int argc, char *argv[])
{
  char *filename;

  if(argc != 2)
  {
    cout << "Usage is: findEuler filename\nwhere filename is a file containing A." << endl;
    exit(1);
  }

  filename = argv[1];
  Complex *theComplex = new Complex(filename);
  theComplex->print();
  cout << "The Euler characteristic is " << theComplex->getEuler() << "." << endl;

  return(0);  
}
