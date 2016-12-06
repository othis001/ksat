#include "smoosher.h"

int main(int argc, char *argv[])
{
  bool up = false, psimp = false, mem = false, matrix = false;

  if(argc < 2)
  {
    cerr << "\nUsage is: smoosher filename options\nwhere filename is a file containing an A.\n"
         << "\nOptions:\n -up    starts smooshing from the bottom.\n"
         << " -simp  displays the simpices remaining.\n"
         << " -matrix saves the resulting boundary matrices.\n"
         << " -mem   switches to program to memory efficient mode at the cost of speed.\n\n";
    exit(1);
  }

  for(int i = 2; i < argc; i++)
  {
    if(argv[i][1] == 'u') up = true;
    if(argv[i][2] == 'e') mem = true;
    if(argv[i][1] == 's') psimp = true;
    if(argv[i][2] == 'a') matrix = true;
  }
  if(mem) cerr << "\nRunning in Memory Efficient mode.\n\n";
  else cerr << "\nRunning in fast mode. Add the -mem parameter to reduce speed and save memory.\n\n";

  char *filename = argv[1];

  Smooshable *theComplex = new Smooshable(filename);
  theComplex->print();
  theComplex->smoosh(up, mem);
  theComplex->printLengths();
  if(psimp) theComplex->printSimps();
  if(matrix) theComplex->printMatrix();

  return(0);  
}
