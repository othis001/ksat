#include <iostream>
#include <fstream>
#include <cstdlib>

/* The needed Linbox libraries. */
#include "linbox/field/gf2.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/solutions/rank.h"

using namespace LinBox;
using namespace std;

class Complex
{
  public:
    Complex(char *);
    void print();
    int getEuler();
    ~Complex();
  protected:
    long int sizeF, sizeA, width;
    int *CiLength;
    long unsigned int *B;
    void makeCi(unsigned int *, unsigned int *, int , int );
    long unsigned int *intersection(long unsigned int *, long unsigned int *);
};

class Homologicial : public Complex
{
  public:
    Homologicial(char *);
    ~Homologicial();
    void calculateHomology();
    void printHomology();
  protected:
    int *homology;
    int findinArray(unsigned int, unsigned int *, int);
};

