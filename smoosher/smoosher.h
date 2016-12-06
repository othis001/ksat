#include <iostream>
#include <fstream>
#include <cstdlib>

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

class Smooshable : public Complex
{
  public:
    Smooshable(char *);
    void smoosh(bool, bool);
    void printLengths();
    void printMatrix();
    void printSimps();
    ~Smooshable();
  protected:
    char *filename;
    long unsigned int **trashtable;
    int down_pos, up_pos;
    long unsigned int *maketrashrow(int);
    bool getTrash(int, long unsigned int *);
    void setTrash(int, long unsigned int *);
    int findinArray(unsigned int, unsigned int *, int);
  private: // No one else should ever see these functions.
    void smooshDown(unsigned int *, unsigned int *, int , bool *);
    void smooshUp(unsigned int *, unsigned int *, int , bool *);
    void qsmoosh(unsigned int *, unsigned int *, int , int , bool *, char *, int , long unsigned int *, unsigned int *);
};

