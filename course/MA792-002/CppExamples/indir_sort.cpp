// File: C++Examples/indir_sort.cpp
// Sorts array "indirectly" by sorting and index array
// Programmer: E. Kaltofen, Sep 28, 1998
//                          Apr 22, 2008 compliant static template members etc.
// STL items demonstrated
// generic algorithm sort with a comparator "function"
// /usr/lib/gcc/i686-pc-cygwin/3.4.4/include/c++/bits

#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

// must overload comparison on index array indir_comp I[]:
//                I[i] < I[j] if A[I[i]] < A[I[j]]
template <class ElType>
class indir_comp {ElType* Aptr;
           public:
           static int comp; // counts comparisons in sorting
           // overload function call on indir_comp object
           //     (needed in third argument to STL sort function)
           indir_comp(ElType* p) : Aptr(p) {}
           bool operator()(const int& i, const int& j)
                          {comp++; return Aptr[i] < Aptr[j];}
//           friend int indir_sort (ElType[], int, int[]);
          };

// declare static member
template <class ElType>
int indir_comp<ElType>::comp = 0;

template <class ElType>
int indir_sort(ElType A[], int n, int I[])
// set the array I[0],...,I[n-1] such that
// A[I[0]] <= A[I[1]] <= ... <= A[I[n-1]]
// note: one must have ElType& ElType::operator<(ElType&)
{for(int i=0; i<n; i++) I[i] = i;
 indir_comp<ElType> comp_fun(A); comp_fun.comp = 0;
 // stable_sort(I, I+n, comp_fun );
 sort(I, I+n, comp_fun );
 return comp_fun.comp;
}

#include <stdlib.h>
#include <stdio.h> // for sscanf

template // explicitly instantiate (not necessary)
int indir_sort<int>(int*, int, int*);

int main(int argc, char** argv)
// run as "a.out 100 25 11", where 100 is the size to be tested,
// and 25 is the range of entries;
// 11 is the seed for the random number generator;
{
  int *array, *indexarray; int comp, i, size, mod, seed;

  if (argc < 4) 
    {
      cerr << "Call as \"" << argv[0] << " <size> <mod> <seed> \"" << endl; 
      exit(1);
    };

  sscanf(argv[1], "%d", &size); sscanf(argv[2], "%d", &mod);
  array = new int[size];
  indexarray = new int[size];

  sscanf(argv[3], "%d",&seed);
  srand(seed);
  for(i = 0; i < size; i++) array[i] = rand() % mod;
  cout << endl;

  cout << "Input array" << endl;
  for(i=0;i<size;i++) cout << setw(3) << array[i] << " ";
  cout<<endl;
  // print positions also
  for(i=0;i<size;i++) cout << setw(3) << i << " "; cout<<endl;

  comp = indir_sort<int>(array, size, indexarray);

  cout << "STL-sorted array after " << comp << " comparisons" << endl;
  for(i = 0; i < size; i++) cout << setw(3) << array[indexarray[i]] << " ";
  cout << endl;
  for(i=0;i<size;i++) cout << setw(3) << indexarray[i] << " "; cout<<endl;
}
