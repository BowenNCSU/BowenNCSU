#include <iostream>
#include <cstdlib> // for malloc

using namespace std;

class A{public:
  int i;
  A(int i_init = 0) : i(i_init) {}
  ~A() {cout << "Calling ~A()" << endl; }

  // user-defined C-style typecasting
  //    (it is not possible in C++ to provide a user definition
  //     for static_cast, dynamic_cast, reinterpret_cast and const_cast)
  operator double() {return (double)i;}

};// A

int main(int argc, char* argv[])
{A *a_ptr;
 A *b_ptr;

 // allocating an object in a user-supplied memory location
 void* buffer;
 buffer = malloc(100);
 a_ptr = new (buffer) A(5);
 A& a = *a_ptr; // declare an alias reference
 cout << "a.i: " << a.i << " (double)a: " << (double)a << endl;

 // explicit destructor call: remove *a_ptr from C++'s object
 // pool but do not deallocate the space that *a_ptr occupied.
 a_ptr->~A();

 // reuse that space
 b_ptr = new (buffer) A(6);
 cout << "b_ptr->i: " << b_ptr->i << endl;
 b_ptr->~A();
}
