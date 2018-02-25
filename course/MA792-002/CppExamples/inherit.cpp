// File: C++Example/inherit.C
// demo of constructors, destructors, inheritance, etc.
#include <iostream>

using namespace std;

class base {
public:
   int size;
   int *data_arr; // gets allocated by constructor
   base (int n) : size(n)
     { cout << "base constructor with arg: " << n << endl;
       data_arr = new int[n];
     }
   ~base()
     { cout << "base destructor" << endl;
       delete [] data_arr;
     }
   base() : size(10)
     { cout << "base constructor without arg" << endl;
       data_arr = new int[10];
     }
   base(const base& b)
   // this constructor makes a new object as a copy of b
     { cout << "base copy constructor" << endl;
       size = b.size;
       data_arr = new int[size];
       for(int i=0; i<size; i++) data_arr[i] = b.data_arr[i];
     } 

   base& operator=(const base& b)
   // this operator "copies" into *this, which is an existing object
   // if not explicitly defined, each member is assigned
   // this causes 2 problems: a) this->data_arr and b.data_arr will point to
   //                            the same array (which may not be desired)
   //                         b) on destruction, the array is deleted twice
   //                            (which may cause a run-time error)
   //                         Note: if you really want to share objects, you must
   //                               keep track of a reference count before deleting them
     { cout << "base assignment operator: COPIES ARRAY" << endl;
       if (this != &b ) // guard against self assignment
          {size = b.size;
           delete [] data_arr;
           data_arr = new int[size];
           for(int i=0; i<size; i++) data_arr[i] = b.data_arr[i];
          }
        return *this;
     }
   friend ostream& operator<<(ostream&, const base&);
}; // end class base

ostream& operator<<(ostream& os, const base& b)
   {os << "write base obj to stream: size = " << b.size;
    return os;
   }

class item {
public:
   int x;
   item (int m) : x(m)
     { cout << "item constructor with arg: " << m << endl;
     }
   ~item()
     { cout << "item destructor" << endl;
     }
   item() : x(0)
     { cout << "item constructor without arg" << endl;
     }
   item(const item& it)
     { cout << "item copy constructor" << endl;
       x = it.x;
     } 
   item& operator=(const item& it)
     { cout << "item assignment operator" << endl;
       if(this != &it)
         x = it.x;
       return *this;
     }
   friend ostream& operator<<(ostream&, item&);
};

ostream& operator<<(ostream& os, const item& it)
   {os << "write item obj to stream: x = " << it.x;
    return os;
   }

class derive : public base {
public:
   item it;
   derive (int m, int n) : base(n), it(m)
     { cout << "derive constructor with args: "
            << m << ", " << n << endl;
     }
   ~derive()
     { cout << "derive destructor" << endl;
     }
   derive()
     { cout << "derive constructor without arg" << endl;
     }
   derive(const derive& d) : base(static_cast<const base&>(d)), it(d.it)
   // implements the default copy constructor
   // only written to print on its invocation
     { cout << "derive copy constructor" << endl;
     } 
   derive& operator=(const derive& d)
   // implements default assignment operator
   // only written to print on its invocation
     { cout << "derive assignment operator" << endl;
       if( this != &d)
         {static_cast<base&>(*this) = static_cast<const base&>(d);
          // alternatively:   base::operator=(static_cast<const base&>d);
          it = d.it;
         }
       return *this;
     }
   friend ostream& operator<<(ostream&, const derive&);
};

ostream& operator<<(ostream& os, const derive& d)
   {os << "write derive obj to stream: " << endl
       << "   " << static_cast<const base&>(d) << endl
       << "   " << d.it;
    return os;
   }
int main(int argc, char *argv[])
{ // 1. explicity create a derive object
  //    (this object must be explicity deleted)
  cout << "1.  new derive(5,20);" << endl;
  derive* d1_ptr;  d1_ptr = new derive(5,20);
  cout << *d1_ptr << endl << endl;
  
  // 2. create a derive object in an automatic variable
  //    (which will be deleted as soon as the var goes out of scope)
  cout << "2.  derive d2;" << endl;
  derive d2; cout << d2 << endl << endl;
  
  // 3. assign a derive object to another
  cout << "3.  d2 = *d1_ptr;" << endl;
  d2 = *d1_ptr; cout << d2 << endl << endl;
  
  // 4. create a derive object by copy through "initialization"
  cout << "4.  derive d3 = d2;" << endl;
  derive d3 = d2; // Note: same as  derive d3(d2);  [NOT operator=]
  cout << d3 << endl << endl;
  
  // 5. assign to base object an expression whose type matches
  //    the type of a constructor parameter
  cout << "5.  base b; b = 4;" << endl;
  base b1; b1 = 4; // Note: such an implicit conversion can be disallowed by
                   //       declaring in class base:  explicit base(int);
  cout << b1 << endl << endl;
  
  // 6. invoke the derive destructor
  cout << "6.  delete d1_ptr;" << endl;
  delete d1_ptr; cout << endl;

  // 7. various implict destructors will be invoked before returning
  cout << "The end" << endl;
  return 0;
}
