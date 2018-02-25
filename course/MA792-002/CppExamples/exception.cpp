// File: exception.cpp
//       runs under g++ version 2.8.1 ("add gcc281")

#include <iostream>
#include <stdexcept>

using namespace std;

// class header info
class Integer
  {int x;
   public:
   Integer() : x(0) {}
   Integer(const int i) : x(i) {}
   inline Integer operator/(const Integer j)
      throw(std::invalid_argument); // indicate what the fun may throw
                         // this is not necessary but helps the caller
   friend ostream& operator<<(ostream& os, const Integer& i);
  };
// implementation
inline Integer Integer::operator/(const Integer j)
      throw(std::invalid_argument)
      {if(j.x == 0)
         throw(std::invalid_argument(
                 "in Integer::operator/: division by zero")
              );
       else return this->x/j.x;
      }
ostream& operator<<(ostream& os, const Integer& i)
  {os << i.x; return os; }
int main(void)
{int x=1, y=0;
 // cout << x/y << endl;
 try
 { Integer three = 3; Integer two = 2; Integer zero = 0;
   cout << three/two << endl;
   cout << three/zero << endl;
 }
 catch (std::invalid_argument& m)
     { cout << "Exception: " << m.what() << endl;
     }
}
