#include <iostream>
class A
{public:
 char C;
 A() : C('A'){}
 operator int() {return 5;}
};

int main(void)
{
A a;
cout << (int)a << endl;
}// end main
