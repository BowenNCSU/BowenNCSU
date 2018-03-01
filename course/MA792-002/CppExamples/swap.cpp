#include <iostream>
using namespace std;

void swap(int& i, int& j)
{int& t=i; // removing & properly swaps refs
 i=j; j=t;
}

int main(void)
{int a = 7, b = 6;
 cout << a << ", " << b << endl;
 swap(a,b);
 cout << a << ", " << b << endl;
}
