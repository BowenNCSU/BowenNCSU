#include <iostream>
#include <typeinfo>

using namespace std;

template <class T>
class A
{public: T x;
 A(T y) : x(y) {}
 typedef T elemT;
};

template <> // template specialization
class A<float>
{public: float z; // new data member
         float x;
A<float>(float y) : z(y+10.0), x(y) {}
};

template< template<class> class C >
// template-template parameters
class B
{public:
 C<int> ci;
 C<float> cf;
 B(int a, float b) : ci(a), cf(b) {}
};

template<class T> // template function: overloaded <<
ostream& operator<<(ostream& out, const A<T>& a)
{out << "An A<" << typeid(T).name() << ">: x member: " << a.x << endl;
 return out; }
 

class A<double>; // put decl. in namespace scope
int main(void)
{// decl.  class A<double>;  not in namespace scope
 A<int> a(5);
 A<char> b('A');
 A<float> c(1.5);
 B<A> u(5, 1.5);

 cout << a << b << c << u.ci << u.cf << u.cf.z << endl;

}//main
