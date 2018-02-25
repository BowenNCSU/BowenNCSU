// Co-routines in C++ using pointers to members
// See main() first
// Programmer: Dmitriy Morozov (foxcub@foxcub.org)

#include <iostream>

using namespace std;

class A
{
    private:
	//typedef of the pointer to the routine	
	typedef void (A::* fType) ();
    
	//The pointer to which part of the routine is actually going to be called
	fType fPointer;
    
    public:
	//Constructor initializes the pointer, so that the first time the function is called
	//the first part is executed
        A()
	{
	    fPointer = &A::f_1;
        }
    
	//The actual function we are going to be calling each time (see main)
	inline void f()
	{
	    (this->*fPointer) ();
	}
	
	//Part 1 of the co-routine
        void f_1()
        {
	    cout << "f_1() is called." << endl;
    	
	    fPointer = &A::f_2;
	}
    
	//Part 2 of the co-routine
	void f_2()
	{
	    cout << "f_2() is called." << endl;
	
	    fPointer = &A::f_3;
	}
    
	//Part 3 of the co-routine
	void f_3()
	{
	    cout << "f_3() is called." << endl;
	}
};

int main(void)
{
    A myA;
    
    //f_1() is called
    myA.f();

    //f_2() is called
    myA.f();

    //f_3() is called
    myA.f();

    //f_3 was our last function in the series so it's called again,
    //this behavior obviously could be altered
    myA.f();
}
