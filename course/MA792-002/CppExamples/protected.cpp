#include <iostream>


struct Toto {
        Toto(int b = 0) : a(b) {}
protected:
        int a;
};

struct Spy {

    class InHeritsToto : Toto {
	    int x;
    protected:
       friend struct Spy;
    };
    int getprotected( const Toto & r ) {
//          return r.a ;  // Error : `int Toto::a' is protected
            return static_cast<const InHeritsToto&>(r).x ;
            return static_cast<const InHeritsToto&>(r).a ;
    }
    int& assignprotected( Toto & r, int b ) {
//          return r.a = b; // Error : `int Toto::a' is protected
            return static_cast<InHeritsToto&>(r).a = b;
    }
};

int main(int argc, char * argv []) {

        Toto t(3);

        std::cout << Spy().getprotected( t )  << std::endl;

        Spy().assignprotected( t, 7 ) ;

        std::cout << Spy().getprotected( t )  << std::endl;

        return 0;
}
