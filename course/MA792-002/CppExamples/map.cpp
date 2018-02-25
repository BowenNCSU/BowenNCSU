#include <map>
#include <iostream>

using namespace std;

int main(void)
{typedef map<char, int, less<char> > mymap;
 typedef pair<const char, int> mypair;
 typedef pair<mymap::iterator, bool> mymap_insert_t;

 mymap M;
 M.insert(mypair('A', 25));
 M.insert(mypair('B', 10));
 mypair x = mypair('A', 26);
 mymap_insert_t p = M.insert(x);
 if (p.second == false) cout << x.first << "," << x.second << " not inserted" << endl;
 cout << M['A'] << endl;
 cout << M['B'] << endl;
 cout << M['C'] << endl;

 typedef multimap<char, int, less<char> > mymultimap;
 typedef mymultimap::iterator mymultimap_insert_t;
 typedef pair<mymultimap::const_iterator,mymultimap::const_iterator> mymultimap_equal_range_t;

 mymultimap MM;
 MM.insert(mypair('A', 25));
 MM.insert(mypair('B', 10));
 MM.insert(mypair('A', 26));

 mymultimap_equal_range_t r = MM.equal_range('A');
 mymultimap::const_iterator pp = r.first;
 cout << (*pp).first << ", " << (*pp).second << endl;
 ++pp;
 cout << (*pp).first << ", " << (*pp).second << endl;
}
