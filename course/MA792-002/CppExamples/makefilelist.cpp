#include <iostream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{basic_string<char> line;
 cout << "<html>" << endl;
 cout << "<h2>Directory listing</h2>" << endl;
 while(!cin.eof())
   {
    getline(cin, line);
    cout << "<a href=\"" << line << "\">"
         << line << "</a><br>" << endl;
   }
 cout << "</html>" << endl;
 return 0;
}
