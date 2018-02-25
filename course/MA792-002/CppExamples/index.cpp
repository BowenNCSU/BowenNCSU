// Return-Path: <hitzm@cs.rpi.edu>
// Received: from vega.cs.rpi.edu by cs.rpi.edu (5.67a/1.4-RPI-CS-Dept)
// 	id AA04155; Wed, 28 Jun 1995 14:29:59 -0400 (hitzm from vega.cs.rpi.edu)
// Message-Id: <199506281829.AA04155@cs.rpi.edu>
// X-Mailer: exmh version 1.6.1 5/23/95
// To: kaltofen
// Cc: hitzm
// Subject: Re: STL program for index generation 
// In-Reply-To: Your message of "Wed, 28 Jun 1995 12:32:25 EDT."
//              <9506281632.AA13343@blackbox.cs.rpi.edu> 
// Mime-Version: 1.0
// Content-Type: text/plain; charset=us-ascii
// Date: Wed, 28 Jun 1995 14:29:57 -0400
// From: Markus Hitz <hitzm@cs.rpi.edu>
// Status: RO
// 
// 
// Erich:
// 
//   the program is attached below. You have to compile it with
// Apogee's "apCC". In order to get the STL files included, you
// can use the following one-line script:
// 
// 	apCC $1 $2 $3 $4 $5 $6 $7 $8 $9 -I/dept/cs/include/stl
// 
// (Write it to a file "astl" and chmod it +x). I think the program
// is too crude to be useful for other people right know. It could
// be easily improved by adding some "kill set" which would allow to
// filter out non-keywords like "the", as well as TeX and LaTeX commands.
// 
// Markus
// 
// --------------------------- C++ program --------------------------------- 
// 
//
// Generates "compressed" index of text files
//
// interface:	a.out	<filename>
//
// output:	file named "filename.index"
//

#include <ctype.h>
#include <iostream.h>
#include <fstream.h>
// using basic_string (EK 9/22/98)
#include <std/bastring.h>
// #include "mstring.h"
#include <set.h>

int main (int argc, char *argv[])
{
  char  buf[128];
  int   c, i=0;

  if (argc < 2)
  {
    cerr << argv[0]
         << ": missing argument" << endl;
    return 1;
  }

  ifstream  text (argv[1], ios::in);

  if (!text)
  {
    cerr << argv[0]
         << ": cannot open "
         << argv[1] << endl;
    return 2;
  }

  strcpy (buf, argv[1]);
  strcat (buf, ".index");

  ofstream  idxfile (buf, ios::out);

  if (!idxfile)
  {
    cerr << argv[0]
         << ": cannot open output file" << endl;
    return 3;
  }

  set <basic_string<char>,
       less<basic_string<char> >
      > idxset;

  while ((c = text.get ()) != EOF)
  {
    c = tolower (c);

    if (isalpha(c))
    {
      buf[i++] = c;
    }
    else if (i)
    {
      buf[i] = '\0';
      if (i > 2) idxset.insert (basic_string<char>(buf));
      i = 0;
    }
  }

  ostream_iterator<basic_string<char> > out (idxfile, "\n");
  copy (idxset.rbegin(), idxset.rend(), out);

  return 0;
}



