/* File: CExamples/first_pgm_cgi_get.c
   Compile as  "gcc -Wall -ansi -o first_pgm_cgi_get first_pgm_cgi_get.c"
   and run via "first_pgm_get.html" in your Internet browser (GET method)
*/

#include <stdio.h>
#include <stdlib.h>
/* look in /usr/include for the contents of this header file */

int main(void)
{int i; /* loop counter */
 int c; /* number of times the text is to be printed */
 char *query_string; /* returned by getenv */

 query_string = getenv("QUERY_STRING");
 sscanf(query_string, "c=%d", &c);

 /* must give proper header information to browser */
 printf ("Content-type: text/plain\n\n");

 for(i=1; i<=c; i++)
    printf("%3d. I've been bad.\n",i);

 return 0;
}
