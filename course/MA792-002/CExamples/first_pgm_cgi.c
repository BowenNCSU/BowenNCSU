/* File: CExamples/first_pgm_cgi.c
   Compile as  "gcc -Wall -ansi -o first_pgm_cgi first_pgm_cgi.c"
   and run via "first_pgm.html" in your Internet browser (POST method)
*/

#include <stdio.h>
/* look in /usr/include for the contents of this header file */

int main(void)
{int i; /* loop counter */
 int c; /* number of times the text is to be printed */
 scanf("c=%d", &c);

 /* must give proper header information to browser */
 printf ("Content-type: text/plain\n\n");

 for(i=1; i<=c; i++)
    printf("%3d. I've been bad.\n",i);

 return 0;
}
