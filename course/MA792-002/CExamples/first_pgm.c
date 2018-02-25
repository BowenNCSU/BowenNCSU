/* file: CExamples/first_pgm.c
   Programmer: Erich Kaltofen
   Version: 1.0
   Compile as  "gcc -Wall -ansi first_pgm.c"
   and run by  "a.out"           */

#include <stdio.h>
/* look in /usr/include for the contents of this header file */

int main(void)
{int i; /* loop counter */
 int c; /* number of times the text is to be printed */
 printf("Enter a number: ");
 fflush(stdout); /* flushes the output buffer */
 scanf("%d", &i);
 for(c=1; c<=i; c++)
    printf("%3d. I've been bad.\n",c);
 return 0;
}
