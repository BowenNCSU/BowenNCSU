/* File CExamples/mainarg.c */
#include <stdio.h>
#include <string.h> /* for strlen */
#include <ctype.h>  /* for isprint */

int main(int argc, char *argv[])
/* This program manipulates the arguments to main for purpose of
   demonstrating the use of characters, strings, and arrays of strings,
   pointers and their integer values in C */
{int i;
 unsigned char *cptr;
 char *lastarg_firstcptr; /* points to first character of last argument */
 char *lastarg_termincptr; /* points to null char terminating last argument */

 printf("argc == %2d\n", argc);
 /* print the location of the argument string in hexadecimal,
    and then print the string */
 for (i=0; i < argc; i++)
     printf("Address of argv[%1d] == %p, string value == \"%s\"\n",
             i, argv[i], argv[i] );
 printf("\n");

 /* Now print each character in the range of the first character of
    the first argument to the null character terminating the last argument */
 printf("argv == %p, argv[0] == %p\n", argv, argv[0]);
 lastarg_firstcptr = argv[argc-1];
 lastarg_termincptr = lastarg_firstcptr + strlen(lastarg_firstcptr);
 printf("lastarg_termincptr = %p\n", lastarg_termincptr);

 if (*argv <= lastarg_termincptr)
    for(cptr = (unsigned char*)*argv; cptr <= (unsigned char*)lastarg_termincptr; cptr++)
       if(isprint((int)*cptr))
          printf("Address == %p, character value == '%c'\n",
                 cptr, *cptr);
       else
          printf("Address == %p, character value == '\\%03o'\n",
                 cptr, *cptr);
 else
    for(cptr = (unsigned char*)argv[argc-1];
        cptr <= (unsigned char*)(argv[0]+strlen(argv[0]));
        cptr++
       )
       if(isprint((int)*cptr))
          printf("Address == %p, character value == '%c'\n",
                 cptr, *cptr);
       else  
          printf("Address == %p, character value == '\\%03o'\n",
                 cptr, *cptr);
 return 0;
}
