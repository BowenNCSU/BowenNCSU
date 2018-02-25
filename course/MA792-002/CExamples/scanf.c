#include <stdio.h>
#include <string.h>
#include <stdlib.h>

int main(int argc, char** argv)
/* scanf example program: run as "a.out < scanf.i > scanf.p
                              or "a.out -v" (-v for verbose) */
{char string[20];
 char format[200];
 int i, ret_code, count; double x;

 for(;;)
    {/* test if a prompt string is desired */
     if(argc > 1)
       {printf("Enter %%s %%d %%f\\n or ^D:");
        fflush(stdout); /* make sure prompt is written */
       }

     ret_code = scanf("%s%i%lf%n", string, &i, &x, &count);

     if(ret_code == EOF) break;
     if(ret_code < 3)
       {fprintf(stderr, "In main(): scanf unable to decipher! ");
        switch(ret_code){
        case 2: fprintf(stderr, "read integer:%d", i);
        case 1: fprintf(stderr, "read string:%s", string);
        default:fprintf(stderr, "\n");
                exit(1); /* indicate that an error occurred */
        } /* end switch */
       }
     printf(strcat(
               strcpy(format,
                      "Read string:\"%s\", integer:%d, float:%5.2f;"
                     ),
            "%2u char's total\n" ),
            string, i, x, count);
    }/* for(;;) */
 return 0;
}/*main*/
