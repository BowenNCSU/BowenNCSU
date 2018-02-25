#include <math.h> /* for powl() etc. */
#include <stdio.h>
#include <stdint.h> /* for typedef uint64_t */
#include <limits.h>
/* compile this program as "gcc -Wall -ansi -std=c99 math.c -lm" */
int main(void)
{uint64_t two_to_the_e;
 long double two_to_the_64_minus_one;
 size_t nr_of_bytes;
 int e, exp=64;

 nr_of_bytes = sizeof(uint64_t);

 /* convert size_t to int for %d */
 printf("sizeof(uint64_t) = %d\n", (int)nr_of_bytes);

 for(e=31; e <= exp; e++)
    {two_to_the_e = (uint64_t)(powl(2,e)+0.5);
     printf("2 to the %d == %Lf, %lu\n",
             e, powl(2,e), two_to_the_e);
    }/*for(e=31...)*/

 two_to_the_64_minus_one = powl(2,64)-1.0L;
 printf("(2 to the %d)-1 == %Lf, %lu\n",
        64, two_to_the_64_minus_one, (uint64_t)two_to_the_64_minus_one);

 printf("ULLONG_MAX = %llu\n", ULLONG_MAX);
 return 0;
}
