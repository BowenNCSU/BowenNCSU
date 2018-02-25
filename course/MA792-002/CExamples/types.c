#include <limits.h>
#include <stdio.h>

/* The following preprocessor macro defines a user set type integer
   that is guaranteed to have sizeof(integer) == 4, i.e. 32 bits */
#if SHORT_MAX == 2147483647
    #define integer short
#elif INT_MAX == 2147483647
    #define integer int
#elif LONG_MAX == 2147483647
    #define integer long
#else
    #error "No integer type with 4 bytes"
    #define integer void
#endif

/* Now define the range of type integer */
#define INTEGER_MAX 2147483647  /* exp2(31) - 1; see math.h */
#define INTEGER_MIN -2147483648 /* exp2(31)     */
#define UINTEGER_MAX 4294967295u /* exp2(32) - 1 */

/* This program exhibits some subtle points for the basic C types */

int main(void)
{short s = SHRT_MAX; int i = INT_MIN; unsigned u = INT_MAX;
 long L = LONG_MAX; unsigned long uL = ULONG_MAX;
 unsigned integer uinteg = UINTEGER_MAX;
 char c = CHAR_MAX;

 /* short, int, long, unsigned, unsigned long */
 printf("sizeof(short) == %lu, SHRT_MAX == %d, ++s == %d\n\n",
        sizeof(short), SHRT_MAX, ++s);
 printf("sizeof(int) == %lu, INT_MIN == %d, INT_MAX == %d, --i == %d\n\n",
        sizeof(int), INT_MIN, INT_MAX, --i);
 printf("sizeof(unsigned) == %lu, UINT_MAX == %u, ++u == %u\n\n",
        sizeof(unsigned), UINT_MAX, ++u);
 printf("sizeof(long) == %lu, LONG_MAX == %li, ++L == %li\n\n",
        sizeof(unsigned), LONG_MAX, ++L);
 printf("sizeof(unsigned long) == %lu, ULONG_MAX == %lu, ++uL == %lu\n\n",
        sizeof(unsigned long), ULONG_MAX, ++uL);

 /* Test macro for type integer */
 printf("sizeof(integer) == %lu, uinteg == %u\n\n", sizeof(integer), uinteg);

 /* char type */
 printf("sizeof(char) == %lu, CHAR_MAX == %d\n", sizeof(char), CHAR_MAX);
 printf("(int) c == %d, (int)((signed char)(c << 1)) == %d, (unsigned)((unsigned char)(c << 1)) == %u\n",
        c, (signed char)(c << 1), (unsigned char)(c << 1));
 printf("(int) (c >> 1) == %d, (int) ++c == %d, (unsigned)((unsigned char)c) == %u\n",
        c >> 1, ++c, (unsigned char)c);
 return 0;
}
