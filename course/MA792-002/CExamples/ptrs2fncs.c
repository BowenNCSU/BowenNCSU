/* this example is taken from the Unix qsort man page */
#include <stdlib.h>
#include <stdio.h>

/* qsort declaration (in stdlib.h) */
extern void qsort(void*  A,  /* array */
                  size_t n, /* nr. of elements */
                  size_t s, /* size of each element */
                  int (*comp) /* pointer to comparison function */
                     (const void* ampA_i, /* pointer to array elem */
                      const void* ampA_j  /* pointer to array elem */
                     )
                 );

static  int intcompare(int *i, int *j)
{
        if (*i < *j)
                return (1);
        if (*i > *j)
                return (-1);
        return (0);
}
int main(void)
{
        int a[10] = {9,8,7,6,5,4,3,2,1,0};
        int i;

        qsort((void *) a, 10, sizeof(int),
              (int (*)(const void*, const void*)) intcompare
             );

        for (i=0; i<10; i++) printf(" %d",a[i]);
        printf("\n");
        return 0;
}
