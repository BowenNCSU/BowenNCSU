/* File CExamples/arrays.c */
#include <stdio.h>
#include <stdlib.h> /* for calloc */
int main(void)
{int digits[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
 /* the dimension is optional and will be determined by the compiler
    note that char name[] = "Erich";      is a shorthand for
              char name[6] = {'E', 'r', 'i', 'c', 'h', '\0'}; */
 float* xptr[5]; /* array of pointers to float */
 float xy[5][6] = {{0.0, 0.1, 0.2, 0.3, 0.4, 0.5},
                   {1.0, 1.1, 1.2, 1.3, 1.4, 1.5},
                   {2.0, 2.1, 2.2, 2.3, 2.4, 2.5},
                   {3.0, 3.1, 3.2, 3.3, 3.4, 3.5},
                   {4.0, 4.1, 4.2, 4.3, 4.4, 4.5}};
 int i, j;

 printf("digits[4] == %d, *(digits+4) == %d\n\n", digits[4], *(digits+4) );

 printf("&xy[0][0] == %p,     &xy[2][0] == %p\n", &xy[0][0], &xy[2][0] );
 printf("(float*)xy == %p,    (float*) (xy+2) == %p\n",
        (float*)xy, (float*)(xy+2) );
 printf("(float*)xy[2] == %p, (float*)*(xy+2) == %p\n\n",
        (float*)xy[2], (float*)*(xy+2) );

 printf("1. xy[2][3]                == %4.1f\n\n", xy[2][3]);
 printf("2. *(xy[2]+3)              == %4.1f\n\n", *(xy[2]+3));
 printf("3. *(*(xy + 2) + 3)        == %4.1f\n\n", *(*(xy + 2) + 3));
 printf("4. *((float*)xy + 2*6 + 3) == %4.1f\n\n", *((float*)xy + 2*6 + 3));
 printf("5. (*(xy+2))[3]            == %4.1f\n\n",
        (*(xy+2))[3]); /* Wen-shin's way */
 printf("6. *(&xy[0][0] + 2*6 + 3)  == %4.1f\n\n",
        *(&xy[0][0] + 2*6 + 3)); /* the 6th version */ 
 printf("7. *(float*)((unsigned char*)xy + sizeof(float)*(2*6 + 3)) == %4.1f\n\n",
        *(float*)((unsigned char*)xy + sizeof(float)*(2*6 + 3))); /* the 7th ver. */
 printf("8. *(float*)((unsigned char*)xy[2] + sizeof(float)*3)      == %4.1f\n\n\n",
        *(float*)((unsigned char*)xy[2] + sizeof(float)*3) );

 /* initialize the array of pointers */
 for(i=0; i<5; i++)
    {xptr[i] = (float *)calloc(6, sizeof(float));
     for(j=0; j<6; j++)
        *(xptr[i] + j) = i + 0.1 * j;
    } /* end for i */

 printf("1. xptr[2][3]         == %4.1f\n\n", xptr[2][3] );
 printf("2. *(xptr[2]+3)       == %4.1f\n\n", *(xptr[2]+3) );
 printf("3. *(*(xptr + 2) + 3) == %4.1f\n\n", *(*(xptr + 2) + 3) );
 printf("5. (*(xptr+2))[3]     == %4.1f\n\n", (*(xptr+2))[3] ); /* Wen-shin's way */
 printf("8. *(float*)((unsigned char*)xptr[2] + sizeof(float)*3) == %4.1f\n\n",
        *(float*)((unsigned char*)xptr[2] + sizeof(float)*3) );

 for(i=0; i<5; i++) free(xptr[i]);
 return 0;
}
