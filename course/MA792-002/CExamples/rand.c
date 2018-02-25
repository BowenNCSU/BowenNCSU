#include <stdio.h>
#include <stdlib.h>


int main(int argc, char* argv[]){
 int i;
 int seed;
 if(argc < 2)
   {fprintf(stderr,"Error, call: %s <seed>\n",
            argv[0]);
    return 1;
   }/*end if*/
 sscanf(argv[1],"%d",&seed);
 srand(seed);
 for (;;)
   {for (i=0; i<49; i++)
        /* print a random bit */
        printf("%1d", 1 & rand());
    printf("\n");
   }
 return 0;
}
