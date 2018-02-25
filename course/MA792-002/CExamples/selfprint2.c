#include <stdio.h>
#include <string.h>

char* S[] = {"};",
"/* the array of strings S is a",
" * representation of the body",
" * of the program from '0'",
" * to the end",
" */",
"\n",
"int main(void)",
"{int i, j;",
" printf(\"#include <stdio.h>\\n#include <string.h>\\n\");",
" printf(\"\\nchar* S[] = {\");",
" for(i=0; S[i]; i++)",
"   {printf(\"\\\"\");",
"    for(j=0; j <= strnlen(S[i],100)-1; j++)",
"       /* treat \", \\, \\n special for the",
"        * declaration of S",
"        */",
"       if(S[i][j] == '\\n') printf(\"\\\\n\");",
"       else if(S[i][j] == '\"') printf(\"\\\\\\\"\");",
"       else if(S[i][j] == '\\\\') printf(\"\\\\\\\\\");",
"       else printf(\"%c\",S[i][j]);",
"    printf(\"\\\",\\n\");",
"   }/*for(i)*/",
" printf(\"0\\n\");",
"\n",
" for(i=0; S[i]; i++)",
"    if (S[i][0] == '\\n') printf(\"\\n\");",
"    else printf(\"%s\\n\",S[i]);",
"}",
0
};
/* the array of strings S is a
 * representation of the body
 * of the program from '0'
 * to the end
 */

int main(void)
{int i, j;
 printf("#include <stdio.h>\n#include <string.h>\n");
 printf("\nchar* S[] = {");
 for(i=0; S[i]; i++)
   {printf("\"");
    for(j=0; j <= strnlen(S[i],100)-1; j++)
       /* treat ", \, \n special for the
        * declaration of S
        */
       if(S[i][j] == '\n') printf("\\n");
       else if(S[i][j] == '"') printf("\\\"");
       else if(S[i][j] == '\\') printf("\\\\");
       else printf("%c",S[i][j]);
    printf("\",\n");
   }/*for(i)*/
 printf("0\n");

 for(i=0; S[i]; i++)
    if (S[i][0] == '\n') printf("\n");
    else printf("%s\n",S[i]);
}
