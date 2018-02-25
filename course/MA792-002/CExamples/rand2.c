#include <stdio.h>
#include <stdlib.h>
/* more on John's question why the program cannot probe if
   the keyboard was struck:

   The GUI library is Unix's X-windows, and we run our program
   in an xterm.  The X-windows library provides an interupt
   handler for the keyboard, and we will use it in Java.

   However, the xterm program does not pass the keyboard interrupts
   to the currently running program.  Instead, it displays
   the struck keys on the line, and when a newline is struck
   it places them in a buffer that the program can read.  This
   is why you can edit your input line before hitting return.
   In no way does it signal new buffer contents to the running 
   program.

   xterm passes certain interrupts through to your program:
   if you strike ^C (control-C), an interrupt signal is sent.
                 ^Z (control-Z), a terminal stop (suspend)
                                 signal is sent.

   In the code below, I catch both signals before exiting.
*/

#include <signal.h>
void my_SIGINT_handler(int sig) /* interrupt by ^C */
  {/* psignal(sig, "in main()"); /* print message; not standard C */
   fprintf(stderr,"In my main() in rand2.c: %d\n\n",sig);
   exit(1);
  }

void my_SIGTSTP_handler(int sig) /* terminal stop by ^Z */
  {/* psignal(sig, "in main()"); /* print message; not standard C */
   fprintf(stderr,"In my main() in rand2.c: %d\n\n",sig);
   exit(1); /* should suspend here and wait for SIGCONT signal */
  }

int main(void){
 int i;

 /* variables that hold pointers to functions */
 void (*sys_SIGINT_handler_ptr)();
 void (*sys_SIGTSTP_handler_ptr)();

 /* activate signal handlers */
 sys_SIGINT_handler_ptr = signal(SIGINT, &my_SIGINT_handler);
 if( sys_SIGINT_handler_ptr == SIG_ERR)
   fprintf(stderr, "?Couldn't establish new SIGINT handler.\n");

 sys_SIGTSTP_handler_ptr = signal(SIGTSTP, &my_SIGTSTP_handler);
 if( sys_SIGTSTP_handler_ptr == SIG_ERR)
   fprintf(stderr, "?Couldn't establish new SIGTSTP handler.\n");

 srand(101);
 for (;;)
   {for (i=0; i<30; i++)
        /* print a random bit */
        printf("%1d", 1 & rand());
    printf("\n");
   }
 return 0;
}
