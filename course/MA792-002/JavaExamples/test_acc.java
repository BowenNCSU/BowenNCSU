// -----------------------------------------------

import java.lang.*; // this is not necessary
//import account;

public class test_acc {
  public static void main (String[] args)
   {
   account  acc1 = new account();
   account  acc2 = new account(5000.0);
   
   System.out.print( "Balance (default constructor): ");
   System.out.println( acc1.getBalance() );
   
   acc1.deposit(1000);
   System.out.print( "New balance acc1: " );
   System.out.println( acc1.getBalance() );
   
   System.out.print( "Balance (second constructor):  ");
   System.out.println( acc2.getBalance() );
   
   acc2.withdraw(1000);
   System.out.print( "New balance acc2: ");
   System.out.println( acc2.getBalance () );
   }
}
// -----------------------------------------------
