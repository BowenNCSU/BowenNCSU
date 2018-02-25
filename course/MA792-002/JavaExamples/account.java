// File: account.java
// -----------------------------------------------
//
// account class example
//
// -----------------------------------------------

import java.lang.System; // for System.err  (imported by default)

public class account
{
    double balance;

    public account () { balance = 0.0; }
    public account (double start) { balance = start; }

    void    deposit  (double amount) { balance += amount; }
    double  withdraw (double amount)
      { if (balance >= amount)
        {
          balance -= amount;
          return amount;
        }
        else
        {
          System.err.println("You are broke!");
          amount = balance;
          balance = 0.0;
          return amount;
        }
      }

    double  getBalance () { return balance; }
}
