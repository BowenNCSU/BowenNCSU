// File: account.cpp
// -----------------------------------------------
//
// account class example
//
// -----------------------------------------------

#include <iostream>
#include <stdio.h>

// using namespace std;

// -----------------------------------------------

class account
{
    double balance;

  public:

    account (void) : balance(0) {}
    account (double start) : balance(start) {}

    void    deposit  (double);
    double  withdraw (double);

    double  getBalance (void);
};


// -----------------------------------------------
//
// implementation of member functions
// (definition)
//
// -----------------------------------------------


void account::deposit (double amount)
{
  this->balance += amount;
}


double account::getBalance (void)
{
  return balance;
}


double account::withdraw (double amount)
{
  if (balance >= amount)
  {
    balance -= amount;
    return amount;
  }
  else
  {
    std::cerr << "You are broke" << std::endl;
    amount = balance;
    balance = 0.0;
    return amount;
  }
}



// -----------------------------------------------

int main (void)
{
account  acc1;
account  acc2(5000.0);

  std::cout << "Balance (default constructor): ";
  std::cout << acc1.getBalance () << std::endl;

  acc1.deposit(1000);
  std::cout << "New balance acc1: " << acc1.account::getBalance () << std::endl;

  std::cout << "Balance (second constructor):  ";
  std::cout << acc2.getBalance () << std::endl;

  acc2.withdraw(10000);
  std::cout << "New balance acc2: " << acc2.getBalance () << std::endl;
}// end main

// -----------------------------------------------

