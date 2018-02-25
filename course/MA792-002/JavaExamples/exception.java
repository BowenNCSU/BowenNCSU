public class exception {

static public class Integer {
// encapsulated to save a file, use Integer name
   private int x;
   public Integer() {x = 0;}
   public Integer(int i) {x = i;}
   public Integer div(Integer j) throws ArithmeticException
         {if(j.x == 0)
            throw new ArithmeticException(
                 "in Integer.div: division by zero");
          else return new Integer(this.x/j.x);
         }
   public String toString()
     {return String.valueOf(x); }
} // end Integer

public static void main(String[] args)
{try
 { Integer three = new Integer(3);
   Integer two = new Integer(2);
   Integer zero = new Integer();
   System.out.println(three.div(two).toString());
   System.out.println(three.div(zero).toString());
 }
 catch (ArithmeticException m)
     { System.out.print("Exception: ");
       System.out.println(m.getMessage());
     }
 finally
     { System.out.println("Continue...");
     }
}// end main
}// end exception
