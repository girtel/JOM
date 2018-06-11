/**
 * @file Ctr.java  Interface for the `Ctr` class.
 * <table border="2">
 * <tr>
 *   <td>  <b>Author</b> </td>  <td> Nicolai N. Pisaruk </td>
 * </tr> 
 * <tr>
 *  <td> <b> e-mail </b> </td>  <td> nicolaipisaruk@gmail.com </td>
 * </tr>
 * <tr>
 *  <td> <b>home page</b> </td> <td> wwww.mipcl-cpp.appspot.com </td>                |
 * </tr>
 * <tr>
 *    <td colspan="2"> <b>&copy; 2017</b> Nicolai N. Pisaruk </td>
 * </tr>
 * </table>
 */

package by.bsu.JVmipshell;

import java.util.List;
import java.util.ArrayList;

import by.bsu.JVmipcl.*;
// import by.bsu.JVmipshell.Var;
import by.bsu.JVmipshell.LinSum;

/**
 *     A `Ctr` object is a linear sum (`LinSum` object) with lower and upper bounds on its values.
 */
public class Ctr
{
    /**
     * Linear sum (linear function).
     */
	LinSum lsum;

    /**
     * Left hand side.
     */
	double lhs;

    /**
     * Right hand side.
     */
        double rhs;

    /**
     * The handle that uniquely identifies the constraint.
     */
        int hd;

    /**
     * Shadow price value.
     */
	double price;

    /**
     * @return reference to the `Linsum` object of this constraint. 
     */
	public LinSum getlsum()
	{
		return lsum;
	}
	
    /**
     * The function sets the reference to a `LinSum` object, which is the middle part of the constraint.
     * @param lsum `LinSum` object.
     */
	public void setlsum(LinSum lsum)
	{
		this.lsum = lsum;
	}

    /**
     * @return left hand side of this constraint. 
     */
	public double getlhs()
	{
		return lhs;
	}

    /**
     * The function sets the left hand side for this constraint.
     * @param lhs left hand side value. 
     */
	public void setlhs(double lhs)
	{
		this.lhs = lhs;
	}

    /**
     * @return right hand side of this constraint. 
     */
	public double getrhs()
	{
		return rhs;
	}

    /**
     * The function sets the right hand side for this constraint.
     * @param rhs right hand side value. 
     */
	public void setrhs(double rhs)
	{
		this.rhs = rhs;
	}

    /**
     * @return handle of this constraint. 
     */
	public int gethd()
	{
		return hd;
	}

    /**
     * The function sets the handle (unique identifier) for this constraint.
     * @param hd handle value. 
     */
	public void sethd(int hd)
	{
		this.hd = hd;
	}

    /**
     * @return shsdow price for this constraint. 
     */
	public double getprice()
	{
		return price;
	}

    /**
     * The function sets the value of shadow price for this constraint.
     * @param price shadow price. 
     */
	public void setprice(double price)
	{
		this.price = price;
	}

    /**
     * The constructor initializes a newly created `Ctr` object.
     *
     * Usually, this constructor is called from `LPshell::addCtr()`.
     *  
     * @param[in] hd handle;
     * @param[in] lsum represents a linear function which values be restricted; 
     * @param[in] lhs left hand side;
     * @param[in] rhs right hand side.
     * @param[in] lhs left hand side;
     */
	public Ctr(int hd, LinSum lsum, double lhs, double rhs)
	{
	        this.hd = hd;
	        this.lsum = lsum;
	        this.lhs = lhs;
	        this.rhs = rhs;
	}

    /**
     * Compute a string representation of the `this` object.
     *
     * @param[in] left if `true`, the left hand side part of `this` constraint is represented by the return string;
     * @param[in] right if `true`, the right hand side part of 'this' constraint is represented by the return string;
     * @return string representing `this` constraint.
     */ 
	public String toString(boolean left, boolean right)
        {
		StringBuilder s = new StringBuilder(64);
        	if (rhs-lhs < 1.0e-10) {
            		s.append(lsum);
			s.append(" == ");
			s.append(rhs);
		}
        	else {
            		if (left && (lhs > -LP.INF+1.0)) {
				s.append(lhs);
				s.append(" <= ");
			}
            		s.append(lsum);
            		if (right && (rhs < LP.INF-1.0)) {
				s.append(" <= ");
				s.append(rhs);
			}
		}
        	return s.toString();
	}

	
    /**
     * The function computes string representations of `Ctr` objects.
     *
     * @return string representing `this` constraint.
     */ 
	public String toString()
        {
		return toString(true,true); 
	}
}    

