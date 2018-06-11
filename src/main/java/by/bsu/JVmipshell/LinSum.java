/**
 * @file LinSum.java  Interface for the `LinSum` class
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
import by.bsu.JVmipshell.Var;

/**
 * A term is a pair of a real coefficient and a variable.
 */
class Term
{
	double coeff;
	Var var;

	public Term(double coeff, Var var)
	{
		this.coeff = coeff;
		this.var = var;
	}

	public double getcoeff()
	{
		return coeff;
	}

	public void setcoeff(double coeff)
	{
		this.coeff = coeff;
	}

	public Var getvar()
	{
		return var;
	}

	public void setvar(Var var)
	{
		this.var = var;
	}

	public String toString()
	{
		StringBuilder s = new StringBuilder(64);
		s.append(coeff);
		s.append('*');
		s.append(var.getname());
		return s.toString();		
	}
}

/**
 * __LinSum__ objects represent linear function.
 *
 *   Each _LinSum_ object is just a list of terms.
 */
public class LinSum
{
    /**
     *  List of Term objects that is a linear sum.
     */
	List<Term> terms;

    /**
     * @return list of terms in `this` linear sum.
     */
	public List<Term> getterms()
	{
		return terms;
	}

    /**
     * Initialize an empty LinSum object.
     */
	public LinSum()
	{
		terms = null;
	}

    /**
     * The functions makes this object ready for building another linear sum.
     *
     * This function is usually called any time an object referencing to
     * `this` `LinSum` object is added to a `LPshel` or `MIPshell` problem.  
     */
	public void reset()
	{
		terms = null;
	}

    /**
     * The function computes a string representation of the `this` object.
     *
     * @return string representation of the `this` object.
     */
	public String toString()
	{
		if (terms != null) {
        		StringBuilder s = new StringBuilder(64);
        		int k = 0;
        		for (Term t : terms) {
				if ((k > 0) && (t.getcoeff() >= 0.0))
                			s.append('+');
            			s.append(t.toString());
            			++k;
			}
        		return s.toString();
		}
		else
			return "";
	}

    /**
     * The function changes the signs of coefficients in all terms of `this` linear sum.
     *
     * @return reference to `this` object.
     */
	public LinSum negate()
	{
        	for (Term t : terms) {
            		t.setcoeff(-t.getcoeff());
		}
		return this;
	}

    /**
     * The function multiplies the coefficients of this `LinSum` object by a give number.
     *
     * @param[in] a multiplier.
     * @return reference to `this` object.
     */ 
	public LinSum multiply(double a)
	{
	        for ( Term t : terms)
            		t.coeff *= a;
		return this;
	}

    /**
     * The function computes the sum of `this` and a given `LinSum` objects.
     *
     * @param[in] lsum reference to `LinSum` object.
     * @return reference to `this` object.
     */
	public LinSum add(LinSum lsum)
	{
		if (terms == null)
			terms = new ArrayList<>();
		terms.addAll(lsum.terms);
		lsum.terms = null;
		return this;
	}

    /**
     * The function adds a term to `this` `LinSum` objects.
     *
     * @param[in] a coefficient;
     * @param[in] var reference to variable or `null`.
     * @return reference to `this` object.
     */
	public LinSum add(double a, Var var)
	{
		if (terms == null)
			terms = new ArrayList<>();
		terms.add(new Term(a,var));
		return this;
	}

    /**
     * The function adds a term to `this` `LinSum` objects.
     *
     * @param[in] term reference to `Term` object.
     * @return reference to `this` object.
     */
	public LinSum add(Term term)
	{
		if (terms == null)
			terms = new ArrayList<>();
		terms.add(term);
		term = null;
		return this;
	}

    /**
     * The function adds a given number to computes `this` LinSum object.
     *
     * @param[in] a real number.
     * @return reference to `this` object.
     */
	public LinSum add(double a)
	{
		if (terms == null)
			terms = new ArrayList<>();
		terms.add(new Term(a,null));
		return this;
	}

    /**
     * The function substracts a given `LinSum` object from `this` object.
     *
     * @param[in] lsum reference to LinSum object.
     * @return reference to `this` object.
     */
	public LinSum sub(LinSum lsum)
	{
		if (terms == null)
			terms = new ArrayList<>();
		for (Term t : lsum.terms) {
	                t.coeff = -t.coeff;
	                terms.add(t);
		}
		lsum.terms = null;
		return this;
	}
}

