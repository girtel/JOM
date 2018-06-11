/**
 * @file Var.java  Interface for the `Var` class.
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

import by.bsu.JVmipshell.LPshell;
import by.bsu.JVmipshell.MIPshell;


/**
 *   In MIPCL-JAVA models, objects of type `Var` represent variables.
 */
public class Var
{

// Attributes:
	/**
	 * Names of variables are used when printing a solution.
	 */
	String name;

	/**
         * Lower and upper bounds.
	 */
	double lb , ub;

	/**
	 * Value, which is computed only if a solution has been found.
	 */
	double val;

	/**
	 *  Handle of this variable, which must unique.
	 */
        int hd;

    /**
     * Type of this variable.
     */
        int type;

    /**
     * Priority of this variable (in range from -50 to 50).
     */
	int priority;

    /**
     * The constructor that creates a new variable, and adds that variable to a specified problem.
     *
     * @param[in] prob LP or MIP object this var to be added;
     * @param[in] name name, which is used when printing a solution;
     * @param[in] type type (`VAR_REAL`, `VAR_INT`, or `VAR_BIN`) of variable being created;
     * @param[in] lb lower bound;
     * @param[in] ub upper bound.
     */
	public Var(LPshell prob, String name, int type, double lb, double ub)
	{
        	this.name = name;
        	this.type = type;
        	if ((type & MIPshell.VAR_BIN) != 0) {
            		this.lb = 0.0;
			this.ub = 1.0;
		}
        	else {
            		this.lb = lb;
			this.ub = ub;
		}
        	priority = 0;
		prob.addVar(this);
	}

    /**
     * This constructor is used to create ether binary variables,
     * or variables (real or integer) that takes values from `0` to infinity.
     *
     * @param[in] prob LP or MIP object this var to be added;
     * @param[in] name name, which is used when printing a solution;
     * @param[in] type type (`VAR_REAL`, `VAR_INT`, or `VAR_BIN`) of variable being created.
     */
	public Var(LPshell prob, String name, int type)
	{
        	this.name = name;
        	this.type = type;
		this.lb =0.0;
        	ub = ((type & MIPshell.VAR_BIN) != 0)? 1.0: LPshell.VAR_INF;
        	priority = 0;
		prob.addVar(this);
	}


    /**
     * @return name of `this` variable.
     */
	public String getname()
	{
		return name;
	}

    /**
     * Set the name for `this` variable.
     *
     * @param[in] name new name.
     */
	public void setname(String name)
	{
		this.name = name;
	}

    /**
     * @return handle of `this` variable.
     */
	public int gethd()
	{
		return hd;
	}

    /**
     * Set the handle for `this` variable.
     *
     * @param[in] hd handle.
     */
	public void sethd(int hd)
	{
		this.hd = hd;
	}

    /**
     * @return type of `this` variable.
     */
	public int gettype()
	{
		return type;
	}

    /**
     * Set the type for `this` variable.
     *
     * @param[in] type new type.
     */
	public void settype(int type)
	{
		this.type = type;
	}

    /**
     * @return priority of `this` variable.
     */
	public int getpriority()
	{
		return priority;
	}

    /**
     * Set the priority for `this` variable.
     *
     * @param[in] priority new priority (integer in the gange from `-50` to `50`).
     */
	public void setpriority(int priority)
	{
		this.priority = priority;
	}

    /**
     * @return lower bound of `this` variable.
     */
	public double getlb()
	{
		return lb;
	}

    /**
     * Set the lower for `this` variable.
     *
     * @param[in] lb lower bound.
     */
	public void setlb(double lb)
	{
		this.lb = lb;
	}

    /**
     * @return upper bound of `this` variable.
     */
	public double getub()
	{
		return ub;
	}

    /**
     * Set the upper bound for `this` variable.
     *
     * @param[in] ub upper bound.
     */
	public void setub(double ub)
	{
		this.ub = ub;
	}

    /**
     * Set the lower and upper bounds for `this` variable.
     *
     * @param[in] lb lower bound;
     * @param[in] ub upper bound.
     */
	public void setbounds(double lb, double ub)
	{
		this.lb = lb;
		this.ub = ub;
	}

    /**
     * @return value of `this` variable.
     */
	public double getval()
	{
		return val;
	}


    /**
     * Set the value of `this` variable.
     *
     * @param[in] val value of `this` variable.
     */
	public void setval(double val)
	{
		this.val = val;
	}

    /**
     * The function increments the value of `this` variable.
     *
     * $param[in] q number to be added to `this` variable.
     */
	public void incval(double q)
	{
		val += q;
	}

    /**
     * The function computes a string representation of the `this` object.
     *
     * @return string representation of the `this` object.
     */
	public String toString()
	{
		StringBuilder s = new StringBuilder(64);
	        if ((type & MIPshell.VAR_BIN) != 0) {
			s.append(name);
			if (ub - lb > 0.5)
                		s.append(" in {0,1}");
            		else {
                		s .append(" = ");
				s.append(ub);
			}
		}
        	else {
            		if (lb > -LPshell.VAR_INF+1.0)
                		s.append(lb);
            		else
                		s.append("-inf");
			s.append(" <= ");
			s.append(name);
			s.append(" <= ");			
            		if (ub < LPshell.VAR_INF-1.0)
                		s.append(ub);
            		else
                		s.append("inf");
            		if ((type & MIPshell.VAR_INT) != 0)
                		s.append(" : integer");
		}
        	return s.toString();
	}
}

