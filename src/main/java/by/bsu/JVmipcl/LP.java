package by.bsu.JVmipcl;

import java.io.IOException;

/**
 * @file LP.java  Interface for the `CLP` class
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
 
/**
 * This class has been designed for solving Linear Programs (LPs).
 *
 * \f{align*}{
 *     & c^Tx \to \max,\\
 *     & b_1 \le Ax \le b_2,\\
 *     & l \le  x \le u,
 * \f}
 * by both, prime and dual, simplex methods.
 */

public class LP
{
    /**
     * Infinity value.
     */
	public static final double INF = 1.0e20;

    /**
     * Infinity value for variables.
     */
	public static final double VAR_INF = 1.0e12;

    /**
     * Stores the pointer to a CLP object.
     */
	long nativeHd;

    /**
     *	The arrays `X`, `rdCost`, `shPrice`, `varHd` and `ctrHd` represent an optimal solution, where 
     *      - `X[i]` is the value of the variable with handle `varHd[i]`, `i=0,...,varHd.length`;
     *      - `rdCost[i]` is the reduced cost of the variable with handle `varHd[i]`, `i=0,...,varHd.length`;
     *      - `shPrice[i]` is the shadow price of the constraint (row) with handle `ctrHd[i]`, `i=0,...,ctrHd.length`.
     */
	double[] X, rdCost, shPrice;
	int[] varHd, ctrHd;

     /**
      * @return reference to the array storing the handles of variables.
      */
	public int[] getVarHandles()
	{
		return varHd;
	}

     /**
      * @return reference to the array storing the handles of constraints.
      */
	public int[] getCtrHandles()
	{
		return ctrHd;
	}

     /**
      * @return reference to the array storing the values of variables.
      */
	public double[] getX()
	{
		return X;
	}

     /**
      * @return reference to the array storing the reduced costs of variables.
      */
	public double[] getReducedCosts()
	{
		return rdCost;
	}

     /**
      * @return reference to the array storing the shadow prices (optimal values of dual variables).
      */
	public double[] getShadowPrices()
	{
		return shPrice;
	}

     /**
      * The defaul constructor.
      */
	public LP() {
        	nativeHd = 0;
		X = shPrice = rdCost = null;
		varHd = ctrHd = null;
   	}

    /**
     * The constructor that calls `createLP()` to initialize an "empty" CLP object.
     * @param name name of the LP to be solved.
     */
	public LP(String name) {
        	nativeHd = createLP(name);
		X = shPrice = rdCost = null;
		varHd = ctrHd = null;
   	}

    /**
     * The procedure creates an "empty" CLP object.
     * @param name name of the LP to be solved.
     */
	public native long createLP(String name);

    /**
     * Allocates memory for an LP of the required  size.
     * @param[in] m number of rows;
     * @param[in] n number of columns;
     * @param[in] nz number of non-zeroes in the matrix.
     * @throws CMemoryException lack of memory.
     */
	public native void openMatrix(int m, int n, int nz);
 
    /**
     * Does all preparations for solving the problem:
     * allocates memory, does preprocessing, scales the matrix.
     */
	public native void closeMatrix();

    /**
     * The function adds a variable (empty column) to the matrix.
     * @param[in] hd handle of variable;
     * @param[in] type type of variable;
     * @param[in] cost objective coefficient;
     * @param[in] l,u lower and upper bounds of variable;
     * @return index of added column (variable).
     * @throws CMemoryException lack of memory.
     */
 	public native int addVar(int hd, int type, double cost, double l, double u);

     /**
     * The function adds column to the matrix.
     * @param[in] hd handle of the variable;
     * @param[in] type type of variable;
     * @param[in] cost objective coefficient;
     * @param[in] l,u lower and upper bounds of variable;
     * @param[in] Val, Row arrays; `Val[i]` is coefficient in row `Row[i]`;
     * @param[in] bSort if `true`, column entries are sorted in increasing order of their row indices.
     * @return index of added column (variable).
     * @throws CMemoryException lack of memory.
     */
   	public native int addColumn(int hd, int type, double cost, double l, double u,
                                    double[] Val, int[] Row, boolean bSort);
 
    /**
     * The function adds a constraint (empty row) to the matrix.
     * @param[in] hd handle of the constraint;
     * @param[in] type type of the constraint;
     * @param[in] lhs,rhs  respectively, left (LHS) and right (RHS) hand sides of constraint.
     * @return constraint index.
     * @throws CMemoryException lack of memory.
     * @sa addRow().
     */
	public native int addCtr(int hd, int type, double lhs, double rhs);

    /**
     * This function is used to add a row (constraint) to the closed matrix.
     * @param[in] hd handle of row (constraint);
     * @param[in] type type of constraint;
     * @param[in] lhs,rhs respectively, left (LHS) and right (RHS) hand sides of constraint;
     * @param[in] Val, Col arrays, `Val[i]` is the coefficient at variable `Col[i]`;
     * @param[in] bSort if `true`, row entries are sorted in increasing order of their column indices.
     * @return constraint index.
     * @throws CMemoryException lack of memory.
     * @sa addCtr().
     */
	public native int addRow(int hd, int type, double lhs, double rhs,
                                 double[] Val, int[] Col, boolean bSort);

    /**
     * The functions adds to the matrix the entry of value `val` into row `i` and column `j`;
     * if `i < 0`, the function sets the cost of variable `j` to `val`.
     * @param[in] val entry value;
     * @param[in] i row index;
     * @param[in] j column index.
     * @throws CMemoryException lack of memory.
     */
	public native void addEntry(double val, int i, int j);

    /**
     * The function sets (or changes) the cost of a variable.
     * @param[in] j index of a variable.
     * @param[in] val cost of variable `j`.
     * \sa addEntry().
     */
	public native void setObjCoeff(int j, double val);

    /**
     * @param[in] sense  if `true`, the goal is to _maximize_ the objective function;
     * otherwise, to _minimize_.
     */
	public native void setObjSense(boolean sense);

    /**
     * @param[in] flag if set to `true`, all run time messages are blocked.
     */
	public native void beSilent(boolean flag);

    /**
     * Call this function to switch off all preprocessing actions, excluding scaling the matrix.
     */
	public native void preprocOff();

    /**
     * The procedure solves LPs.
     * @sa MIP.optimize().
     */
	public native void optimize();

    /**
     * @return `true` if the LP in the memory has been solved and an optimal solution found.
     * @sa MIP.isSolution().
     */
	public native boolean isSolution();
	
	/**
     * @return `true` if LP in memory has been solved and proven to be _infeasible_.
     */
    public native boolean isLpInfeasible();
	
	/**
     * @return `true` if LP in memory has been solved and proven to be _unbounded_.
     */
    public native boolean isLpUnbounded();
	
    /**
     * This function is usually called when the problem has been solved already.
     * @return optimal objective value of solved LP.
     */
	public native double getObjVal();

    /**
     * The function gets from the `CLP` object
     *     - variable handles (stores them in `varHd`);
     *     - constraint handles (stores them in `ctrHd`);
     *     - values of variables (stores them in `X`);
     *     - reduced costs (stores them in `rdCost`);
     *     - shadow prices (stores them in `shPrice`).
     * @remark The order in which the variables are listed in `varHd`
     *     is not necessary the same they were added to the matrix by `addVar()` or `addColumn()`.
     *   Similarly, the order in which the constraints are listed in `ctrHd`
     *     is not necessary the same they were added to the matrix by `addCtr()` or `addRow()`.
     * @sa MIP.getSolution().
     */
	public native void getSolution();

    /**
     * When the problem being solved has no solution,
     * this function gets from the `CLP` object a certificate of inconsistency
     *     - row multipliers (stores then in `shPrice`);
     *     - column multipliers (stores then in `rdCost`).
     */
	public native void getInconsistencyCertificate();

	/**
	 * The function writes the solution found by the __MIPCL__ solver  into a file.
	 * The user can overload this function to store solutions in an appropriate way.
	 * @param[in] filename name of the file to store solutions.
         * @throws CFileException.
	 */
	public native void printSolution(String filename) throws IOException;

    /**
     * The function destroys the CLP object which reference is stored in `nativeHd`.
    */
	public native void dispose();
}

