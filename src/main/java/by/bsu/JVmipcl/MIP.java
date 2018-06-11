package by.bsu.JVmipcl;

import java.io.IOException;

/**
 * @file MIP.java  Interface for the `CMIP` class
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
 * This class has been designed for solving Mixed Integer Programs (MIPs).
 *
 * \f{align*}{
 *     & c^Tx \to \max,\\
 *     & b_1 \le Ax \le b_2,\\
 *     & l \le  x \le u,\\
 *     & x_i\in\mathbb{Z}\quad\forall\; i\in I,
 * \f}
 * by the branch-and-cut method.
 *
 * __Examples__ of how to use the `MIP` class:
 *    - \ref JVprimer;
 *    - \ref JVfcnf;
 *    - \ref JVinfLP.
 */
public class MIP extends LP
{
    /**
     * In __MIPCL__ this flag indicates that a variables takes only integer values. 
     */
	public static final int VAR_INT = 0x00001000;

    /**
     * In __MIPCL__ this flag indicates that a variables is _binary_, i.e., it takes only two values, `0` or `1`. 
     */
	public static final int VAR_BIN = 0x00002000;

    /**
     * The constructor that calls `createMIP()` to initialize an "empty" CMIP object.
     * @param name name of the MIP to be solved.
     */
	public MIP(String name) {
		super();
        	nativeHd = createCMIP(name);
   	}

    /**
     * The procedure creates an "empty" CMIP object.
     * @param name name of the MIP to be solved.
     */
	public native long createCMIP(String name);

    /**
     * Allocates memory for a MIP of the required  size.
     * @param[in] m number of rows;
     * @param[in] n number of columns;
     * @param[in] nz number of non-zeroes in the matrix.
     * @throws CMemoryException lack of memory.
     */
	@Override
        public native void openMatrix(int m, int n, int nz);

    /**
     * Does all preparations for solving the problem:
     * allocates memory, does preprocessing, scales the matrix.
     */
	@Override
        public native void closeMatrix();

    /** 
     * Integer variables can be assigned priorities, i.e.,
     * integer numbers from `CMIP::VAR_PRI_MIN` to `CMIP::VAR_PRI_MAX`,
     * which are used when selecting a variable for branching.
     * The higher priority the more chances for a fractional variable to be chosen.
     * @param[in] j index of variable which priority to be changed
     * @param[in] pr new priority
     * @sa incVarPriority(), getVarPriority().
     */
	public native void setVarPriority(int j, int pr);

    /**
     * The solver automatically generates cuts for a given number of initially processed nodes,
     * and also for nodes of height less or equal than another given number.
     * @param[in] nodeNum number of nodes for which cuts are generated;
     * @param[in] height cuts are generated for nodes of height `<= height`.
     * @sa setCutTypePattern(), setAutoCutRounds().
     */
	public native void setAutoCutPattern(int nodeNum, int height);

    /**
     * The procedure solves MIPs.
     * @param[in] timeLimit limit on solution time (in seconds);
     * @param[in] gap integrality gap;
     * @param[in] solFileName pointer to string with file name for storing intermediate solutions.
     * @remark When solving LPs, input parameters are not used.
     * @sa MIP.optimize().
     */
	public native void optimize(long timeLimit, double gap, String solFileName);

    /**
     * The procedure solves LPs, and it is equivalent to the `optimize()`
     * with three arguments when the latter are set to the followiing default values:
     *      - `timeLimit = 10000000`;
     *      - `gap = 0.0`;
     *      - `solFileName = null`.
     * @sa LP.optimize().
     */
	@Override
	public native void optimize();

    /**
     * @return `true` if the MIP in the memory has been solved and an optimal solution found.
     * @sa LP.isSolution().
     */
	@Override
	public native boolean isSolution();

    /**
     * First, `isSolution()` is called; if it returns `true`,
     * we call `isSolutionOptimal()` to verify whether
     * the solution found by __MIPCL__ is optimal.
     * @return `true`, if an optimal solution has been found.
     */
	public native boolean isSolutionOptimal();
	
    /**
     * @return `true`, if this MIP is infeasible (has no solution).
     */
	public native boolean isInfeasible();
	
	/**
     * @return `true`, if solution procedure stopped after exceeding given time limit.
     */
	public native boolean timeLimitStop();

    /**
     * This function is usually called when the problem has been solved already.
     * @return optimal objective value of solved MIP.
     */
	@Override
	public native double getObjVal();

    /**
     * The function gets from the `CMIP` object
     * the variable handles (stores them in `varHd`),
     * and the values of variables (stores them in `X`).
     * @remark The order in which the variables are listed in `varHd`
     *     is not necessary the same they were added to the matrix by `addVar()` or `addColumn()`.
     * @sa LP.getSolution().
     */
	@Override
	public native void getSolution();

	/**
	 * The function writes the solution found by the __MIPCL__ solver into a file.
	 * The user can overload this function to store solutions in any other appropriate way.
	 * @param[in] filename name of the file to store solutions.
         * @throws CFileException.
	 */
	@Override
	public native void printSolution(String filename) throws IOException;

    /**
     * The function destroys the CLP object which reference is stored in `nativeHd`.
    */
	@Override
	public native void dispose();
}

