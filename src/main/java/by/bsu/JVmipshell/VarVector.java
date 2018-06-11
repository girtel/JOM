/**
 * @file VarVector.java  Interface for the `VarVector` class.
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

import by.bsu.JVmipcl.*;
import by.bsu.JVmipshell.Var;
import by.bsu.JVmipshell.LPshell;

/**
 * Class `VarVector`, usually, represents an array of variables of identical type.
 */
public class VarVector
{
    /**
     * Dimention of this vector.
     */
	int dim;

    /**
     * Array of sizes. 
     */
	int[] sizes;

    /**
     * Array of vector variables. 
     */
	Var[] vars;

    /**
     * Return the dimention of this vector. 
     */
	public int getdim()
	{
		return dim;
	}

    /**
     * Return the sizes of this vector. 
     */
	public int[] getsizes()
	{
		return sizes;
	}

    /**
     * Return a reference to vars. 
     */
	public Var[] getvars()
	{
		return vars;
	}

    /**
     * The constructor creates a vector (list) of variables.
     *
     * @param[in] name vector name that is used when the values of variables are printed;
     * @param[in] type type of vector variables (`REAL`, `INT`, or `BIN`),
     *                 all vector variables are of the same type;
     * @param[in] lb common lower bound of all vector variables;
     * @param[in] ub common upper bound of all vector variables;
     * @param[in] sizes list of integers, the size of this list is the dimention of the vector.
     * @return reference to the created list of variables.
     * @warning maximal dimension of a vector is 4!
     */
	private void CreateVarVector(LPshell prob, String name, int type, double lb, double ub, int[] sizes) throws OutOfMemoryError
	{
		this.sizes = sizes;
		dim = sizes.length;
		int totalSize = 1;
		for (int sz : sizes) {
			totalSize *= sz;
		}
		vars = new Var[totalSize];

		int k = 0;
		if (dim == 1) {
			for (int i0=0; i0 < sizes[0]; ++i0)
				vars[k++] = new Var(prob,name + "(" + i0 + ")",type,lb,ub);
		}
		else if (dim == 2) {
			for (int i0=0; i0 < sizes[0]; ++i0) {
				for (int i1=0; i1 < sizes[1]; ++i1)
					vars[k++] = new Var(prob,name + "(" + i0 + "," + i1 + ")",type,lb,ub);
			}
		}
		else if (dim ==  3) {
			for (int i0=0; i0 < sizes[0]; ++i0) {
				for (int i1=0; i1 < sizes[1]; ++i1) {
					for (int i2=0; i2 < sizes[2]; ++i2)
						vars[k++] = new Var(prob,name + "(" + i0 + "," + i1 + "," + i2 + ")",
										type,lb,ub);
				}
			}
		}
		else if (dim == 4) {
			for (int i0=0; i0 < sizes[0]; ++i0) {
				for (int i1=0; i1 < sizes[1]; ++i1) {
					for (int i2=0; i2 < sizes[2]; ++i2) {
						for (int i3=0; i3 < sizes[3]; ++i3)
							vars[k++] = new Var(prob,name + "(" + i0 + "," + i1 + "," + i2 + "," + i3 + ")",
										type,lb,ub);
					}
				}
			}
		}
	}

    /**
     * General constructor.
     *
     * @param[in] prob problem to which all variables of `this` vector be added;
     * @param[in] name this name appended with the index of a variable from `this` vector is used when printing that variable;
     * @param[in] type common type of all variables in `this` vector;
     * @param[in] lb common lower bound of all variables in `this` vector;
     * @param[in] ub common upper bound of all variables in `this` vector;
     * @param[in] sizes array of sizes.
     */
	public VarVector(LPshell prob, String name, int type, double lb, double ub, int... sizes) throws OutOfMemoryError
	{
		CreateVarVector(prob,name,type,lb,ub,sizes);
	}

    /**
     * The constructor useful for declaring binary variables or non-negative real and integer variables.
     *
     * @param[in] prob problem to which all variables of `this` vector be added;
     * @param[in] name this name appended with the index of a variable from `this` vector is used when printing that variable;
     * @param[in] type common type of all variables in `this` vector;
     * @param[in] sizes array of sizes.
     */
	public VarVector(LPshell prob, String name, int type, int... sizes) throws OutOfMemoryError
	{
		CreateVarVector(prob,name,type,0.0,((type & MIP.VAR_BIN) != 0)? 1.0 : LP.VAR_INF,sizes);
	}

    /**
     * @param[in] ind array (tuple) of integers.
     * @return vector variable indexed by the input tuple.
     * @throws `IndexOutOfBoundsException` if the input tuple is not a valid index.
     */
	public Var get(int... ind)
	{
		boolean flag = true;
		int k = 0;
		if (ind.length == dim) {
			for (int i=0; i < dim; ++i) {
				if (ind[i] < 0 || ind[i] >= sizes[0]) {
					flag = false;
					break;
				}
			}
		}
		else
			flag = false;
		if (flag) {
			int d = dim - 1;
			k = ind[d];
			if (d > 0) {
				int sz = sizes[d];
				k += (ind[--d] * sz);
				if (d > 0) {
					sz *= sizes[d];
					k += (ind[--d] * sz);
					if (d > 0) {
						sz *= sizes[d];
						k += (ind[--d] * sz);
					}
				}
			}
		}
		else {
			throw new IndexOutOfBoundsException("Var.get(): bad index " + ind.toString());
		}

		return vars[k];
	}
}


