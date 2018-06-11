/**
 * @file LPshell.java  Interface for the `LP` class'.
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

import java.io.IOException;
import java.lang.IndexOutOfBoundsException;
import java.lang.OutOfMemoryError;
import java.util.Formatter;
import java.util.Arrays;
import java.util.List;
import java.util.ArrayList;

import by.bsu.JVmipcl.*;
import by.bsu.JVmipshell.Var;
import by.bsu.JVmipshell.Ctr;
import by.bsu.JVmipshell.LinSum;



/**
 * This class connects all components --- variables, constraints, objective --- together to build an LP Problem.
 */
public class LPshell
{
    /**
     * Infinity value.
     */
	public static final double INF = LP.INF;

    /**
     * Infinity value for variables.
     */
	public static final double VAR_INF = LP.VAR_INF;

    /**
     * Value of \"zero\", i.e., cooefficients that less than `ZERO` are treated as zeroes.
     */
	public static final double ZERO = 1.0e-12;

    /**
     * Problem name.
     */
	String name;

    /**
     * If `True`, the objective is maximized; otherwise, the objective is minimized.
     */
	boolean sense;

    /**
     * Reference to a `LinSum` object that represents the objective function.
     */
	LinSum obj;

    /**
     * If `true`, the problem has been solved.
     */
	boolean problemSolved;

    /**
     * When the problem has been solved, this flag is set to `true` if a feasible solution has been found;
       otherwise, the flag is set to `false`.
     */
	boolean is_solution;

    /**
     * This flag is set to `true` only if the problem being solved is infeasible (has no solution).
     */
	boolean is_infeasible;

    /**
     * This flag is set to `true` only if the problem being solved is _unbounded_ (its objective tends to (+/-)infinity).
     */
	boolean is_unbounded;

    /**
     * List of problem variables.
     */
	List<Var> vars;

    /**
     * List of problem constraints.
     */
	List<Ctr> ctrs;

    /**
     * Objective value of the best solution found.
     */
	double objVal;

	LP mp;

    /**
     *  The constructor initializes an empty (without variables, constraints, and objective) problem.
     *
     * @param[in] name problem name.
     */
	public LPshell(String name)
	{
		mp = null;
		this.name = name;
		problemSolved = false;
		sense = true;
		obj = null;
		is_solution = false;
		is_infeasible = false;
		is_unbounded = false;
		vars = new ArrayList<Var>();
		ctrs = new ArrayList<Ctr>();
	}

    /**
     * The function adds a variable to the problem list of variables.
     *
     * @param[in] var variable to be added.
     * @return `true` in case of success.
     */
	public boolean addVar(Var var)
	{
		var.sethd(this.vars.size());
        	return this.vars.add(var);
	}

    /**
     * The function adds a constraint to the problem list of constraints.
     *
     * If the input `LinSum` object is empty (contains no term), then
     * the function returns `null` indicating that no new constraint has been created.
     *
     * @param[in] lsum linear sum object;
     * @param[in] lhs left hand side;
     * @param[in] rhs right hand side.
     * @return reference to newly created constraint, or `null` if lsum is empty.
     */
	public Ctr addCtr(LinSum lsum, double lhs, double rhs)
	{
		if (lsum.getterms() != null) {
			Ctr ctr = new Ctr(ctrs.size(),lsum,lhs,rhs);
			ctrs.add(ctr);
			return ctr;
		}
		else
			return null;
	}

    /**
     * Set the objective to be minimize a linear sum.
     *
     * @param[in] lsum objective linear function.
     */
	public void minimize(LinSum lsum)
	{
		sense = false;
		obj = lsum;
	}

    /**
     * Set the objective to minimize a variable.
     *
     * @param[in] var objective variable.
     */
	public void minimize(Var var)
	{
		sense = false;
		obj = new LinSum();
		obj.add(1.0,var);
	}

    /**
     * Set the objective to be maximize a linear sum.
     *
     * @param[in] lsum objective linear function.
     */
	public void maximize(LinSum lsum)
	{
		sense = true;
		obj = lsum;
	}

    /**
     * Set the objective to maximize a variable.
     *
     * @param[in] var objective variable.
     */
	public void maximize(Var var)
	{
		sense = true;
		obj = new LinSum();
		obj.add(1.0,var);
	}

    /**
     * @return `true` if a solution has been found.
     */
	public boolean isSolution()
 	{return is_solution;}

    /**
     * @return `true` if this LP is _infeasible_ (has no solution).
     */
	public boolean isInfeasible()
 	{return is_infeasible;}

    /**
     * @return `true` if this LP is _unbounded_ (its objective tends to (+/-)infinity).
     */
	public boolean isUnbounded()
 	{return is_unbounded;}

    /**
     * The function initializes either an CLP-instance for the LP.
     */
	public void mipclModel()
	{
		int n = vars.size();
		int m = ctrs.size();

        	int nz = 0;        
		for (Ctr ctr : ctrs) {
			nz += ctr.lsum.getterms().size();
		}

		mp = new LP(name);

		mp.preprocOff();

        	mp.openMatrix(m,n,nz); //,true,false,0,0,0);

		for (Var var : vars) {
			int j = mp.addVar(var.gethd(),var.gettype(),0.0,var.getlb(),var.getub());
            		var.setval(0.0); // is used when adding constraints to handle terms with the same variables
		}

		mp.setObjSense(sense);
		for (Term t : obj.getterms()) {
			mp.setObjCoeff(t.getvar().gethd(),t.getcoeff());
		}

		int row = 0;
		for (Ctr ctr : ctrs) {
			LinSum lsum = ctr.getlsum();
			if (lsum != null) {
				List<Term> terms = lsum.getterms();
				double v = 0.0;
				for (Term t : terms) {
					if (t.getvar() == null)
						v += t.getcoeff();
					else
						t.getvar().incval(t.coeff);
				}
				double lhs = ctr.getlhs();
				double rhs = ctr.getrhs();
				if (v > ZERO || v < -ZERO) {
					if (lhs > -INF+1.0)
						lhs -= v;
					if (rhs < INF-1.0)
						rhs -= v;
				}
				mp.addCtr(row,0,lhs,rhs);
				for (Term t : terms) {
					Var var = t.getvar();
					if (var != null) {
						double val = var.getval();
						if (val > ZERO || val < -ZERO) {
							mp.addEntry(val,row,var.gethd());
                            				var.setval(0.0);
						}
					}
				}
				++row;
			}
		}
            
		mp.closeMatrix();
	} // end of mipclModel()

    /**
     *  The procedure starts solving the problem.
     *
     * @param[in] silent  if set to `true`, no MIPCL info-messages are displayed;
     */
	public void solve(boolean silent)
	{
		mp.beSilent(silent);
		mp.optimize();
		problemSolved = true;
		if (mp.isSolution()) {
			is_solution = true;
			objVal = mp.getObjVal();
			mp.getSolution(); // LP asks CLP for an optimal solution
			double[] y = mp.getShadowPrices(); // LPshell gets shadow prices from LP
			int i = 0;
			for (Ctr ctr : ctrs) {
				ctr.setprice(y[i++]);
			}
			double[] x = mp.getX(); // LPshell gets solution vector from LP
			i = 0;
			for (Var var : vars) {
				var.setval(x[i]);
			}
		}
		else {
			is_solution = false;
			if (mp.isLpInfeasible()) {
				is_infeasible = true;
				mp.getInconsistencyCertificate(); // LP asks CLP for an inconsistency certificate
				double[] y = mp.getShadowPrices(); // LPshell gets row multiplipyers from LP
				int i = 0;
				for (Ctr ctr : ctrs) {
					ctr.setprice(y[i++]);
				}
				double[] x = mp.getReducedCosts(); // LPshell gets column multiplipyers from LP
				i = 0;
				for (Var var : vars) {
					var.setval(x[i]);
				}
			}
			else if (mp.isLpUnbounded())
				is_unbounded = true;
		}
	} // end of solve()

    /**
     * The procedure first calls `mipclModel()` to build the model,
     * and then calls `solve()` to solve that model.
     *
     * @param[in] silent  if set to `true`, MIPCL info-messages will not be displayed;
     */
	public void optimize(boolean silent)
	{
		mipclModel();
		solve(silent);
	}

    /**
     * The procedure just calls `optimize(boolean)` setting its only parameter to `true` (verbose mode).
     */
	public void optimize()
	{
		optimize(true);
	}


    /**
     * Get the objective value.
     *
     * @return the objective value of the solution found.
     */
	public double getobjVal()
	{
        	return objVal;
	}

    /**
     * The procedure prints the problem (useful for debugging).
     */
	public void _print()
	{
		System.out.println(((sense)? "max ": "min ") + obj.toString());
		for (Ctr ctr : ctrs) {
            		System.out.println(ctr);
		}
		for (Var var : vars) {
			System.out.println(var);
		}
	}

    /**
     * The procedure prints the solution.
     *
     * When an LP has no feasible solution, the procudure prints a certificate of inconsistency.
     */
	public void printSolution()
	{
		if (!problemSolved) {
			System.out.println("Please call optimize() first");
			return;
		}

		if (is_solution) {
			System.out.format("Optimal objective value is %.4f%n}",objVal);
			for (Var var : vars) {
				System.out.format("%s = %.4f%n",var.getname(),var.getval());
			}
			System.out.println("=== Constraint shadow prices (dual solution):"); 
			for (Ctr ctr : ctrs) {
				System.out.format("%.4f : %s%n",ctr.getprice(),ctr);
			}
		}
		else { // infeasible
			System.out.println("=== Problem constraints are contradictory.");
			System.out.println("Proof of inconsistency:");
			double q = 0.0;
			for (Ctr ctr : ctrs) {
				double price = ctr.getprice();
				if (price > 1.0e-6) {
					System.out.format("%.4f x (%s)%n",price,ctr.toString(false,true));
					q += ctr.getrhs() * price;
				}
				else if (price < -1.0e-6) {
					System.out.format("%.4f x (%s)%n",price,ctr.toString(true,false));
					q += ctr.getlhs() * price;
				}
			}

			for (Var var : vars) {
				double val = var.getval();
				if (val > 1.0e-6) {
					System.out.format("%.4f x (%s <= %.4f)%n",val,var.getname(),var.getub());
					q += var.getub() * val;
				}
				else if (val < -1.0e-6) {
					System.out.format("%.4f x (%.4f <= %s)%n",val,var.getlb(),var.getname());
                               		q += var.getlb() * val;
				}
			}
			System.out.format("= (0 <= %.4f)%n",q);
		}
	} // end of printSolution()

    /**
     * The procedure writes the solution to a file.
     *
     * @param[in] solfile file name.
     * @throws  IOException
     */
	public void writeSolution(String solfile) throws IOException
	{
		if (!problemSolved) {
			System.out.println("Please call optimize() first");
			return;
		}

		Formatter fout = new Formatter(solfile);

		if (is_solution) {
			fout.format("Optimal objective value is %.4f%n",objVal);
			for (Var var : vars) {
				fout.format("%s = %.4f%n",var.getname(),var.getval());
			}
                        fout.format("=== Constraint shadow prices (dual solution):%n"); 
			for (Ctr ctr : ctrs) {
				fout.format("%.4f : %s",ctr.getprice(),ctr);
			}
		}
		else { // infeasible problem
			fout.format("=== Problem constraints are contradictory.%n");
			fout.format("Proof of inconsistency:%s");
			double q = 0.0;
			for (Ctr ctr : ctrs) {
				double price = ctr.getprice();
				if (ctr.getprice() > 1.0e-6) {
					fout.format("%.4f x (%s)%n",price,ctr.toString(false,true));
					q += ctr.getrhs() * price;
				}
				else if (price < -1.0e-6) {
					fout.format("%.4f x (%s)%n",price,ctr.toString(true,false));
					q += ctr.getlhs() * price;
				}
			}

			for (Var var : vars) {
				double val = var.getval();
				fout.format("%.4f x %n",val);
				if (val > 1.0e-6) {
					fout.format("%.4f x (%s <= %.4f)%n",val,var.getname(),var.getub());
					q += var.getub() * val;
				}
				else if (val < -1.0e-6) {
					fout.format("%.4f x (%.4f <= %s)%n",val,var.getlb(),var.getname());
					q += var.getlb() * val;
				}
			}

			fout.format("= (0 <= %.4f)%n",q);
		}
	}
}


