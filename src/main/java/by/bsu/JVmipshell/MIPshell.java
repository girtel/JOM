/**
 * @file MIPshell.java  Interface for the `MIP` class.
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
import by.bsu.JVmipshell.VarVector;
import by.bsu.JVmipshell.LinSum;
import by.bsu.JVmipshell.Ctr;
import by.bsu.JVmipshell.Function;

/**
 * This class connects all components --- variables, constraints, objective --- together to build a MIP Problem.
 */
public class MIPshell extends LPshell
{
    /**
     * In __MIPCL__ this flag indicates that a variables takes only integer values. 
     */
	public static final int VAR_BIN = MIP.VAR_BIN;

    /**
     * In __MIPCL__ this flag indicates that a variables is _binary_, i.e., it takes only two values, `0` or `1`. 
     */
	public static final int VAR_INT = MIP.VAR_INT;

    /**
     * List of problem function constraints.
     */
	List<Function> funcs;

    /**
     * Math program (reference to a MIP object).
     */
	MIP mp;

    /**
     * This flag is set to `true` only if an optimal solution has been found.
     */
	boolean is_solutionOptimal;

    /**
     * This flag is set to `true` if solution procedure stopped after exceeding given time limit.
     */
	boolean timeLimitStop;
	
    /**
     *  The constructor initializes an empty (without variables, constraints, and objective) problem.
     *
     * @param[in] name problem name.
     */
	public MIPshell(String name)
	{
		super(name);
		timeLimitStop = false;
		is_solutionOptimal = false;
		funcs = null;
	}

    /**
     * The procedure adds a function constraint to the problem list of function constraints.
     *
     * @param[in] func function constraint to be added.
     * @return `true` in case of success.
     */
	public boolean addFunc(Function func)
	{
//		isPureLP = false;
		if (funcs == null)
			funcs = new ArrayList<Function>();
		func.sethd(funcs.size());
		return funcs.add(func);
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
     * The function initializes either a CMIP-instance for the MIP problem.
     */
	@Override
	public void mipclModel()
	{
		int n = vars.size();
		int m = ctrs.size();

        	int nz = 0;        
		for (Ctr ctr : ctrs) {
			nz += ctr.lsum.getterms().size();
		}
		
		int fn=0, funcSize=0;
		if (funcs != null) {
			fn = funcs.size();
			for (Function func : funcs) {
				funcSize += func.getpoints().size();
			}
		}

		mp = new MIP(name);

        	mp.openMatrix(m+3*fn,n+funcSize,nz+3*funcSize+2*fn); //,true,false,0,0,0);

		for (Var var : vars) {
			int j = mp.addVar(var.gethd(),var.gettype(),0.0,var.getlb(),var.getub());
            		var.setval(0.0); // is used when adding constraints to handle terms with the same variables
			if (var.getpriority() != 0)
				mp.setVarPriority(j,var.getpriority());
		}

		if (funcs != null) {		
			int hd = n;
			for (Function func : funcs) {
				int k = func.getpoints().size();
				for (int i=0; i < k; ++i) {
					mp.addVar(hd,0,0.0,0.0,1.0);
				}
				++hd ;
			}
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

		if (funcs != null) {
// model function constraints; each such constraint is represented by 3 equations
			for (Function func : funcs) {
				mp.addCtr(row,0,0.0,0.0);
				mp.addEntry(-1.0,row,func.getx().gethd());
				int hd = n;
				List<Point> points = func.getpoints();
				for (Point point : points) {
					mp.addEntry(point.getx(),row,hd);
					++hd;
				}
				++row;
				mp.addCtr(row,0,0.0,0.0);
				mp.addEntry(-1.0,row,func.y.hd);
				hd = n;
				for (Point point : points) {
					mp.addEntry(point.gety(),row,hd);
					++hd;
				}
				++row;
				mp.addCtr(row,0x00010000,1.0,1.0); // 0x00010000 is SOS2 flag in CMIP
				for (int j=n; j < hd; ++j) {
					mp.addEntry(1.0,row,j);
				}
				++row;
				n = hd;
			}
		}
            
		mp.closeMatrix();
	} // end of mipclModel()

    /**
     *  The procedure starts solving the problem.
     *
     * @param[in] silent  if set to `true`, MIPCL info-messages will not be displayed;
     * @param[in] timeLimit  time limit (in seconds) for the solver.
     */
	public void solve(boolean silent, long timeLimit)
	{
		mp.beSilent(silent);
		mp.optimize(timeLimit,0.0,"");
		timeLimitStop = mp.timeLimitStop();
		problemSolved = true;
		if (mp.isSolution()) {
			is_solution = true;
			is_solutionOptimal = mp.isSolutionOptimal();
			objVal = mp.getObjVal();
			mp.getSolution();
			int[] hd = mp.getVarHandles();
			double[] x = mp.getX();
			int n = hd.length;
			for (int i=0; i < n; ++i) {
				vars.get(hd[i]).setval(x[i]);
			}
		}
		else {
			is_solution = false;
			if (!timeLimitStop && mp.isInfeasible())
				is_infeasible = true;
		}
	} // end of solve()

    /**
     * The procedure first calls `mipclModel()` to build the model,
     * and then calls `solve()` to solve that model.
     *
     * @param[in] silent  if set to `true`, MIPCL info-messages will not be displayed;
     * @param[in] timeLimit  time limit (in seconds) for the solver.
     */
	public void optimize(boolean silent, long timeLimit)
	{
		mipclModel();
		solve(silent, timeLimit);
	}

    /**
     * The procedure just calls `optimize(boolean, long)` setting its parameters to default values:
     *      - `silent` to `true` (verbose mode);
     *      - timeLimit to `10000000l`.
     */
	@Override
	public void optimize()
	{
		optimize(false, 10000000l);
	}

    /**
     * @return true if the solution found by the MIPCL solver is optimal.
     */
	public boolean isSolutionOptimal()
	{return is_solutionOptimal;}

    /**
     * The procedure prints the solution.
     *
     * When an LP has no feasible solution, the procudure prints a certificate of inconsistency.
     */
	@Override
	public void printSolution()
	{
		if (!problemSolved) {
			System.out.println("Please call optimize() first");
			return;
		}

		if (is_solution) {
			System.out.format("Optimal objective value is %.4f%n",objVal);
			for (Var var : vars) {
				System.out.format("%s = %.4f%n",var.getname(),var.getval());
			}
		}
		else // MIP
			System.out.println("Optimal solution has not been found");
	} // end of printSolution()

    /**
     * The procedure writes the solution to a file.
     *
     * @param[in] solfile file name.
     * @throws  IOException
     */
	@Override
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
		}
		else // MIP
			fout.format("Optimal solution has not been found%n");
	}
}


