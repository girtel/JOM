/**
 * Copyright (c) 2015 Pablo Pavon Mari�o.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * <p>
 * Contributors:
 * Pablo Pavon Mari�o - initial API and implementation
 */

/**
 *
 */
package com.jom;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntFactory2D;
import cern.colt.matrix.tint.IntMatrix2D;
import com.sun.jna.Native;
import com.sun.jna.Pointer;

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

class _SOLVERWRAPPER_IPOPT
{
	final static int returnCode_Diverging_Iterates                 = 4;
	final static int returnCode_Error_In_Step_Computation          = -3;
	final static int returnCode_Feasible_Point_Found               = 6;
	final static int returnCode_Infeasible_Problem_Detected        = 2;
	final static int returnCode_Insufficient_Memory                = -102;
	final static int returnCode_Internal_Error                     = -199;
	final static int returnCode_Invalid_Number_Detected            = -13;
	final static int returnCode_Invalid_Option                     = -12;
	final static int returnCode_Invalid_Problem_Definition         = -11;
	final static int returnCode_Maximum_CpuTime_Exceeded           = -4;
	final static int returnCode_Maximum_Iterations_Exceeded        = -1;
	final static int returnCode_NonIpopt_Exception_Thrown          = -101;
	final static int returnCode_Not_Enough_Degrees_Of_Freedom      = -10;
	final static int returnCode_Restoration_Failed                 = -2;
	final static int returnCode_Search_Direction_Becomes_Too_Small = 3;
	final static int returnCode_Solve_Succeeded                    = 0;
	final static int returnCode_Solved_To_Acceptable_Level         = 1;
	final static int returnCode_Unrecoverable_Exception            = -100;
	final static int returnCode_User_Requested_Stop                = 5;
	private final HashMap<String, Object> param;
	private final _INTERNAL_SolverIO      s;
	private final String                  solverLibraryName;

	_SOLVERWRAPPER_IPOPT(_INTERNAL_SolverIO s, HashMap<String, Object> param)
	{
		this.s = s;
		this.solverLibraryName = (String) param.get("solverLibraryName");
		this.param = param;
	}

	static String errorMessage(int code)
	{
		switch (code)
		{
			case 0:
				return "Solve_Succeeded. This message indicates that IPOPT found a (locally) optimal point within the desired tolerances";
			case 1:
				return "Solved_To_Acceptable_Level. This indicates that the algorithm did not converge to the ``desired'' tolerances, but that it "
						+ "was able to obtain a point satisfying the ``acceptable'' tolerance level as specified by acceptable-* options. This may "
						+ "happen if the desired tolerances are too small for the current problem";
			case 2:
				return "Infeasible_Problem_Detected. The restoration phase converged to a point that is a minimizer for the constraint violation (in"
						+ " the -norm), but is not feasible for the original problem. This indicates that the problem may be infeasible (or at least"
						+ " that the algorithm is stuck at a locally infeasible point). The returned point (the minimizer of the constraint "
						+ "violation) might help you to find which constraint is causing the problem. If you believe that the NLP is feasible, it "
						+ "might help to start the optimization from a different point.";
			case 3:
				return "Search_Direction_Becomes_Too_Small. This indicates that IPOPT is calculating very small step sizes and making very little "
						+ "progress. This could happen if the problem has been solved to the best numerical accuracy possible given the current "
						+ "scaling.";
			case 4:
				return "Diverging_Iterates. This message is printed if the max-norm of the iterates becomes larger than the value of the option "
						+ "diverging_iterates_tol. This can happen if the problem is unbounded below and the iterates are diverging. ";
			case 5:
				return "User_Requested_Stop. This message is printed if the user call-back method intermediate_callback returned false";
			case 6:
				return "Feasible_Point_Found. This message is printed if the problem is ``square'' (i.e., it has as many equality constraints as "
						+ "free variables) and IPOPT found a feasible point.";
			case -1:
				return "Maximum_Iterations_Exceeded. This indicates that IPOPT has exceeded the maximum number of iterations as specified by the "
						+ "option max_iter";
			case -2:
				return "Restoration_Failed. This indicates that the restoration phase failed to find a feasible point that was acceptable to the "
						+ "filter line search for the original problem. This could happen if the problem is highly degenerate, does not satisfy the "
						+ "constraint qualification, or if your NLP code provides incorrect derivative information";
			case -3:
				return "Error_In_Step_Computation. This message is printed if IPOPT is unable to compute a search direction, despite several "
						+ "attempts to modify the iteration matrix. Usually, the value of the regularization parameter then becomes too large. One "
						+ "situation where this can happen is when values in the Hessian are invalid (NaN or Inf). You can check whether this is "
						+ "true by using the check_derivatives_for_naninf option. ";
			case -4:
				return "Maximum_CpuTime_Exceeded. ";
			case -10:
				return "Not_Enough_Degrees_Of_Freedom. This indicates that your problem, as specified, has too few degrees of freedom. This can "
						+ "happen if you have too many equality constraints, or if you fix too many variables (IPOPT removes fixed variables). ";
			case -11:
				return "Invalid_Problem_Definition. This indicates that there was an exception of some sort when building the IpoptProblem structure"
						+ " in the C or Fortran interface. Likely there is an error in your model or the main routine. ";
			case -12:
				return "Invalid_Option. This indicates that there was some problem specifying the options. See the specific message for details.";
			case -13:
				return "Invalid_Number_Detected. ";
			case -100:
				return "Unrecoverable_Exception. This indicates that IPOPT has thrown an exception that does not have an internal return code. See "
						+ "the specific message for details. ";
			case -101:
				return "NonIpopt_Exception_Thrown. An unknown exception was caught in IPOPT. This exception could have originated from your model or"
						+ " any linked in third party code. ";
			case -102:
				return "Insufficient_Memory. An error occurred while trying to allocate memory. The problem may be too large for your current memory"
						+ " and swap configuration. ";
			case -199:
				return "Internal_Error. An unknown internal error has occurred. Please notify the authors of IPOPT via the mailing list";
			default:
				throw new JOMException("JOM - IPOPT interface. Unknown IPOPT error message");
		}
	}

	void solve()
	{
		_JNA_IPOPT g = (_JNA_IPOPT) Native.loadLibrary(solverLibraryName, _JNA_IPOPT.class);

		final int n = this.s.in.numDecVariables;
		final double[] x_L = s.in.primalSolutionLowerBound.toArray();
		final double[] x_U = s.in.primalSolutionUpperBound.toArray();
		final int m = s.in.numLpConstraints + s.in.numNlpConstraints;
		final double[] g_L = s.in.constraintLowerBound.toArray();
		final double[] g_U = s.in.constraintUpperBound.toArray();

		LinkedHashMap<Integer, HashSet<Integer>> activeVarIds = (m == 0) ?
				new LinkedHashMap<Integer, HashSet<Integer>>() :
				s.in.lhsMinusRhsAccumulatedConstraint.getActiveVarIds();

		int numNonZerosActiveVarIdsMatrix = 0;
		if (m > 0) for (HashSet<Integer> colIds : activeVarIds.values()) numNonZerosActiveVarIdsMatrix += colIds.size();
		final int nele_jac = (m == 0) ? 0 : numNonZerosActiveVarIdsMatrix;
		final int numActiveCoordinates = nele_jac;
		final int nele_hess = 1; // we will use automatic hessian calculation
		final int index_style = 0; // C style

		int[] aux_rowsArrayActiveVarIds = new int[numNonZerosActiveVarIdsMatrix];
		int[] aux_colsArrayActiveVarIds = new int[numNonZerosActiveVarIdsMatrix];
		int counter = 0;
		for (Entry<Integer, HashSet<Integer>> e : activeVarIds.entrySet())
			for (Integer dv : e.getValue())
			{
				aux_rowsArrayActiveVarIds[counter] = e.getKey();
				aux_colsArrayActiveVarIds[counter++] = dv;
			}
		final int[] rowsArrayActiveVarIds = (s.in.lhsMinusRhsAccumulatedConstraint == null) ? new int[0] : aux_rowsArrayActiveVarIds;
		final int[] colsArrayActiveVarIds = (s.in.lhsMinusRhsAccumulatedConstraint == null) ? new int[0] : aux_colsArrayActiveVarIds;
		if ((rowsArrayActiveVarIds.length != numActiveCoordinates) || (colsArrayActiveVarIds.length != numActiveCoordinates))
			throw new JOMException("JOM - IPOPT interface. Unexpected error");
		IntMatrix2D coord = IntFactory2D.dense.make(nele_jac, 2);
		if (s.in.lhsMinusRhsAccumulatedConstraint != null)
		{
			coord.viewColumn(0).assign(rowsArrayActiveVarIds);
			coord.viewColumn(1).assign(colsArrayActiveVarIds);
		}
		final int[] indexesOfGradientActiveCoordinates = (s.in.lhsMinusRhsAccumulatedConstraint == null) ?
				new int[0] :
				DoubleMatrixND.sub2ind(coord, IntFactory1D.dense.make(new int[]{s.in.numConstraints, s.in.numDecVariables})).toArray();

		_JNA_IPOPT_CallBack_Eval_F eval_f = new _JNA_IPOPT_CallBack_Eval_F()
		{
			@Override
			public boolean callback(int n, Pointer x, boolean new_x, Pointer obj_value, Pointer user_data)
			{
				double[] values = x.getDoubleArray(0, n);
				obj_value.write(0, new double[]{s.in.objectiveFunction.evaluate_internal(values).get(0)}, 0, 1);
				return true;
			}
		};
		_JNA_IPOPT_CallBack_Eval_G eval_g = new _JNA_IPOPT_CallBack_Eval_G()
		{
			@Override
			public boolean callback(int n, Pointer x, boolean new_x, int m, Pointer g, Pointer user_data)
			{
				double[] values = x.getDoubleArray(0, n);
				double[] gOut = s.in.lhsMinusRhsAccumulatedConstraint.evaluate_internal(values).elements().toArray();
				if (gOut.length != m) throw new JOMException("JOM - IPOPT interface. Unexpected error");
				g.write(0, gOut, 0, m);
				return true;
			}
		};
		_JNA_IPOPT_CallBack_Eval_Grad_F eval_grad_f = new _JNA_IPOPT_CallBack_Eval_Grad_F()
		{
			@Override
			public boolean callback(int n, Pointer x, boolean new_x, Pointer grad_f, Pointer user_data)
			{
				double[] values = x.getDoubleArray(0, n);
				grad_f.write(0, s.in.objectiveFunction.evaluateJacobian_internal(values).viewRow(0).toArray(), 0, n);
				return true;
			}
		};
		_JNA_IPOPT_CallBack_Eval_Jac_G eval_jac_g = new _JNA_IPOPT_CallBack_Eval_Jac_G()
		{
			@Override
			public boolean callback(int n, Pointer x, boolean new_x, int m, int nele_jac, Pointer iRow, Pointer jCol, Pointer values, Pointer
					user_data)
			{
				if (values == null)
				{
					/* Fill in the sparsity pattern --> fill in iRow and jCol arrays */
					iRow.write(0, rowsArrayActiveVarIds, 0, nele_jac);
					jCol.write(0, colsArrayActiveVarIds, 0, nele_jac);
					return true;
				} else
				{
					/* Return the values */
					double[] xValues = x.getDoubleArray(0, n);
					DoubleMatrix2D res = s.in.lhsMinusRhsAccumulatedConstraint.evaluateJacobian_internal(xValues);
					double[] resSelectedIndexes = res.vectorize().viewSelection(indexesOfGradientActiveCoordinates).toArray();
					values.write(0, resSelectedIndexes, 0, nele_jac);
					return true;
				}
			}
		};
		_JNA_IPOPT_CallBack_Eval_H eval_h = new _JNA_IPOPT_CallBack_Eval_H()
		{
			@Override
			public boolean callback(int n, Pointer x, boolean new_x, double obj_factor, int m, Pointer lambda, boolean new_lambda, int nele_hess,
					Pointer iRow, Pointer jCol, Pointer values, Pointer user_data)
			{
				// return false;
				if (!s.isLinearProblem()) // hessians are computed automatically
					return false;

				if (values == null) // problem is linear, hessian is zero
				{
					int[] auxVar = new int[nele_hess];
					java.util.Arrays.fill(auxVar, 0);
					/* Fill in the sparsity pattern --> fill in iRow and jCol arrays */
					iRow.write(0, auxVar, 0, nele_hess);
					jCol.write(0, auxVar, 0, nele_hess);
					return true;
				} else
				{
					/* Return the values */
					double[] auxVar = new double[nele_hess];
					java.util.Arrays.fill(auxVar, 0);
					values.write(0, auxVar, 0, nele_hess);
					return true;
				}
			}
		};

		Pointer nlp = g.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style, eval_f, eval_g, eval_grad_f, eval_jac_g,
				eval_h);

		if (!s.in.toMinimize) g.AddIpoptNumOption(nlp, "obj_scaling_factor", -1);
		
		/* Set specific parameters specified by the user */
		for (Entry<String, Object> entry : this.param.entrySet())
		{
			String key = entry.getKey();
			if (key.equalsIgnoreCase("solverLibraryName")) continue; // a non-IPOPT specific option
			if (key.equals("maxSolverTimeInSeconds")) 
			{
				final Double val_maxSolverTimeInSeconds = ((Number) entry.getValue()).doubleValue();
				if (val_maxSolverTimeInSeconds > 0)
					g.AddIpoptNumOption(nlp, "max_cpu_time", val_maxSolverTimeInSeconds);
				continue;
			}
			Object value = entry.getValue();
			if (value instanceof String)
				g.AddIpoptStrOption(nlp, key, (String) value);
			else if (value instanceof Integer)
				g.AddIpoptIntOption(nlp, key, ((Integer) value).intValue());
			else if (value instanceof Double)
				g.AddIpoptNumOption(nlp, key, ((Double) value).doubleValue());
			else
				throw new JOMException("JOM - IPOPT interface. Unknown value type in parameters");
		}

		if (s.in.lhsMinusRhsAccumulatedConstraint != null)
			if (s.in.lhsMinusRhsAccumulatedConstraint.isLinear())
			{
				g.AddIpoptStrOption(nlp, "jac_c_constant", "yes"); // equality constraints
				// are linear => the
				// jacobian is called
				// once
				g.AddIpoptStrOption(nlp, "jac_d_constant", "yes"); // inequality
				// constraints are
				// linear => the
				// jacobian is called
				// once
			}
		if (s.isLinearProblem())
		{
			g.AddIpoptStrOption(nlp, "hessian_approximation", "exact"); // the hessian
			// is given
			// exactly,
			// and will be
			// equal to 0
			g.AddIpoptStrOption(nlp, "hessian_constant", "yes"); // the hessian is
			// constant
			// (objective
			// function and
			// constraints are
			// quadratic or
			// linear) => the
			// hessian is called
			// once
		} else
		{
			g.AddIpoptStrOption(nlp, "hessian_approximation", "limited-memory");
			g.AddIpoptStrOption(nlp, "hessian_constant", "no");
		}

		double[] x0 = s.in.primalInitialSolution.toArray();
		double[] gArray = new double[m];
		double[] obj_val = new double[1];
		double[] mult_g = new double[m];
		double[] mult_x_L = new double[n];
		double[] mult_x_U = new double[n];

		int status = g.IpoptSolve(nlp, x0, gArray, obj_val, mult_g, mult_x_L, mult_x_U, Pointer.NULL);

		s.problemAlreadyAttemptedTobeSolved = true;

		if (status == 0)
		{ // Solve_Succeeded = 0
		}

		s.out.statusCode = status;
		s.out.statusMessage = errorMessage(status);
		s.out.primalSolution = DoubleFactory1D.dense.make(x0);
		s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(mult_x_L);
		s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(mult_x_U);
		s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(gArray);
		s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(mult_g);
		s.out.primalCost = obj_val[0];
		s.out.bestOptimalityBound = s.in.toMinimize? -Double.MAX_VALUE : Double.MAX_VALUE; //p.getDblAttrib(XPRSconstants.BESTBOUND); // pablo: this may not be the one for this

		switch (status)
		{
			case returnCode_Solve_Succeeded:
			case returnCode_Solved_To_Acceptable_Level:
			case returnCode_Feasible_Point_Found:
				s.out.solutionIsFeasible = true;
				s.out.solutionIsOptimal = true;
				s.out.feasibleSolutionDoesNotExist = false;
				s.out.foundUnboundedSolution = false;
				break;
			case returnCode_Infeasible_Problem_Detected:
				s.out.solutionIsFeasible = false;
				s.out.solutionIsOptimal = false;
				s.out.feasibleSolutionDoesNotExist = true;
				s.out.foundUnboundedSolution = false;
				break;
			case returnCode_Search_Direction_Becomes_Too_Small:
			case returnCode_Maximum_Iterations_Exceeded:
				s.out.solutionIsFeasible = true;
				s.out.solutionIsOptimal = false;
				s.out.feasibleSolutionDoesNotExist = false;
				s.out.foundUnboundedSolution = false;
				break;
			case returnCode_Diverging_Iterates:
				s.out.solutionIsFeasible = true;
				s.out.solutionIsOptimal = false;
				s.out.feasibleSolutionDoesNotExist = false;
				s.out.foundUnboundedSolution = true;
				break;
			case returnCode_User_Requested_Stop:
			case returnCode_Restoration_Failed:
			case returnCode_Error_In_Step_Computation:
			case returnCode_Maximum_CpuTime_Exceeded:
			case returnCode_Not_Enough_Degrees_Of_Freedom:
			case returnCode_Invalid_Problem_Definition:
			case returnCode_Invalid_Option:
			case returnCode_Invalid_Number_Detected:
			case returnCode_Unrecoverable_Exception:
			case returnCode_NonIpopt_Exception_Thrown:
			case returnCode_Insufficient_Memory:
			case returnCode_Internal_Error:
				s.out.solutionIsFeasible = false;
				s.out.solutionIsOptimal = false;
				s.out.feasibleSolutionDoesNotExist = false;
				s.out.foundUnboundedSolution = false;
				break;
		}
	}

}
