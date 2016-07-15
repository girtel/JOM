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

import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import org.coinor.Ipopt;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntFactory2D;
import cern.colt.matrix.tint.IntMatrix2D;

class _SOLVERWRAPPER_JIPOPT
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

	_SOLVERWRAPPER_JIPOPT(_INTERNAL_SolverIO s, HashMap<String, Object> param)
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
//		_JNA_IPOPT g = (_JNA_IPOPT) Native.loadLibrary(solverLibraryName, _JNA_IPOPT.class);
//
//		final int n = this.s.in.numDecVariables;
//		final double[] x_L = s.in.primalSolutionLowerBound.toArray();
//		final double[] x_U = s.in.primalSolutionUpperBound.toArray();
//		final int m = s.in.numLpConstraints + s.in.numNlpConstraints;
//		final double[] g_L = s.in.constraintLowerBound.toArray();
//		final double[] g_U = s.in.constraintUpperBound.toArray();
//
//		LinkedHashMap<Integer, HashSet<Integer>> activeVarIds = (m == 0) ?
//				new LinkedHashMap<Integer, HashSet<Integer>>() :
//				s.in.lhsMinusRhsAccumulatedConstraint.getActiveVarIds();
//
//		int numNonZerosActiveVarIdsMatrix = 0;
//		if (m > 0) for (HashSet<Integer> colIds : activeVarIds.values()) numNonZerosActiveVarIdsMatrix += colIds.size();
//		final int nele_jac = (m == 0) ? 0 : numNonZerosActiveVarIdsMatrix;
//		final int numActiveCoordinates = nele_jac;
//		final int nele_hess = 1; // we will use automatic hessian calculation
//		final int index_style = 0; // C style
//
//		int[] aux_rowsArrayActiveVarIds = new int[numNonZerosActiveVarIdsMatrix];
//		int[] aux_colsArrayActiveVarIds = new int[numNonZerosActiveVarIdsMatrix];
//		int counter = 0;
//		for (Entry<Integer, HashSet<Integer>> e : activeVarIds.entrySet())
//			for (Integer dv : e.getValue())
//			{
//				aux_rowsArrayActiveVarIds[counter] = e.getKey();
//				aux_colsArrayActiveVarIds[counter++] = dv;
//			}
//		final int[] rowsArrayActiveVarIds = (s.in.lhsMinusRhsAccumulatedConstraint == null) ? new int[0] : aux_rowsArrayActiveVarIds;
//		final int[] colsArrayActiveVarIds = (s.in.lhsMinusRhsAccumulatedConstraint == null) ? new int[0] : aux_colsArrayActiveVarIds;
//		if ((rowsArrayActiveVarIds.length != numActiveCoordinates) || (colsArrayActiveVarIds.length != numActiveCoordinates))
//			throw new JOMException("JOM - IPOPT interface. Unexpected error");
//		IntMatrix2D coord = IntFactory2D.dense.make(nele_jac, 2);
//		if (s.in.lhsMinusRhsAccumulatedConstraint != null)
//		{
//			coord.viewColumn(0).assign(rowsArrayActiveVarIds);
//			coord.viewColumn(1).assign(colsArrayActiveVarIds);
//		}
//		final int[] indexesOfGradientActiveCoordinates = (s.in.lhsMinusRhsAccumulatedConstraint == null) ?
//				new int[0] :
//				DoubleMatrixND.sub2ind(coord, IntFactory1D.dense.make(new int[]{s.in.numConstraints, s.in.numDecVariables})).toArray();
//		_JNA_IPOPT_CallBack_Eval_F eval_f = new _JNA_IPOPT_CallBack_Eval_F()
//		{
//			@Override
//			public boolean callback(int n, Pointer x, boolean new_x, Pointer obj_value, Pointer user_data)
//			{
//				double[] values = x.getDoubleArray(0, n);
//				obj_value.write(0, new double[]{s.in.objectiveFunction.evaluate_internal(values).get(0)}, 0, 1);
//				return true;
//			}
//		};
//		_JNA_IPOPT_CallBack_Eval_G eval_g = new _JNA_IPOPT_CallBack_Eval_G()
//		{
//			@Override
//			public boolean callback(int n, Pointer x, boolean new_x, int m, Pointer g, Pointer user_data)
//			{
//				double[] values = x.getDoubleArray(0, n);
//				double[] gOut = s.in.lhsMinusRhsAccumulatedConstraint.evaluate_internal(values).elements().toArray();
//				if (gOut.length != m) throw new JOMException("JOM - IPOPT interface. Unexpected error");
//				g.write(0, gOut, 0, m);
//				return true;
//			}
//		};
//		_JNA_IPOPT_CallBack_Eval_Grad_F eval_grad_f = new _JNA_IPOPT_CallBack_Eval_Grad_F()
//		{
//			@Override
//			public boolean callback(int n, Pointer x, boolean new_x, Pointer grad_f, Pointer user_data)
//			{
//				double[] values = x.getDoubleArray(0, n);
//				grad_f.write(0, s.in.objectiveFunction.evaluateJacobian_internal(values).viewRow(0).toArray(), 0, n);
//				return true;
//			}
//		};
//		_JNA_IPOPT_CallBack_Eval_Jac_G eval_jac_g = new _JNA_IPOPT_CallBack_Eval_Jac_G()
//		{
//			@Override
//			public boolean callback(int n, Pointer x, boolean new_x, int m, int nele_jac, Pointer iRow, Pointer jCol, Pointer values, Pointer
//					user_data)
//			{
//				if (values == null)
//				{
//					/* Fill in the sparsity pattern --> fill in iRow and jCol arrays */
//					iRow.write(0, rowsArrayActiveVarIds, 0, nele_jac);
//					jCol.write(0, colsArrayActiveVarIds, 0, nele_jac);
//					return true;
//				} else
//				{
//					/* Return the values */
//					double[] xValues = x.getDoubleArray(0, n);
//					DoubleMatrix2D res = s.in.lhsMinusRhsAccumulatedConstraint.evaluateJacobian_internal(xValues);
//					double[] resSelectedIndexes = res.vectorize().viewSelection(indexesOfGradientActiveCoordinates).toArray();
//					values.write(0, resSelectedIndexes, 0, nele_jac);
//					return true;
//				}
//			}
//		};
//		_JNA_IPOPT_CallBack_Eval_H eval_h = new _JNA_IPOPT_CallBack_Eval_H()
//		{
//			@Override
//			public boolean callback(int n, Pointer x, boolean new_x, double obj_factor, int m, Pointer lambda, boolean new_lambda, int nele_hess,
//					Pointer iRow, Pointer jCol, Pointer values, Pointer user_data)
//			{
//				// return false;
//				if (!s.isLinearProblem()) // hessians are computed automatically
//					return false;
//
//				if (values == null) // problem is linear, hessian is zero
//				{
//					int[] auxVar = new int[nele_hess];
//					java.util.Arrays.fill(auxVar, 0);
//					/* Fill in the sparsity pattern --> fill in iRow and jCol arrays */
//					iRow.write(0, auxVar, 0, nele_hess);
//					jCol.write(0, auxVar, 0, nele_hess);
//					return true;
//				} else
//				{
//					/* Return the values */
//					double[] auxVar = new double[nele_hess];
//					java.util.Arrays.fill(auxVar, 0);
//					values.write(0, auxVar, 0, nele_hess);
//					return true;
//				}
//			}
//		};

//		Pointer nlp = g.CreateIpoptProblem(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style, eval_f, eval_g, eval_grad_f, eval_jac_g,
//				eval_h);
//
		
		IPOPTSolverAccess nlp = new IPOPTSolverAccess();
		if (!s.in.toMinimize) 
			if (!nlp.addNumOption("obj_scaling_factor", -1)) throw new JOMException ("The scaling option to -1, to force a maximization problem, could not be set in IPOPT");
		
		/* Set specific parameters specified by the user */
		for (Entry<String, Object> entry : this.param.entrySet())
		{
			String key = entry.getKey();
			if (key.equalsIgnoreCase("solverLibraryName")) continue; // a non-IPOPT specific option
			Object value = entry.getValue();
			if (value instanceof String)
				if (!nlp.addStrOption(key, (String) value)) throw new JOMException ("The option " + key + " could not be set as an String in IPOPT solver");
			else if (value instanceof Integer)
				if (!nlp.addIntOption(key, (Integer) value)) throw new JOMException ("The option " + key + " could not be set as an Integer in IPOPT solver");
			else if (value instanceof Double)
				if (!nlp.addNumOption(key, (Double) value)) throw new JOMException ("The option " + key + " could not be set as a Double in IPOPT solver");
			else
				throw new JOMException("JOM - IPOPT interface. Unknown value type in parameters");
		}

		if (s.in.lhsMinusRhsAccumulatedConstraint != null)
			if (s.in.lhsMinusRhsAccumulatedConstraint.isLinear())
			{
				/* equality constraints */
				if (!nlp.addStrOption("jac_c_constant", "yes")) throw new JOMException ("The option " + "jac_c_constant" + " could not be set in IPOPT solver"); 
				/* inequality constraints */
				if (!nlp.addStrOption("jac_d_constant", "yes")) throw new JOMException ("The option " + "jac_d_constant" + " could not be set in IPOPT solver"); 
			}
		if (s.isLinearProblem())
		{
			/* the hessian is given exactly, and will be equal to 0 */
			if (!nlp.addStrOption("hessian_approximation", "exact")) throw new JOMException ("The option " + "hessian_approximation" + " could not be set in IPOPT solver"); 
			if (!nlp.addStrOption("hessian_constant", "yes")) throw new JOMException ("The option " + "hessian_constant" + " could not be set in IPOPT solver"); 
			/* the hessian is constant (objective function and constraints are quadratic or linear) => the hessian is called once */
		} else
		{
			if (!nlp.addStrOption("hessian_approximation", "limited-memory")) throw new JOMException ("The option " + "hessian_approximation" + " could not be set in IPOPT solver"); 
			if (!nlp.addStrOption("hessian_constant", "no")) throw new JOMException ("The option " + "hessian_constant" + " could not be set in IPOPT solver"); 
		}

		double[] x0 = s.in.primalInitialSolution.toArray();
//		double[] gArray = new double[m];
//		double[] obj_val = new double[1];
//		double[] mult_g = new double[m];
//		double[] mult_x_L = new double[n];
//		double[] mult_x_U = new double[n];

		int status = nlp.solve (x0); //IpoptSolve(nlp, x0, gArray, obj_val, mult_g, mult_x_L, mult_x_U, Pointer.NULL);

		s.problemAlreadyAttemptedTobeSolved = true;

		if (status == 0)
		{ // Solve_Succeeded = 0
		}

		s.out.statusCode = status;
		s.out.statusMessage = errorMessage(status);
		s.out.primalSolution = DoubleFactory1D.dense.make(x0);
		s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(nlp.getMultLowerBounds());
		s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(nlp.getMultUpperBounds());
		double [] gArray = new double [s.in.numConstraints];
		if (!nlp.eval_g(s.in.numDecVariables, x0, true, s.in.numConstraints, gArray)) throw new JOMException ("Could not retrieve the value of the constraints at the returned solution in IPOPT");
		s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(gArray);
		s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(nlp.getMultConstraints());
		s.out.primalCost = nlp.getObjVal();
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
		if (s.out.solutionIsOptimal) s.out.bestOptimalityBound = s.out.primalCost;
	}

	class IPOPTSolverAccess extends Ipopt
	{
		final int n;
		final double[] x_L;
		final double[] x_U;
		final int m;
		final double[] g_L;
		final double[] g_U;
		final int nele_jac;
		final int[] rowsArrayActiveVarIds;
		final int[] colsArrayActiveVarIds;
		final int[] indexesOfGradientActiveCoordinates;
		
		public IPOPTSolverAccess () 
		{
			super ("C:\\Dropbox\\eclipseWorkspaceMars1\\st-jom\\src\\main\\resources", "jipopt");
			n = s.in.numDecVariables;
			x_L = s.in.primalSolutionLowerBound.toArray();
			x_U = s.in.primalSolutionUpperBound.toArray();
			m = s.in.numLpConstraints + s.in.numNlpConstraints;
			g_L = s.in.constraintLowerBound.toArray();
			g_U = s.in.constraintUpperBound.toArray();

			LinkedHashMap<Integer, HashSet<Integer>> activeVarIds = (m == 0) ?
					new LinkedHashMap<Integer, HashSet<Integer>>() :
					s.in.lhsMinusRhsAccumulatedConstraint.getActiveVarIds();

			int numNonZerosActiveVarIdsMatrix = 0;
			if (m > 0) for (HashSet<Integer> colIds : activeVarIds.values()) numNonZerosActiveVarIdsMatrix += colIds.size();
			nele_jac = (m == 0) ? 0 : numNonZerosActiveVarIdsMatrix;
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
			rowsArrayActiveVarIds = (s.in.lhsMinusRhsAccumulatedConstraint == null) ? new int[0] : aux_rowsArrayActiveVarIds;
			colsArrayActiveVarIds = (s.in.lhsMinusRhsAccumulatedConstraint == null) ? new int[0] : aux_colsArrayActiveVarIds;
			if ((rowsArrayActiveVarIds.length != numActiveCoordinates) || (colsArrayActiveVarIds.length != numActiveCoordinates))
				throw new JOMException("JOM - IPOPT interface. Unexpected error");
			IntMatrix2D coord = IntFactory2D.dense.make(nele_jac, 2);
			if (s.in.lhsMinusRhsAccumulatedConstraint != null)
			{
				coord.viewColumn(0).assign(rowsArrayActiveVarIds);
				coord.viewColumn(1).assign(colsArrayActiveVarIds);
			}
			indexesOfGradientActiveCoordinates = (s.in.lhsMinusRhsAccumulatedConstraint == null) ?
					new int[0] :
					DoubleMatrixND.sub2ind(coord, IntFactory1D.dense.make(new int[]{s.in.numConstraints, s.in.numDecVariables})).toArray();

			create(n, x_L, x_U, m, g_L, g_U, nele_jac, nele_hess, index_style);
		}
		
		@Override
		protected boolean eval_f(int n, double[] x, boolean new_x, double[] obj_value) 
		{
			assert n == this.n;
			obj_value[0] = s.in.objectiveFunction.evaluate_internal(x).get(0);
			return true;
		}

		@Override
		protected boolean eval_grad_f(int n, double[] x, boolean new_x, double[] grad_f) 
		{
			assert n == this.n;
			System.arraycopy(s.in.objectiveFunction.evaluateJacobian_internal(x).viewRow(0).toArray() , 0 , grad_f , 0 , n);
			return true;
		}
		
		@Override
		protected boolean eval_g(int n, double[] x, boolean new_x, int m, double[] g) 
		{
			assert n == this.n;
			assert m == this.m;
			System.arraycopy(s.in.lhsMinusRhsAccumulatedConstraint.evaluate_internal(x).elements().toArray() , 0 , g , 0 , m);
			return true;
		}

		@Override
		protected boolean eval_jac_g(int n, double[] x, boolean new_x, int m, int nele_jac, int[] iRow, int[] jCol, double[] values) 
		{
			assert n == this.n;
			assert m == this.m;

			if (values == null) 
			{
				/* return the structure of the jacobian */
				/* Fill in the sparsity pattern --> fill in iRow and jCol arrays */
				System.arraycopy(rowsArrayActiveVarIds, 0, iRow , 0 , nele_jac);
				System.arraycopy(colsArrayActiveVarIds, 0, jCol , 0 , nele_jac);
			}
			else 
			{
				/* return the values of the jacobian of the constraints */
				DoubleMatrix2D res = s.in.lhsMinusRhsAccumulatedConstraint.evaluateJacobian_internal(x);
				double[] resSelectedIndexes = res.vectorize().viewSelection(indexesOfGradientActiveCoordinates).toArray();
				System.arraycopy(resSelectedIndexes , 0 , values , 0 , nele_jac);
			}
			return true;
		}
		@Override
		protected boolean eval_h(int n, double[] x, boolean new_x, double obj_factor, int m, double[] lambda,
				boolean new_lambda, int nele_hess, int[] iRow, int[] jCol, double[] values) 
		{
			if (!s.isLinearProblem()) // hessians are computed automatically
			return false;

			if (values == null) // problem is linear, hessian is zero
			{
				int[] auxVar = new int[nele_hess];
				java.util.Arrays.fill(auxVar, 0);
				/* Fill in the sparsity pattern --> fill in iRow and jCol arrays */
				System.arraycopy(auxVar , 0 , iRow , 0 , nele_hess);
				System.arraycopy(auxVar , 0 , jCol , 0 , nele_hess);
			} else
			{
				/* Return the values */
				double[] auxVar = new double[nele_hess];
				java.util.Arrays.fill(auxVar, 0);
				System.arraycopy(auxVar , 0 , values , 0 , nele_hess);
			}
			return true;
		}
		
	}
	
}
