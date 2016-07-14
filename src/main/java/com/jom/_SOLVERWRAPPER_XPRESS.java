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

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import com.dashoptimization.XPRS;
import com.dashoptimization.XPRSconstants;
import com.dashoptimization.XPRSenumerations;
import com.dashoptimization.XPRSprob;
import com.sun.jna.Native;
import com.sun.jna.Pointer;

import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.jet.math.tdouble.DoubleFunctions;

/** @author Pablo */
class _SOLVERWRAPPER_XPRESS
{
	private final HashMap<String, Object> param;
	private final _INTERNAL_SolverIO      s;
	private final String                  solverLibraryName;

	_SOLVERWRAPPER_XPRESS(_INTERNAL_SolverIO s, HashMap<String, Object> param)
	{
		this.s = s;
		this.solverLibraryName = (String) param.get("solverLibraryName");
		this.param = param;
	}

	int solve()
	{
		XPRSprob p = null;

		try
		{
			XPRS.init("C:\\xpressmp\\xpauth.xpr"); // PABLO: use hare solverLibraryName
			p = new XPRSprob();
            final int ncols = s.in.numDecVariables;
            final int nrows = s.in.numConstraints;
            final byte[] _srowtypes = new byte [s.in.numConstraints]; Arrays.fill(_srowtypes , (byte) 'R');
            final double[] _drhs = s.in.constraintUpperBound.toArray();
            final double [] _drange = s.in.constraintUpperBound.copy ().assign(s.in.constraintLowerBound,DoubleFunctions.minus).toArray();
            final double[] _dobj = s.in.objectiveFunction.getAffineExpression().getCellLinearCoefsFull(0);
			_INTERNAL_AffineExpressionCoefs constraintsMatrix = s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression();
			int[][] return_mnel = new int [0][]; 
			int[][] return_mrwind = new int [0][];
			int[][] return_mstart = new int [0][];
			double[][] return_dmatval = new double [0][];
			constraintsMatrix.getNonZerosRowColValForXpressSolver(return_mnel, return_mrwind, return_mstart, return_dmatval);
			final double [] _dlb = new double [ncols];
			final double [] _dub = new double [ncols];
			for (int col = 0 ; col < ncols ; col ++)
			{
				_dlb [col] = s.in.primalSolutionLowerBound.get(col) == -Double.MAX_VALUE? XPRSconstants.MINUSINFINITY : s.in.primalSolutionLowerBound.get(col);
				_dub [col] = s.in.primalSolutionUpperBound.get(col) == Double.MAX_VALUE? XPRSconstants.PLUSINFINITY : s.in.primalSolutionLowerBound.get(col);
			}
			final int ngents = s.in.primalSolutionIsInteger.zSum();
			final int nsets = 0;
			final byte[] _qgtype = new byte [ngents]; 
			int[] _mgcols = new int [ngents]; 
			if (ngents > 0)
			{
				Arrays.fill(_qgtype , (byte) 'I');
				IntArrayList nonZerosCols = new IntArrayList ();
				s.in.primalSolutionIsInteger.getNonZeros(nonZerosCols , new IntArrayList ());
				_mgcols = nonZerosCols.elements();
			}
			final double[] _dlim = null;
			final byte[] _stype = null;
			final int[] _msstart = null;
			final int[] _mscols = null;
			final double[] _dref = null;
			
			p.loadGlobal("" , ncols, nrows, _srowtypes, _drhs, _drange, _dobj, return_mstart [0],
					return_mnel [0], return_mrwind [0], return_dmatval [0], _dlb, _dub, ngents, nsets, _qgtype, _mgcols, _dlim,
					_stype, _msstart, _mscols, _dref);
			
			 p.chgObjSense(s.in.toMinimize? XPRSenumerations.ObjSense.MINIMIZE : XPRSenumerations.ObjSense.MAXIMIZE);

			 p.mipOptimize();

			 s.out.bestOptimalityBound = p.getDblAttrib(XPRSconstants.BESTBOUND);
			 s.out.statusCode = p.getIntAttrib(XPRSconstants.ERRORCODE);
			 s.out.statusMessage = p.getLastError();
			 
			 // pag 454: LPSTATUS --> infreasible, optimal, unbounded... 
			 // pag 457: MIPSTATUS --> infreasible, optimal, unbounded... 
			 // pag 467: STOPSTATUS --> stopped by time limit...
			 
			 // pag. 369: BARGAPSTOP, BARPRIMALSTOP, FEASTOL, FEASTOLTARGET: stop conditions
			 // PAG. 410: MAXTIME: maximum solver time
			 // MIPABSSTOP, MIPRELSTOP: stop criterion based on absolute gap
			 
			 s.out.primalCost = p.getDblAttrib(XPRSconstants.MIPBESTOBJVAL);
			 s.out.

			 
			 double [] primalSolution = new double [2];
			 double [] slackSolution = new double [2];
			 p.getMipSol(primalSolution , slackSolution);

			s.out.statusCode = solstat;
				s.out.statusMessage = getMsg_CPXgetstat(solstat);
				s.out.solutionIsOptimal = ((solstat == CPLEXConstants.CPXMIP_OPTIMAL) || (solstat == CPLEXConstants.CPXMIP_OPTIMAL_TOL));
				s.out.solutionIsFeasible = s.out.solutionIsOptimal || (solstat == CPLEXConstants.CPXMIP_NODE_LIM_FEAS) || (solstat == CPLEXConstants
						.CPXMIP_TIME_LIM_FEAS) || (solstat == CPLEXConstants.CPXMIP_FAIL_FEAS) || (solstat == CPLEXConstants.CPXMIP_MEM_LIM_FEAS) ||
						(solstat == CPLEXConstants.CPXMIP_ABORT_FEAS) || (solstat == CPLEXConstants.CPXMIP_FAIL_FEAS_NO_TREE);
				s.out.feasibleSolutionDoesNotExist = (solstat == CPLEXConstants.CPXMIP_INFEASIBLE);
				s.out.foundUnboundedSolution = (solstat == CPLEXConstants.CPXMIP_UNBOUNDED);

				if (!s.out.solutionIsFeasible)
				{
					s.problemAlreadyAttemptedTobeSolved = true;
					return s.out.statusCode;
				} // throw new JOMException("JOM - CPLEX interface. The problem could not be solved. CPLEX error message: " + s.out.statusMessage);


				/* Retrieve the optimal cost */
				double[] objval = new double[1];
				res = g.CPXgetobjval(env, lp, objval);
				if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetobjval");
				s.out.primalCost = objval[0];

				/* Check the number of constraitns and variables */
				int cur_numrows = g.CPXgetnumrows(env, lp);
				int cur_numcols = g.CPXgetnumcols(env, lp);
				if (cur_numrows != s.in.numConstraints) throw new JOMException("JOM - CPLEX interface. Unexpected error");
				if (cur_numcols != s.in.numDecVariables) throw new JOMException("JOM - CPLEX interface. Unexpected error");

				/* Retrieve the optimal primal solution */
				double[] x = new double[s.in.numDecVariables];
				res = g.CPXgetx(env, lp, x, 0, cur_numcols - 1);
				if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetx");
				s.out.primalSolution = DoubleFactory1D.dense.make(x);

				/* Retrieve the values of the constraints in the solution */
				double[] slack = new double[s.in.numConstraints];
				res = g.CPXgetslack(env, lp, slack, 0, cur_numrows - 1);
				if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetslack");
				s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(rhsCplex).assign(DoubleFactory1D.dense.make(slack), DoubleFunctions.minus);
				s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(s.in.numConstraints);
				s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
				s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);

				s.problemAlreadyAttemptedTobeSolved = true;

				/* Free memory */
				res = g.CPXfreeprob(env, lp);
				if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXfreeprob");

				return s.out.statusCode;
			 
			 

//			/* Set the environment parameters */
//			for (Entry<String, Object> entry : this.param.entrySet())
//			{
//				String keyStr = entry.getKey();
//				if (keyStr.equalsIgnoreCase("solverLibraryName")) continue;
//				int key = new Integer(entry.getKey());
//				Object value = entry.getValue();
//				if (value instanceof String)
//					aux = g.CPXsetstrparam(env, key, ((String) value).getBytes());
//				else if (value instanceof Integer)
//					aux = g.CPXsetintparam(env, key, ((Integer) value).intValue());
//				else if (value instanceof Double)
//					aux = g.CPXsetdblparam(env, key, ((Double) value).doubleValue());
//				else
//					throw new JOMException("JOM - CPLEX interface. Unknown value type in parameters");
//				if (aux != 0) throw new JOMException("JOM - CPLEX interface. Failure to set parameter: key = " + key);
//			}

			
//			p.loadGlobal("", 2, 2, _srowtypes , new double [] {1,1.5} , new double [] {1,1.5},new double [] {1,1} , 
//					 new int [] {0,2} , new int [] {2,2} , new int [] {0,1,0,1}, new double [] {1,1,-1,1}, new double [] {0,0} , 
//					 new double [] {XPRSconstants.PLUSINFINITY , XPRSconstants.PLUSINFINITY } , 
//					 2 , 0 , new byte [] {'I' , 'I'} , 
//					 new int [] {0,1} , new double [2] , new byte [2] , null ,  null , null);

			
//			PABLOOOO
//			- METER SOLUCION INICIAL
//			- METER ATRIBUTOS ETC.
//			- METER MAXSOLVERTIME


	
		} catch (RuntimeException e)
		{
			XPRS.free(); // frees the license
			if (p != null) p.destroy(); // frees any memory associated to the object
			e.printStackTrace(); throw e;
		} 
		

		/* Set maximization or minimization problem */
		g.CPXchgobjsen(env, lp, (s.in.toMinimize) ? _JNA_CPLEX.CPX_MIN : _JNA_CPLEX.CPX_MAX); /* Problem is maximization */

		/* Set the constraint bounds */
		byte[] sense = new byte[s.in.numConstraints];
		for (int contC = 0; contC < sense.length; contC++)
		{
			double lb = s.in.constraintLowerBound.get(contC); // 0 or -Double.MAX_VALUE
			double ub = s.in.constraintUpperBound.get(contC); // 0 or Double.MAX_VALUE
			if (lb == ub)
				sense[contC] = 'E';
			else if ((lb == -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
				throw new JOMException("JOM - CPLEX interface. Unbounded contraints are not allowed"); // cType.set(contC,GLP_FR);
			else if ((lb != -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
				sense[contC] = 'G';
			else if ((lb == -Double.MAX_VALUE) && (ub != Double.MAX_VALUE))
				sense[contC] = 'L';
			else
				throw new JOMException("JOM - CPLEX interface. Double bounded contraints are supposed not to exist in JOM");
		}

		/* Sets the problem columns: objective coef, lower and upper bounds */
		double[] coefObjectives = s.in.objectiveFunction.getAffineExpression().getCellLinearCoefsFull(0);
		byte[] xctype = new byte[s.in.numDecVariables];
		for (int dv = 0; dv < s.in.numDecVariables; dv++)
			if (s.in.primalSolutionIsInteger.get(dv) == 0)
				xctype[dv] = 'C';
			else
				xctype[dv] = ((s.in.primalSolutionLowerBound.get(dv) == 0) && (s.in.primalSolutionUpperBound.get(dv) == 1)) ? (byte) 'B' : (byte)
						'I';

		if (!s.in.hasIntegerVariables) xctype = null; // all variables are continuous and we can call CPXlpopt to solve the problem

		int res = g.CPXnewcols(env, lp, s.in.numDecVariables, coefObjectives, s.in.primalSolutionLowerBound.toArray(), s.in.primalSolutionUpperBound
				.toArray(), xctype, null);
		if (res != 0) throw new JOMException("JOM - CPLEX interface. Failure to set the problem columns (error code: " + res + ")");


		/* Sets the matrix coefficients */
		double[] rhsCplex = (s.in.numConstraints == 0) ?
				new double[0] :
				Arrays.copyOf(s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray(), s.in
						.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray().length);
		if (s.in.numConstraints > 0)
		{
			for (int cont = 0; cont < rhsCplex.length; cont++) rhsCplex[cont] = -rhsCplex[cont];
			res = g.CPXnewrows(env, lp, s.in.numConstraints, rhsCplex, sense, null, null);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failure to set the constraint bounds (error code: " + res + ")");

			_INTERNAL_AffineExpressionCoefs constraintsMatrix = s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression();
			int[][] rows = new int[1][];
			int[][] cols = new int[1][];
			double[][] vals = new double[1][];
			constraintsMatrix.getNonZerosRowColVal(rows, cols, vals);
			int[] ia = rows[0];
			int[] ja = cols[0];
			double[] ar = vals[0];

			res = g.CPXchgcoeflist(env, lp, ia.length, ia, ja, ar);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failure to set the constraint coefficients (error code: " + res + ")");
		}

		/* Execute the solver */
		if (s.in.hasIntegerVariables)
		{
			// integer problem
			res = g.CPXmipopt(env, lp);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed to optimize MIP (error code: " + res + ")");

			int solstat = g.CPXgetstat(env, lp);
			if (solstat == 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetstat");

			s.out.statusCode = solstat;
			s.out.statusMessage = getMsg_CPXgetstat(solstat);
			s.out.solutionIsOptimal = ((solstat == CPLEXConstants.CPXMIP_OPTIMAL) || (solstat == CPLEXConstants.CPXMIP_OPTIMAL_TOL));
			s.out.solutionIsFeasible = s.out.solutionIsOptimal || (solstat == CPLEXConstants.CPXMIP_NODE_LIM_FEAS) || (solstat == CPLEXConstants
					.CPXMIP_TIME_LIM_FEAS) || (solstat == CPLEXConstants.CPXMIP_FAIL_FEAS) || (solstat == CPLEXConstants.CPXMIP_MEM_LIM_FEAS) ||
					(solstat == CPLEXConstants.CPXMIP_ABORT_FEAS) || (solstat == CPLEXConstants.CPXMIP_FAIL_FEAS_NO_TREE);
			s.out.feasibleSolutionDoesNotExist = (solstat == CPLEXConstants.CPXMIP_INFEASIBLE);
			s.out.foundUnboundedSolution = (solstat == CPLEXConstants.CPXMIP_UNBOUNDED);

			if (!s.out.solutionIsFeasible)
			{
				s.problemAlreadyAttemptedTobeSolved = true;
				return s.out.statusCode;
			} // throw new JOMException("JOM - CPLEX interface. The problem could not be solved. CPLEX error message: " + s.out.statusMessage);


			/* Retrieve the optimal cost */
			double[] objval = new double[1];
			res = g.CPXgetobjval(env, lp, objval);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetobjval");
			s.out.primalCost = objval[0];

			/* Check the number of constraitns and variables */
			int cur_numrows = g.CPXgetnumrows(env, lp);
			int cur_numcols = g.CPXgetnumcols(env, lp);
			if (cur_numrows != s.in.numConstraints) throw new JOMException("JOM - CPLEX interface. Unexpected error");
			if (cur_numcols != s.in.numDecVariables) throw new JOMException("JOM - CPLEX interface. Unexpected error");

			/* Retrieve the optimal primal solution */
			double[] x = new double[s.in.numDecVariables];
			res = g.CPXgetx(env, lp, x, 0, cur_numcols - 1);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetx");
			s.out.primalSolution = DoubleFactory1D.dense.make(x);

			/* Retrieve the values of the constraints in the solution */
			double[] slack = new double[s.in.numConstraints];
			res = g.CPXgetslack(env, lp, slack, 0, cur_numrows - 1);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetslack");
			s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(rhsCplex).assign(DoubleFactory1D.dense.make(slack), DoubleFunctions.minus);
			s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(s.in.numConstraints);
			s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
			s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);

			s.problemAlreadyAttemptedTobeSolved = true;

			/* Free memory */
			res = g.CPXfreeprob(env, lp);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXfreeprob");

			return s.out.statusCode;
		} else
		{
			// continuous problem
			res = g.CPXlpopt(env, lp);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed to optimize LP (error code: " + res + ")");

			int solstat = g.CPXgetstat(env, lp);
			if (solstat == 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetstat");

			s.out.statusCode = solstat;
			s.out.statusMessage = getMsg_CPXgetstat(solstat);
			s.out.solutionIsOptimal = (solstat == CPLEXConstants.CPX_STAT_OPTIMAL);
			s.out.solutionIsFeasible = s.out.solutionIsOptimal || (solstat == CPLEXConstants.CPX_STAT_NUM_BEST);
			s.out.feasibleSolutionDoesNotExist = (solstat == CPLEXConstants.CPX_STAT_INFEASIBLE);
			s.out.foundUnboundedSolution = (solstat == CPLEXConstants.CPX_STAT_UNBOUNDED);

			if (!s.out.solutionIsFeasible)
				throw new JOMException("JOM - CPLEX interface. CPLEX: The problem could not be solved. CPLEX error message: " + s.out.statusMessage);

			/* Retrieve the optimal cost */
			double[] objval = new double[1];
			res = g.CPXgetobjval(env, lp, objval);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetobjval");
			s.out.primalCost = objval[0];

			/* Check the number of constraitns and variables */
			int cur_numrows = g.CPXgetnumrows(env, lp);
			int cur_numcols = g.CPXgetnumcols(env, lp);
			if (cur_numrows != s.in.numConstraints) throw new JOMException("JOM - CPLEX interface. Unexpected error");
			if (cur_numcols != s.in.numDecVariables) throw new JOMException("JOM - CPLEX interface. Unexpected error");

			/* Retrieve the optimal primal solution */
			double[] x = new double[s.in.numDecVariables];
			res = g.CPXgetx(env, lp, x, 0, cur_numcols - 1);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetx");
			s.out.primalSolution = DoubleFactory1D.dense.make(x);

			/* Retrieve the values of the constraints in the solution */
			double[] slack = new double[s.in.numConstraints];
			res = g.CPXgetslack(env, lp, slack, 0, cur_numrows - 1);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetslack");
			s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(rhsCplex).assign(DoubleFactory1D.dense.make(slack), DoubleFunctions.minus);

			/* Retrieve multipliers of constraints */
			double[] pi = new double[s.in.numConstraints];
			res = g.CPXgetpi(env, lp, pi, 0, pi.length - 1);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetpi");
			s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(pi);

			/* Retrieve multipliers of variable bounds */
			double[] dj = new double[s.in.numDecVariables];
			res = g.CPXgetdj(env, lp, dj, 0, dj.length - 1);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXgetdj");
			s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables, 0.0);
			s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables, 0.0);
			for (int dv = 0; dv < s.in.numDecVariables; dv++)
			{
				double lb = s.in.primalSolutionLowerBound.get(dv);
				double ub = s.in.primalSolutionUpperBound.get(dv);
				double val = s.out.primalSolution.get(dv);
				if (lb == ub)
				{
					s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(dv, dj[dv]);
					s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(dv, -dj[dv]);
				} else if (Math.abs(val - lb) > Math.abs(val - ub))
					s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(dv, dj[dv]);
				else
					s.out.multiplierOfLowerBoundConstraintToPrimalVariables.set(dv, dj[dv]);
			}

			s.problemAlreadyAttemptedTobeSolved = true;

			/* Free memory */
			res = g.CPXfreeprob(env, lp);
			if (res != 0) throw new JOMException("JOM - CPLEX interface. Failed in call to CPXfreeprob");

			return s.out.statusCode;
		}

	}

	private final String getMsg_CPXgetstat(int s)
	{
		switch (s)
		{
			case CPLEXConstants.CPXMIP_OPTIMAL:
				return "Optimal integer solution has been found";
			case CPLEXConstants.CPXMIP_OPTIMAL_TOL:
				return "Optimal solution with the tolerance defined by epgap or epagap has been found";
			case CPLEXConstants.CPXMIP_INFEASIBLE:
				return "Solution is integer infeasible ";
			case CPLEXConstants.CPXMIP_SOL_LIM:
				return "The limit on mixed integer solutions has been reached ";
			case CPLEXConstants.CPXMIP_NODE_LIM_FEAS:
				return "Node limit has been exceeded but integer solution exists ";
			case CPLEXConstants.CPXMIP_NODE_LIM_INFEAS:
				return "Node limit has been reached; no integer solution ";
			case CPLEXConstants.CPXMIP_TIME_LIM_FEAS:
				return "Time limit exceeded, but integer solution exists";
			case CPLEXConstants.CPXMIP_TIME_LIM_INFEAS:
				return "Time limit exceeded; no integer solution ";
			case CPLEXConstants.CPXMIP_FAIL_FEAS:
				return "Terminated because of an error, but integer solution exists ";
			case CPLEXConstants.CPXMIP_FAIL_INFEAS:
				return "Terminated because of an error; no integer solution ";
			case CPLEXConstants.CPXMIP_MEM_LIM_FEAS:
				return "Limit on tree memory has been reached, but an integer solution exists ";
			case CPLEXConstants.CPXMIP_MEM_LIM_INFEAS:
				return "Limit on tree memory has been reached; no integer solution ";
			case CPLEXConstants.CPXMIP_ABORT_FEAS:
				return "Stopped, but an integer solution exists ";
			case CPLEXConstants.CPXMIP_ABORT_INFEAS:
				return "Stopped; no integer solution ";
			case CPLEXConstants.CPXMIP_OPTIMAL_INFEAS:
				return "Problem is optimal with unscaled infeasibilities ";
			case CPLEXConstants.CPXMIP_FAIL_FEAS_NO_TREE:
				return "Out of memory, no tree available, integer solution exists ";
			case CPLEXConstants.CPXMIP_FAIL_INFEAS_NO_TREE:
				return "Out of memory, no tree available, no integer solution ";
			case CPLEXConstants.CPXMIP_UNBOUNDED:
				return "Problem has an unbounded ray ";
			case CPLEXConstants.CPXMIP_INForUNBD:
				return "Problem has been proved either infeasible or unbounded";
			case CPLEXConstants.CPXMIP_FEASIBLE_RELAXED_SUM:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_MIN_SUM (or MinSum) on a mixed integer "
						+ "problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that "
						+ "relaxation) was installed. The relaxation is minimal.";
			case CPLEXConstants.CPXMIP_OPTIMAL_RELAXED_SUM:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_OPT_SUM (or OptSum) on a mixed integer "
						+ "problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that "
						+ "relaxation) was installed. The relaxation is optimal.";
			case CPLEXConstants.CPXMIP_FEASIBLE_RELAXED_INF:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_MIN_INF (or MinInf) on a mixed integer "
						+ "problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that "
						+ "relaxation) was installed. The relaxation is minimal.";
			case CPLEXConstants.CPXMIP_OPTIMAL_RELAXED_INF:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_OPT_INF (or OptInf) on a mixed integer "
						+ "problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that "
						+ "relaxation) was installed. The relaxation is optimal.";
			case CPLEXConstants.CPXMIP_FEASIBLE_RELAXED_QUAD:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_MIN_QUAD (or MinQuad) on a mixed integer "
						+ "problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that "
						+ "relaxation) was installed. The relaxation is minimal. ";
			case CPLEXConstants.CPXMIP_OPTIMAL_RELAXED_QUAD:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_OPT_QUAD (or OptQuad) on a mixed integer "
						+ "problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that "
						+ "relaxation) was installed. The relaxation is optimal.";
			case CPLEXConstants.CPXMIP_ABORT_RELAXED:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt), "
						+ "when the algorithm terminates prematurely, for example after reaching a limit. This status means that a relaxed solution "
						+ "is available and can be queried ";
			case CPLEXConstants.CPXMIP_FEASIBLE:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "on a MIP problem. The problem under consideration was found to be feasible after phase 1 of FeasOpt. A feasible solution "
						+ "is available. This status is also used in the status field of solution and mipstart files for solutions from the solution"
						+ " pool.";
			case CPLEXConstants.CPXMIP_POPULATESOL_LIM:
				return "This status occurs only after a call to the Callable Library routine CPXpopulate (or the Concert Technology method populate)"
						+ " on a MIP problem. The limit on mixed integer solutions generated by populate, as specified by the parameter "
						+ "CPX_PARAM_POPULATELIM, has been reached.";
			case CPLEXConstants.CPXMIP_OPTIMAL_POPULATED:
				return "This status occurs only after a call to the Callable Library routine CPXpopulate (or the Concert Technology method populate)"
						+ " on a MIP problem. Populate has completed the enumeration of all solutions it could enumerate.";
			case CPLEXConstants.CPXMIP_OPTIMAL_POPULATED_TOL:
				return "This status occurs only after a call to the Callable Library routine CPXpopulate (or the Concert Technology method populate)"
						+ " on a MIP problem. Populate has completed the enumeration of all solutions it could enumerate whose objective value fit "
						+ "the tolerance specified by the parameters CPX_PARAM_SOLNPOOLAGAP and CPX_PARAM_SOLNPOOLGAP";

			// codes for LP programs (non integer)
			case CPLEXConstants.CPX_STAT_OPTIMAL:
				return "Optimal solution is available";
			case CPLEXConstants.CPX_STAT_UNBOUNDED:
				return "Problem has an unbounded ray; see the concept Unboundedness for more information about infeasibility and unboundedness as a "
						+ "solution status";
			case CPLEXConstants.CPX_STAT_INFEASIBLE:
				return "Problem has been proven infeasible; see the topic Interpreting Solution Quality in the CPLEX User's Manual for more details";
			case CPLEXConstants.CPX_STAT_INForUNBD:
				return "Problem has been proven either infeasible or unbounded; see the topic Effect of Preprocessing on Feasibility in the CPLEX "
						+ "User's Manual for more detail. ";
			case CPLEXConstants.CPX_STAT_OPTIMAL_INFEAS:
				return "Optimal solution is available, but with infeasibilities after unscaling";
			case CPLEXConstants.CPX_STAT_NUM_BEST:
				return "Solution is available, but not proved optimal, due to numeric difficulties during optimization";
			case CPLEXConstants.CPX_STAT_ABORT_IT_LIM:
				return "Stopped due to limit on number of iterations.";
			case CPLEXConstants.CPX_STAT_ABORT_TIME_LIM:
				return "Stopped due to a time limit.";
			case CPLEXConstants.CPX_STAT_ABORT_OBJ_LIM:
				return "Stopped due to an objective limit.";
			case CPLEXConstants.CPX_STAT_ABORT_USER:
				return "Stopped due to a request from the user.";
			case CPLEXConstants.CPX_STAT_FEASIBLE_RELAXED_SUM:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_MIN_SUM (or MinSum) on a continuous problem."
						+ " A relaxation was successfully found and a feasible solution for the problem. (if relaxed according to that relaxation) "
						+ "was installed. The relaxation is minimal.";
			case CPLEXConstants.CPX_STAT_OPTIMAL_RELAXED_SUM:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_OPT_SUM (or OptSum) on a continuous problem."
						+ " A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that relaxation) "
						+ "was installed. The relaxation is optimal";
			case CPLEXConstants.CPX_STAT_FEASIBLE_RELAXED_INF:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_MIN_INF (or MinInf) on a continuous problem."
						+ " A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that relaxation) "
						+ "was installed. The relaxation is minimal.";
			case CPLEXConstants.CPX_STAT_OPTIMAL_RELAXED_INF:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_OPT_INF (or OptInf) on a continuous problem."
						+ " A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that relaxation) "
						+ "was installed. The relaxation is optimal.";
			case CPLEXConstants.CPX_STAT_FEASIBLE_RELAXED_QUAD:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_MIN_QUAD (or MinQuad) on a continuous "
						+ "problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that "
						+ "relaxation) was installed. The relaxation is minimal.";
			case CPLEXConstants.CPX_STAT_OPTIMAL_RELAXED_QUAD:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) "
						+ "with the parameter CPX_PARAM_FEASOPTMODE (or FeasOptMode) set to CPX_FEASOPT_OPT_QUAD (or OptQuad) on a continuous problem. A relaxation was successfully found and a feasible solution for the problem (if relaxed according to that relaxation) was installed. The relaxation is optimal.";
			case CPLEXConstants.CPX_STAT_FEASIBLE:
				return "This status occurs only after a call to the Callable Library routine CPXfeasopt (or the Concert Technology method feasOpt) on a continuous problem. The problem under consideration was found to be feasible after phase 1 of FeasOpt. A feasible solution is available. ";
		}
		return "Unknown code";
	}

	class CPLEXConstants
	{
		static public final int CPX_BASIC_SOLN                 = 1;
		/* Solution type return values from CPXsolninfo */
		static public final int CPX_NO_SOLN                    = 0;
		static public final int CPX_NONBASIC_SOLN              = 2;
		static public final int CPX_PRIMAL_SOLN                = 3;
		static public final int CPX_STAT_ABORT_IT_LIM          = 10;
		static public final int CPX_STAT_ABORT_OBJ_LIM         = 12;
		static public final int CPX_STAT_ABORT_TIME_LIM        = 11;
		static public final int CPX_STAT_ABORT_USER            = 13;
		static public final int CPX_STAT_FEASIBLE              = 23;
		static public final int CPX_STAT_FEASIBLE_RELAXED_INF  = 16;
		static public final int CPX_STAT_FEASIBLE_RELAXED_QUAD = 18;
		static public final int CPX_STAT_FEASIBLE_RELAXED_SUM  = 14;
		static public final int CPX_STAT_INFEASIBLE            = 3;
		static public final int CPX_STAT_INForUNBD             = 4;
		static public final int CPX_STAT_NUM_BEST              = 6;
		/* Values returned for 'stat' by solution */
		static public final int CPX_STAT_OPTIMAL               = 1;
		static public final int CPX_STAT_OPTIMAL_INFEAS        = 5;
		static public final int CPX_STAT_OPTIMAL_RELAXED_INF   = 17;
		static public final int CPX_STAT_OPTIMAL_RELAXED_QUAD  = 19;
		static public final int CPX_STAT_OPTIMAL_RELAXED_SUM   = 15;
		static public final int CPX_STAT_UNBOUNDED             = 2;
		static public final int CPXMIP_ABORT_FEAS              = 113;
		static public final int CPXMIP_ABORT_INFEAS            = 114;
		static public final int CPXMIP_ABORT_RELAXED           = 126;
		static public final int CPXMIP_FAIL_FEAS               = 109;
		static public final int CPXMIP_FAIL_FEAS_NO_TREE       = 116;
		static public final int CPXMIP_FAIL_INFEAS             = 110;
		static public final int CPXMIP_FAIL_INFEAS_NO_TREE     = 117;
		static public final int CPXMIP_FEASIBLE                = 127;
		static public final int CPXMIP_FEASIBLE_RELAXED_INF    = 122;

		static public final int CPXMIP_FEASIBLE_RELAXED_QUAD = 124;
		static public final int CPXMIP_FEASIBLE_RELAXED_SUM  = 120;
		static public final int CPXMIP_INFEASIBLE            = 103;
		static public final int CPXMIP_INForUNBD             = 119;

		static public final int CPXMIP_MEM_LIM_FEAS          = 111;
		static public final int CPXMIP_MEM_LIM_INFEAS        = 112;
		static public final int CPXMIP_NODE_LIM_FEAS         = 105;
		static public final int CPXMIP_NODE_LIM_INFEAS       = 106;
		/* MIP status codes */
		static public final int CPXMIP_OPTIMAL               = 101;
		static public final int CPXMIP_OPTIMAL_INFEAS        = 115;
		static public final int CPXMIP_OPTIMAL_POPULATED     = 129;
		static public final int CPXMIP_OPTIMAL_POPULATED_TOL = 130;
		static public final int CPXMIP_OPTIMAL_RELAXED_INF   = 123;
		static public final int CPXMIP_OPTIMAL_RELAXED_QUAD  = 125;
		static public final int CPXMIP_OPTIMAL_RELAXED_SUM   = 121;
		static public final int CPXMIP_OPTIMAL_TOL           = 102;
		static public final int CPXMIP_POPULATESOL_LIM       = 128;
		static public final int CPXMIP_SOL_LIM               = 104;
		static public final int CPXMIP_TIME_LIM_FEAS         = 107;
		static public final int CPXMIP_TIME_LIM_INFEAS       = 108;
		static public final int CPXMIP_UNBOUNDED             = 118;
	}

}
