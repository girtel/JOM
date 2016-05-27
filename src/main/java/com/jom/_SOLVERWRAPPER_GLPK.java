/*******************************************************************************
 * Copyright (c) 2015 Pablo Pavon Mariño.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * 
 * Contributors:
 *     Pablo Pavon Mariño - initial API and implementation
 ******************************************************************************/



 




/**
 * 
 */
package com.jom;

import java.util.HashMap;
import java.util.Map.Entry;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;

import com.sun.jna.Native;
import com.sun.jna.Pointer;

class _SOLVERWRAPPER_GLPK
{
	private final HashMap<String, Object> param;
	private final _INTERNAL_SolverIO s;
	private final String solverLibraryName;
	
	
	private static int GLP_BV = 3;


	private static int GLP_CV = 1;

	private static int GLP_DB = 4;

	private static int GLP_EBADB = 0x01;

	private static int GLP_EBOUND = 0x04;

	
	private static int GLP_ECOND = 0x03; 
	private static int GLP_EFAIL = 0x05;
	private static int GLP_EINSTAB = 0x11;
	private static int GLP_EITLIM = 0x08;
	private static int GLP_EMIPGAP = 0x0e;
	private static int GLP_ENOCVG = 0x10;
	private static int GLP_ENODFS = 0x0B;
	private static int GLP_ENOPFS = 0x0a;

	private static int GLP_EOBJLL  = 0x06; 
	private static int GLP_EOBJUL = 0x07; 
	private static int GLP_EROOT = 0x0c;
	private static int GLP_ESING = 0x02;
	private static int GLP_ESTOP = 0x0d;
	private static int GLP_ETMLIM = 0x09;
	private static int GLP_FEAS = 2; 
	private static int GLP_FR = 1; 
	
	private static int GLP_FX = 5;
	private static int GLP_INFEAS = 3;
	private static int GLP_IV = 2;
	private static int GLP_LO = 2;
	private static int GLP_MAX = 2;
	private static int GLP_MIN = 1;

	private static int GLP_NOFEAS = 4;
	private static int GLP_OPT = 5;

	private static int GLP_UNBND = 6;
	private static int GLP_UNDEF = 1;
	private static int GLP_UP = 3;

	_SOLVERWRAPPER_GLPK(_INTERNAL_SolverIO s, HashMap<String, Object> param)
	{
		this.s = s;
		this.solverLibraryName = (String) param.get("solverLibraryName");
		this.param = param;
	}
	String glpk_errorMessage_solve(_INTERNAL_SolverIO s, int errorCode)
	{
		if (s.in.hasIntegerVariables)
		{
			if (errorCode == 0) return "The MIP problem instance has been successfully solved. (This code does not necessarily mean that the solver has found optimal solution. It only means that the solution process was successful.)";
			if (errorCode == GLP_EBOUND) return "Unable to start the search, because some double-bounded variables have incorrect bounds or some integer variables have non-integer (fractional) bounds.";
			if (errorCode == GLP_EROOT) return "Unable to start the search, because optimal basis for initial LP relaxation is not provided. (This code may appear only if the presolver is disabled.)";
			if (errorCode == GLP_ENOPFS) return "Unable to start the search, because LP relaxation of the MIP problem instance has no primal feasible solution. (This code may appear only if the presolver is enabled.)";
			if (errorCode == GLP_ENODFS) return "Unable to start the search, because LP relaxation of the MIP problem instance has no dual feasible solution. In other word, this code means that if the LP relaxation has at least one primal feasible solution, its optimal solution is unbounded, so if the MIP problem has at least one integer feasible solution, its (integer) optimal solution is also unbounded. (This code may appear only if the presolver is enabled.)";
			if (errorCode == GLP_EFAIL) return "The search was prematurely terminated due to the solver failure.";
			if (errorCode == GLP_EMIPGAP) return "The search was prematurely terminated, because the relative mip gap tolerance has been reached.";
			if (errorCode == GLP_ETMLIM) return "The search was prematurely terminated, because the time limit has been exceeded.";
			if (errorCode == GLP_ESTOP) return "The search was prematurely terminated by application. (This code may appear only if the advanced solver interface is used.)";
		} else if (param.get("glpkSolverType").equals("simplex"))
		{
			if (errorCode == 0) return "The LP problem instance has been successfully solved. (This code does not necessarily mean that the solver has found optimal solution. It only means that the solution process was successful.)";
			if (errorCode == GLP_EBADB) return "Unable to start the search, because the initial basis specified in the problem object is invalid|the number of basic (auxiliary and structural) variables is not the same as the number of rows in the problem object.";
			if (errorCode == GLP_ESING) return "Unable to start the search, because the basis matrix corresponding to the initial basis is singular within the working precision. ";
			if (errorCode == GLP_ECOND) return "Unable to start the search, because the basis matrix corresponding to the initial basis is ill-conditioned, i.e. its condition number is too large.";
			if (errorCode == GLP_EBOUND) return "Unable to start the search, because some double-bounded (auxiliary or structural) variables have incorrect bounds. ";
			if (errorCode == GLP_EFAIL) return "The search was prematurely terminated due to the solver failure.";
			if (errorCode == GLP_EOBJLL) return "The search was prematurely terminated, because the objective function being maximized has reached its lower limit and continues decreasing (the dual simplex only).";
			if (errorCode == GLP_EOBJUL) return "The search was prematurely terminated, because the objective function being minimized has reached its upper limit and continues increasing (the dual simplex only).";
			if (errorCode == GLP_EITLIM) return "The search was prematurely terminated, because the simplex iteration limit has been exceeded.";
			if (errorCode == GLP_ETMLIM) return "The search was prematurely terminated, because the time limit has been exceeded.";
			if (errorCode == GLP_ENOPFS) return "The LP problem instance has no primal feasible solution (only if the LP presolver is used).";
			if (errorCode == GLP_ENODFS) return "The LP problem instance has no dual feasible solution (only if the LP presolver is used).";
		} else if (param.get("glpkSolverType").equals("interior-point"))
		{
			if (errorCode == 0) return "The LP problem instance has been successfully solved. (This code does not necessarily mean that the solver has found optimal solution. It only means that the solution process was successful.)";
			if (errorCode == GLP_EFAIL) return "The problem has no rows/columns.";
			if (errorCode == GLP_ENOCVG) return "Very slow convergence or divergence.";
			if (errorCode == GLP_EITLIM) return "Iteration limit exceeded.";
			if (errorCode == GLP_EINSTAB) return "Numerical instability on solving Newtonian system.";
		}
		return "Unkown message!!!";
	}
	String glpk_statusMessage_solve(_INTERNAL_SolverIO s, int status)
	{
		if (s.in.hasIntegerVariables)
		{
			if (status == GLP_UNDEF) return "MIP solution is undefined.";
			if (status == GLP_OPT) return "MIP solution is integer optimal.";
			if (status == GLP_FEAS) return "MIP solution is integer feasible, however, its optimality (or non-optimality) has not been proven, perhaps due to premature termination of the search.";
			if (status == GLP_NOFEAS) return "problem has no integer feasible solution (proven by the solver).";
		} else if (param.get("glpkSolverType").equals("simplex"))
		{
			if (status == GLP_OPT) return "Solution is optimal";
			if (status == GLP_FEAS) return "Solution is feasible";
			if (status == GLP_INFEAS) return "Solution is infeasible";
			if (status == GLP_NOFEAS) return "Problem has no feasible solution";
			if (status == GLP_UNBND) return "Problem has unbounded solution";
			if (status == GLP_UNDEF) return "Solution is undefined";
		} else if (param.get("glpkSolverType").equals("interior-point"))
		{
			if (status == GLP_UNDEF) return "interior-point solution is undeFIned.";
			if (status == GLP_OPT) return "interior-point solution is optimal";
			if (status == GLP_INFEAS) return "interior-point solution is infeasible";
			if (status == GLP_NOFEAS) return "no feasible primal-dual solution exists";
		}
		return "Unknown message!!!";
	}
	int solve()
	{
		_JNA_GLPK g = (_JNA_GLPK) Native.loadLibrary(this.solverLibraryName , _JNA_GLPK.class);
		Pointer lp;

		/* Compute number of decision variables */
		lp = g.glp_create_prob();

		/* Set objective direction: min / max */
		if (s.in.toMinimize)
			g.glp_set_obj_dir(lp, GLP_MIN);
		else
			g.glp_set_obj_dir(lp, GLP_MAX);
		double [] objectiveFunctionCoefs = s.in.objectiveFunction.getAffineExpression().getCellLinearCoefsFull(0);

		/* Initialize columns (decision variables) of the problem */
		/* At the same time, initialize objective function coefficients */
		g.glp_add_cols(lp, s.in.numDecVariables);
		IntMatrix1D dvType = IntFactory1D.dense.make(s.in.numDecVariables);
		for (int contDV = 0; contDV < s.in.numDecVariables; contDV++)
		{
			if (s.in.primalSolutionIsInteger.get(contDV) == 1)
				g.glp_set_col_kind(lp, 1 + contDV, GLP_IV);
			else
				g.glp_set_col_kind(lp, 1 + contDV, GLP_CV);
			double lb = s.in.primalSolutionLowerBound.get(contDV);
			double ub = s.in.primalSolutionUpperBound.get(contDV);
			if (lb == ub)
				dvType.set(contDV, GLP_FX);
			else if ((lb == -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
				dvType.set(contDV, GLP_FR);
			else if ((lb != -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
				dvType.set(contDV, GLP_LO);
			else if ((lb == -Double.MAX_VALUE) && (ub != Double.MAX_VALUE))
				dvType.set(contDV, GLP_UP);
			else
				dvType.set(contDV, GLP_DB);
			g.glp_set_col_bnds(lp, 1 + contDV, dvType.get(contDV), lb, ub);

			/* Set value of the objective function coefficients */
			g.glp_set_obj_coef(lp, 1 + contDV, objectiveFunctionCoefs [contDV]);
		}
		/* the constant value in the objective */
		g.glp_set_obj_coef(lp, 0, s.in.objectiveFunction.getAffineExpression().getCellConstantCoef(0));

		/* Create constraint bounds */
		if (s.in.lhsMinusRhsAccumulatedConstraint != null)
			g.glp_add_rows(lp, s.in.lhsMinusRhsAccumulatedConstraint.getNumScalarExpressions());
		IntMatrix1D cType = IntFactory1D.dense.make(s.in.numConstraints);
		for (int contC = 0; contC < s.in.numConstraints; contC++)
		{
			double lb = s.in.constraintLowerBound.get(contC); // 0 or -Double.MAX_VALUE
			double ub = s.in.constraintUpperBound.get(contC); // 0 or Double.MAX_VALUE
			if (lb == ub)
			{
				cType.set(contC, GLP_FX);
				lb = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(contC); // equal to minus the constant value
				ub = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(contC); // equal to minus the constant value
			} else if ((lb == -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
			{
				throw new JOMException("JOM - GLPK interface. Unbounded contraints are not allowed"); // cType.set(contC,GLP_FR);
			} else if ((lb != -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
			{
				cType.set(contC, GLP_LO);
				lb = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(contC); // equal to minus the constant value
			} else if ((lb == -Double.MAX_VALUE) && (ub != Double.MAX_VALUE))
			{
				cType.set(contC, GLP_UP);
				ub = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(contC); // equal to minus the constant value
			} else
			{
				throw new JOMException("JOM - GLPK interface. Double bounded constraints are supposed not to exist in JOM");
				// cType.set(contC, GLP_DB);
			}
			g.glp_set_row_bnds(lp, 1 + contC, cType.get(contC), lb, ub);
		}

		/* Load the constraints matrix */
		if (s.in.lhsMinusRhsAccumulatedConstraint != null)
		{
			_INTERNAL_AffineExpressionCoefs constraintsMatrix = s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression();//.linear_getLinearCoefs();

			int [][] rows = new int [1][]; int [][] cols = new int [1][]; double [][] vals = new double [1][];  
			constraintsMatrix.getNonZerosRowColVal (rows , cols , vals);
			
			int [] ia = rows [0]; int [] ja = cols [0]; double [] ar = vals [0];
			// the indexes must be 1-indexed and not 0-indexes. Also the arrays should
			// have one more cell than required, firt cell is not used
			int [] ia_added1 = new int [1 + ia.length]; // extra cell, since the index 0 is not used (GLPK things...)
			int [] ja_added1 = new int [1 + ia.length]; // extra cell, since the index 0 is not used (GLPK things...)
			double [] ar_extraCell = new double [1 + ia.length]; // extra cell, since the index 0 is not used (GLPK things...)
			System.arraycopy(ia, 0, ia_added1, 1, ia.length);
			System.arraycopy(ja, 0, ja_added1, 1, ja.length);
			System.arraycopy(ar, 0, ar_extraCell, 1, ar.length);
			
			for (int cont = 1 ; cont < ia_added1.length ; cont ++) { ia_added1 [cont] += 1; ja_added1 [cont] += 1; } // sum 1 since it is 1-indexed  
			g.glp_load_matrix(lp, ia.length , ia_added1, ja_added1, ar_extraCell);
		}
		s.out.primalSolution = DoubleFactory1D.dense.make(s.in.numDecVariables);
		s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
		s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
		s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(s.in.numConstraints);
		s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(s.in.numConstraints);

		/* If MILP use the MILP solver */
		if (s.in.hasIntegerVariables)
		{
			_JNA_GLPK.GLP_IOCP parm = new _JNA_GLPK.GLP_IOCP();
			g.glp_init_iocp(parm.getPointer());
			parm.writeField("presolve", _JNA_GLPK.GLP_ON);
			/* Set specific parameters specified by the user */
			for (Entry<String, Object> entry : this.param.entrySet())
			{
				String key = entry.getKey();
				if (key.equalsIgnoreCase("solverLibraryName")) continue; // a non-GPLK specific option
				if (key.equalsIgnoreCase("glpkSolverType")) continue; // a non-GPLK specific option
				Object value = entry.getValue();
				if (value instanceof Integer)
					parm.writeField(key , ((Integer) value).intValue());
				else if (value instanceof Double)
					parm.writeField(key , ((Double) value).doubleValue());
				else
					throw new JOMException("JOM - GLPK interface. Unknown value type in parameters");
			}

			int ret = g.glp_intopt(lp, parm.getPointer());
			// Retrieve solution
			if ((ret == 0) || (ret == GLP_ETMLIM) || (ret == GLP_EMIPGAP)) // optimum, time limit, or optimality gap limit
			{
				int status = g.glp_mip_status(lp);
				s.out.statusCode = status;
				s.out.statusMessage = glpk_statusMessage_solve(s, status);
				s.out.solutionIsFeasible = (status == GLP_OPT) || (status == GLP_FEAS);
				s.out.solutionIsOptimal = (status == GLP_OPT);
				s.out.feasibleSolutionDoesNotExist = (status == GLP_NOFEAS);
				s.out.foundUnboundedSolution = false;
				if (!s.out.solutionIsFeasible) throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_statusMessage_solve(s, status));

				s.out.primalCost = g.glp_mip_obj_val(lp);
				int nCols = g.glp_get_num_cols(lp);
				if (nCols != s.in.numDecVariables) throw new JOMException("JOM - GLPK interface. Unexptected error");
				for (int cont = 1; cont <= nCols; cont++)
					s.out.primalSolution.set(cont - 1, g.glp_mip_col_val(lp, cont));

				int nRows = g.glp_get_num_rows(lp);
				if (nRows != s.in.numConstraints) throw new JOMException("JOM - GLPK interface. Unexpected error");
				for (int cont = 1; cont <= nRows; cont++)
					s.out.primalValuePerConstraint.set(cont - 1, g.glp_mip_row_val(lp, cont));
			} else
				throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_errorMessage_solve(s, ret));
			s.problemAlreadyAttemptedTobeSolved = true;
			g.glp_delete_prob(lp);
			return ret;
		}

		if ("interior-point".equals(param.get("glpkSolverType")))
		{
			/* use interior point method */
			_JNA_GLPK.GLP_IPTCP parm = new _JNA_GLPK.GLP_IPTCP();
			g.glp_init_iptcp(parm.getPointer());
			/* Set specific parameters specified by the user */
			for (Entry<String, Object> entry : this.param.entrySet())
			{
				String key = entry.getKey();
				if (key.equalsIgnoreCase("solverLibraryName")) continue; // a non-GPLK specific option
				if (key.equalsIgnoreCase("glpkSolverType")) continue; // a non-GPLK specific option
				if (key.equalsIgnoreCase("presolve")) continue; // not valid for interior-point
				Object value = entry.getValue();
				if (value instanceof Integer)
					parm.writeField(key , ((Integer) value).intValue());
				else if (value instanceof Double)
					parm.writeField(key , ((Double) value).doubleValue());
				else
					throw new JOMException("JOM - GLPK interface. Unknown value type in parameters");
			}

			int ret = g.glp_interior(lp, parm.getPointer());
			if (ret == 0)
			{
				int status = g.glp_ipt_status(lp);
				s.out.statusCode = status;
				s.out.statusMessage = glpk_statusMessage_solve(s, status);
				s.out.solutionIsFeasible = (status == GLP_OPT);
				s.out.solutionIsOptimal = (status == GLP_OPT);
				s.out.feasibleSolutionDoesNotExist = (status == GLP_NOFEAS);
				s.out.foundUnboundedSolution = false;
				if (!s.out.solutionIsFeasible) throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_statusMessage_solve(s, status));
				;

				if (status != GLP_OPT) throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_statusMessage_solve(s, status));
				;
				s.out.solutionIsFeasible = true;
				s.out.solutionIsOptimal = true;
				s.out.primalCost = g.glp_ipt_obj_val(lp);
				int nCols = g.glp_get_num_cols(lp);
				if (nCols != s.in.numDecVariables) throw new JOMException("JOM - GLPK interface. Unexpected error");
				for (int cont = 1; cont <= nCols; cont++)
				{
					double multiplier = g.glp_ipt_col_dual(lp, cont);
					s.out.primalSolution.set(cont - 1, g.glp_ipt_col_prim(lp, cont));
					double[] bounds = glpk_obtainMultipliersLbUb(dvType.get(cont - 1), multiplier);
					s.out.multiplierOfLowerBoundConstraintToPrimalVariables.set(cont - 1, bounds[0]);
					s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(cont - 1, bounds[1]);
				}

				int nRows = g.glp_get_num_rows(lp);
				if (nRows != s.in.numConstraints) throw new JOMException("JOM - GLPK interface. Unexpected error");
				for (int cont = 1; cont <= nRows; cont++)
				{
					s.out.primalValuePerConstraint.set(cont - 1, g.glp_ipt_row_prim(lp, cont));
					s.out.multiplierOfConstraint.set(cont - 1, g.glp_ipt_row_dual(lp, cont));
				}
			} else
				throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_errorMessage_solve(s, ret));

			s.problemAlreadyAttemptedTobeSolved = true;
			g.glp_delete_prob(lp);
			return ret;
		}

		if ("simplex".equals(param.get("glpkSolverType")))
		{
			/* use simplex method */
			_JNA_GLPK.GLP_SMCP parm = new _JNA_GLPK.GLP_SMCP();
			g.glp_init_smcp(parm.getPointer());
			/* Set specific parameters specified by the user */
			for (Entry<String, Object> entry : this.param.entrySet())
			{
				String key = entry.getKey();
				if (key.equalsIgnoreCase("solverLibraryName")) continue; // a non-GPLK specific option
				if (key.equalsIgnoreCase("glpkSolverType")) continue; // a non-GPLK specific option
				Object value = entry.getValue();
				if (value instanceof Integer)
					parm.writeField(key , ((Integer) value).intValue());
				else if (value instanceof Double)
					parm.writeField(key , ((Double) value).doubleValue());
				else
					throw new JOMException("JOM - GLPK interface. Unknown value type in parameters: " + key);
			}
			int ret = g.glp_simplex(lp, parm.getPointer());
			if ((ret == 0) || (ret == GLP_ETMLIM))
			{
				int status = g.glp_get_status(lp);
				s.out.statusCode = status;
				s.out.statusMessage = glpk_statusMessage_solve(s, status);
				s.out.solutionIsFeasible = ((status == GLP_OPT) || (status == GLP_FEAS));
				s.out.solutionIsOptimal = (status == GLP_OPT);
				s.out.feasibleSolutionDoesNotExist = (status == GLP_NOFEAS);
				s.out.foundUnboundedSolution = (status == GLP_UNBND);
				if (!s.out.solutionIsFeasible) throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_statusMessage_solve(s, status));
				
				if ((status != GLP_OPT) && (status != GLP_ETMLIM)) throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_statusMessage_solve(s, status));
				s.out.solutionIsFeasible = true;
				s.out.solutionIsOptimal = true;
				s.out.solutionIsFeasible = true;
				s.out.solutionIsOptimal = true;
				s.out.primalCost = g.glp_get_obj_val(lp);
				int nCols = g.glp_get_num_cols(lp);
				if (nCols != s.in.numDecVariables) throw new JOMException("JOM - GLPK interface. Unexpected error"); 
				for (int cont = 1; cont <= nCols; cont++)
				{
					double multiplier = g.glp_get_col_dual(lp, cont);
					s.out.primalSolution.set(cont - 1, g.glp_get_col_prim(lp, cont));
					double[] bounds = glpk_obtainMultipliersLbUb(dvType.get(cont - 1), multiplier);
					s.out.multiplierOfLowerBoundConstraintToPrimalVariables.set(cont - 1, bounds[0]);
					s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(cont - 1, bounds[1]);
				}

				int nRows = g.glp_get_num_rows(lp);
				if (nRows != s.in.numConstraints) throw new JOMException("JOM - GLPK interface. Unexpected error");
				for (int cont = 1; cont <= nRows; cont++)
				{
					s.out.primalValuePerConstraint.set(cont - 1, g.glp_get_row_prim(lp, cont));
					s.out.multiplierOfConstraint.set(cont - 1, g.glp_get_row_dual(lp, cont));
				}
			} else
				throw new JOMException("JOM - GLPK interface. The problem could not be solved. GLPK error message: " + this.glpk_errorMessage_solve(s, ret));
			s.problemAlreadyAttemptedTobeSolved = true;
			g.glp_delete_prob(lp);
			return ret;
		}

		return -1;

	}
	private double[] glpk_obtainMultipliersLbUb(int varOrConstraintType, double multiplier)
	{
		double[] res = new double[2];
		if (varOrConstraintType == GLP_FR)
		{
			res[0] = 0;
			res[1] = 0;
		} else if (varOrConstraintType == GLP_LO)
		{
			res[0] = multiplier;
			res[1] = 0;
		} else if (varOrConstraintType == GLP_UP)
		{
			res[0] = 0;
			res[1] = multiplier;
		} else if (varOrConstraintType == GLP_FX)
		{
			res[0] = multiplier;
			res[1] = -multiplier;
		} else if (varOrConstraintType == GLP_DB) 
		{
			if (this.s.in.toMinimize)
			{
				// if multiplier is positive is of the lower bound 
				res[0] = (multiplier > 0)? multiplier : 0;
				res[1] = (multiplier <= 0)? multiplier : 0;
			}
			else
			{
				// if multiplier is positive is of the upper bound 
				res[0] = (multiplier < 0)? multiplier : 0;
				res[1] = (multiplier >= 0)? multiplier : 0;
			}
				 
		} else
			throw new JOMException("JOM - GLPK interface. Unexpected error");
		return res;

	}

}
