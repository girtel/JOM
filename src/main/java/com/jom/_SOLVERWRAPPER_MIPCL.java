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
import java.util.Map;

import org.apache.commons.lang.SystemUtils;

import by.bsu.JVmipcl.LP;
import by.bsu.JVmipcl.MIP;
import cern.colt.matrix.tdouble.DoubleFactory1D;

/** @author Pablo */
@SuppressWarnings("unchecked")
class _SOLVERWRAPPER_MIPCL
{
	private final Map<String, Object> param;
	private final _INTERNAL_SolverIO      s;

	_SOLVERWRAPPER_MIPCL(_INTERNAL_SolverIO s, Map<String, Object> param)
	{
		this.s = s;
		this.param = param;
	}

	int solve()
	{
		final boolean isIlp = s.in.hasIntegerVariables;
		final LP prob = isIlp? new MIP ("") : new LP(""); // create a new MIP problem
		try
		{
            prob.beSilent(true);

			_INTERNAL_AffineExpressionCoefs constraintsMatrix = s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression();
			int[][] rows = new int[1][];
			int[][] cols = new int[1][];
			double[][] vals = new double[1][];
			constraintsMatrix.getNonZerosRowColVal(rows, cols, vals);
			final int[] ia = rows[0];
			final int[] ja = cols[0];
			final double[] ar = vals[0];
			
            final int ncols = s.in.numDecVariables;
            final int nrows = s.in.numConstraints; 
            final int numNonZeros = ar.length;
            prob.openMatrix(nrows, ncols, numNonZeros);

            /* Set the minimization/maximization sense */
            prob.setObjSense(!s.in.toMinimize);
            
            /* Add the vars with the objective coefficients */
            final double[] objectiveFunctionCoefs = s.in.objectiveFunction.getAffineExpression().getCellLinearCoefsFull(0);
			for (int col = 0 ; col < ncols ; col ++)
			{
				final double dlb = s.in.primalSolutionLowerBound.get(col) == -Double.MAX_VALUE? -LP.VAR_INF : s.in.primalSolutionLowerBound.get(col);
				final double dub = s.in.primalSolutionUpperBound.get(col) == Double.MAX_VALUE? LP.VAR_INF : s.in.primalSolutionUpperBound.get(col);
				final int type = s.in.primalSolutionIsInteger.get(col) != 0? (dlb == 0 && dub == 1? MIP.VAR_BIN : MIP.VAR_INT) : 0;
				prob.addVar(col, type, objectiveFunctionCoefs [col], dlb, dub);
			}
            
            /* Add empty constraints, with LHS and RHS */
			for (int row = 0 ; row < nrows ; row ++)
			{
				double lhs = s.in.constraintLowerBound.get(row); // 0 or -Double.MAX_VALUE
				double rhs = s.in.constraintUpperBound.get(row); // 0 or Double.MAX_VALUE
				if (lhs == rhs)
				{
					lhs = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(row); // equal to minus the constant value
					rhs = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(row); // equal to minus the constant value
				} else if ((lhs == -Double.MAX_VALUE) && (rhs == Double.MAX_VALUE))
				{
					throw new JOMException("JOM - GLPK interface. Unbounded contraints are not allowed"); // cType.set(contC,GLP_FR);
				} else if ((lhs != -Double.MAX_VALUE) && (rhs == Double.MAX_VALUE))
				{
					rhs = LP.INF;
					lhs = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(row); // equal to minus the constant value
				} else if ((lhs == -Double.MAX_VALUE) && (rhs != Double.MAX_VALUE))
				{
					lhs = -LP.INF;
					rhs = -s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getCellConstantCoef(row); // equal to minus the constant value
				} else
				{
					throw new JOMException("JOM - GLPK interface. Double bounded constraints are supposed not to exist in JOM");
				}
				prob.addCtr(row, 0, lhs, rhs);
			}
            
            /* Add the entries of the matrix */
            for (int cont = 0; cont < numNonZeros ; cont ++)
            	prob.addEntry(ar [cont], ia [cont], ja [cont]);
            
            // prob.preprocOff();

            /* Closes the matrix, needed for preparing it for solving */
            prob.closeMatrix(); // close the matrix

            /* solves the problem */
            if (isIlp)
            {
    			/* Set the environment parameters, including the maxSolverTime */
                final long maxSolverTimeInSeconds = ((Number) this.param.getOrDefault("maxSolverTimeInSeconds", new Double (10000000))).longValue();
                final double dualityGapLimitInAbsoluteValue = ((Number) this.param.getOrDefault("dualityGapLimitInAbsoluteValue", new Double (0.0))).doubleValue();
            	((MIP) prob).optimize(maxSolverTimeInSeconds , dualityGapLimitInAbsoluteValue, SystemUtils.IS_OS_WINDOWS? "NUL" : (SystemUtils.IS_OS_LINUX? "/dev/null" : null));
//            	((MIP) prob).optimize(maxSolverTimeInSeconds , dualityGapLimitInAbsoluteValue, null); // this crashes JVM
            	
            }
            else
            	prob.optimize();
            
            /* solves the problem */

			s.problemAlreadyAttemptedTobeSolved = true;
			s.out.bestOptimalityBound = s.in.toMinimize? -Double.MAX_VALUE : Double.MAX_VALUE; 
			s.out.statusCode = 0;
			s.out.statusMessage = "";
			s.out.solutionIsFeasible = prob.isSolution();
			s.out.solutionIsOptimal = isIlp? ((MIP) prob).isSolutionOptimal() : s.out.solutionIsFeasible;
			s.out.feasibleSolutionDoesNotExist = isIlp? ((MIP)prob).isInfeasible() : prob.isLpInfeasible();
			s.out.foundUnboundedSolution = prob.isLpUnbounded();

			/* I may have a bound even if a feasible solution was not found */
			s.out.primalCost = prob.getObjVal();

			if (!s.out.solutionIsFeasible)
			{
				prob.dispose();
				return s.out.statusCode;
			}

			/* Retrieve the optimal primal solution */
			double [] primalSolution = new double [s.in.numDecVariables];
			prob.getSolution(); // makes MIPCL to update internal info, needed before retrieving the solution
			final double[] solutionVectorIndexedByVarHandles = prob.getX(); assert solutionVectorIndexedByVarHandles.length == s.in.numDecVariables;
			final int[] varHandles = prob.getVarHandles(); assert varHandles.length == s.in.numDecVariables;
			for (int i=0; i < s.in.numDecVariables ; ++i) 
				primalSolution [i] = solutionVectorIndexedByVarHandles [varHandles [i]];
			
			s.out.primalSolution = DoubleFactory1D.dense.make(primalSolution);

			/* Retrieve the values of the constraints in the solution */
			double[] rhsCplex = (s.in.numConstraints == 0) ? new double[0] :
					Arrays.copyOf(s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray(), s.in
							.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray().length);
			for (int cont = 0; cont < rhsCplex.length; cont++) rhsCplex[cont] = -rhsCplex[cont];
				
			s.out.primalValuePerConstraint = DoubleFactory1D.sparse.make(s.in.numConstraints); // we do not compute this... too costly
			s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(s.in.numConstraints);
			s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
			s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);

			if (!s.in.hasIntegerVariables)
			{
				final double [] reducedCostsIndexedByVarHandle_var = prob.getReducedCosts();
				if (reducedCostsIndexedByVarHandle_var != null)
					for (int i=0; i < s.in.numDecVariables ; ++i) 
					{
						final double reducedCost = Math.abs(reducedCostsIndexedByVarHandle_var [varHandles [i]]);
						if (Math.abs(s.out.primalSolution.get(i) - s.in.primalSolutionLowerBound.get(i)) < 1e-3)
							s.out.multiplierOfLowerBoundConstraintToPrimalVariables.set(i, s.in.toMinimize? reducedCost : -reducedCost);
						if (Math.abs(s.out.primalSolution.get(i) - s.in.primalSolutionUpperBound.get(i)) < 1e-3)
							s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(i, s.in.toMinimize? -reducedCost : reducedCost);
					}
				final double [] shadowPricesIndexedByCtrHandler = prob.getShadowPrices();
				final int[] constrHandles = prob.getCtrHandles();
				if (shadowPricesIndexedByCtrHandler != null && constrHandles != null)
					for (int i=0; i < s.in.numConstraints ; ++i) 
					{
						final double multiplier = Math.abs(shadowPricesIndexedByCtrHandler [constrHandles [i]]);
						s.out.multiplierOfConstraint.set (i , multiplier);
					}
			}
			prob.dispose();

			return s.out.statusCode;
		} catch (Exception e)
		{
			if (prob != null) prob.dispose();
			e.printStackTrace();
			System.out.println(e.getLocalizedMessage());
			throw new JOMException (e.getLocalizedMessage());
		} 
		
	}

}
