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

/* PABLO: CHANGE POM TO COMPILE: THE JAR SHOULD BE LOCAL? SOME FORM TO COMPILE WITHOUT THE JAR?
 * VEr: http://stackoverflow.com/questions/8325819/how-can-i-create-an-executable-jar-without-dependencies-using-maven
 *   */

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

import com.dashoptimization.XPRS;
import com.dashoptimization.XPRSconstants;
import com.dashoptimization.XPRSenumerations;
import com.dashoptimization.XPRSprob;

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
            final double[] _drhs = (s.in.numConstraints == 0) ? new double[0] : 
    				Arrays.copyOf(s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray(), s.in
    						.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray().length);
            final byte[] _srowtypes = new byte [s.in.numConstraints]; Arrays.fill(_srowtypes , (byte) 'R');
    		for (int contC = 0; contC < _srowtypes.length; contC++)
    		{
    			_drhs[contC] = -_drhs[contC];
    			final double lb = s.in.constraintLowerBound.get(contC); // 0 or -Double.MAX_VALUE
    			final double ub = s.in.constraintUpperBound.get(contC); // 0 or Double.MAX_VALUE
    			if (lb == ub)
    				_srowtypes[contC] = 'E';
    			else if ((lb == -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
    				throw new JOMException("JOM - XPRESS interface. Unbounded contraints are not allowed"); // cType.set(contC,GLP_FR);
    			else if ((lb != -Double.MAX_VALUE) && (ub == Double.MAX_VALUE))
    				_srowtypes[contC] = 'G';
    			else if ((lb == -Double.MAX_VALUE) && (ub != Double.MAX_VALUE))
    				_srowtypes[contC] = 'L';
    			else
    				throw new JOMException("JOM - XPRESS interface. Double bounded contraints are supposed not to exist in JOM");
    		}
            final double [] _drange = null;
            final double[] _dobj = s.in.objectiveFunction.getAffineExpression().getCellLinearCoefsFull(0);
			_INTERNAL_AffineExpressionCoefs constraintsMatrix = s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression();
			int[][] return_mnel = new int [1][]; 
			int[][] return_mrwind = new int [1][];
			int[][] return_mstart = new int [1][];
			double[][] return_dmatval = new double [1][];
			constraintsMatrix.getNonZerosRowColValForXpressSolver(return_mnel, return_mrwind, return_mstart, return_dmatval, ncols);
			final double [] _dlb = new double [ncols];
			final double [] _dub = new double [ncols];
			for (int col = 0 ; col < ncols ; col ++)
			{
				_dlb [col] = s.in.primalSolutionLowerBound.get(col) == -Double.MAX_VALUE? XPRSconstants.MINUSINFINITY : s.in.primalSolutionLowerBound.get(col);
				_dub [col] = s.in.primalSolutionUpperBound.get(col) == Double.MAX_VALUE? XPRSconstants.PLUSINFINITY : s.in.primalSolutionUpperBound.get(col);
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
			
			/* Set the environment parameters, including the maxSolverTime */
			System.out.println(this.param);
			for (Entry<String, Object> entry : this.param.entrySet())
			{
				String keyStr = entry.getKey();
				if (keyStr.equalsIgnoreCase("solverLibraryName")) continue;
				int key = new Integer(entry.getKey());
				Object value = entry.getValue();
				if (value instanceof String)
					p.setStrControl(key , (String) value);
				else if (value instanceof Integer)
					p.setIntControl(key , (Integer) value);
				else if (value instanceof Double)
				{
					System.out.println("here");
					p.setDblControl(key , (double) value);
				}
				else
					throw new JOMException("JOM - XPRESS interface. Unknown value type in parameters");
			}
			
			
			if (s.in.hasIntegerVariables)
			{
				/* load the problem */
				p.loadGlobal("" , ncols, nrows, _srowtypes, _drhs, _drange, _dobj, return_mstart [0],
						return_mnel [0], return_mrwind [0], return_dmatval [0], _dlb, _dub, ngents, nsets, _qgtype, _mgcols, _dlim,
						_stype, _msstart, _mscols, _dref);
				/* Minimize or maximize */
				p.chgObjSense(s.in.toMinimize? XPRSenumerations.ObjSense.MINIMIZE : XPRSenumerations.ObjSense.MAXIMIZE);

				/* Call the solver */
				p.mipOptimize();

				s.problemAlreadyAttemptedTobeSolved = true;
				s.out.bestOptimalityBound = p.getDblAttrib(XPRSconstants.BESTBOUND);
				s.out.statusCode = p.getIntAttrib(XPRSconstants.ERRORCODE);
				s.out.statusMessage = p.getLastError();

				final int mipStatus = p.getIntAttrib(XPRSconstants.MIPSTATUS);
				s.out.solutionIsOptimal = (mipStatus == XPRSconstants.MIP_OPTIMAL);
				s.out.solutionIsFeasible = s.out.solutionIsOptimal || (mipStatus == XPRSconstants.MIP_SOLUTION);
				s.out.feasibleSolutionDoesNotExist = (mipStatus == XPRSconstants.MIP_INFEAS);
				s.out.foundUnboundedSolution = (mipStatus == XPRSconstants.MIP_UNBOUNDED);

				if (!s.out.solutionIsFeasible)
					return s.out.statusCode;

				s.out.primalCost = p.getDblAttrib(XPRSconstants.MIPBESTOBJVAL);

				/* Check the number of constraitns and variables */
				if (p.getIntAttrib(XPRSconstants.ROWS) != s.in.numConstraints) throw new JOMException("JOM - CPLEX interface. Unexpected error");
				if (p.getIntAttrib(XPRSconstants.COLS)  != s.in.numDecVariables) throw new JOMException("JOM - CPLEX interface. Unexpected error");

				/* Retrieve the optimal primal solution */
				double [] primalSolution = new double [s.in.numDecVariables];
				double [] slackSolution = new double [s.in.numConstraints];
				p.getMipSol(primalSolution , slackSolution);
				s.out.primalSolution = DoubleFactory1D.dense.make(primalSolution);

				/* Retrieve the values of the constraints in the solution */
				double[] rhsCplex = (s.in.numConstraints == 0) ? new double[0] :
						Arrays.copyOf(s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray(), s.in
								.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray().length);
				for (int cont = 0; cont < rhsCplex.length; cont++) rhsCplex[cont] = -rhsCplex[cont];
				
				double[] slack = new double[s.in.numConstraints];
				p.calcSlacks(primalSolution , slack);
				s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(rhsCplex).assign(DoubleFactory1D.dense.make(slack), DoubleFunctions.minus);
				s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(s.in.numConstraints);
				s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
				s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
				
				return s.out.statusCode;
			}
			else
			{
				 // pag 454: LPSTATUS --> infreasible, optimal, unbounded... 
				 // pag 457: MIPSTATUS --> infreasible, optimal, unbounded... 
				 // pag 467: STOPSTATUS --> stopped by time limit...
				 
				 // pag. 369: BARGAPSTOP, BARPRIMALSTOP, FEASTOL, FEASTOLTARGET: stop conditions
				 // PAG. 410: MAXTIME: maximum solver time
				 // MIPABSSTOP, MIPRELSTOP: stop criterion based on absolute gap
				throw new RuntimeException ("Only MIPs");
			}
	
		} catch (RuntimeException e)
		{
			XPRS.free(); // frees the license
			if (p != null) p.destroy(); // frees any memory associated to the object
			e.printStackTrace(); throw e;
		} 
		
	}

}
