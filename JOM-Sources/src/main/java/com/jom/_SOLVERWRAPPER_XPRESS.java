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


import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.jet.math.tdouble.DoubleFunctions;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

/** @author Pablo */
@SuppressWarnings("unchecked")
class _SOLVERWRAPPER_XPRESS
{
	private final HashMap<String, Object> param;
	private final _INTERNAL_SolverIO      s;
	private final String                  solverLibraryName;
	private final Class c_XPRSconstants , c_XPRS , c_XPRSenumerations_ObjSense, c_XPRSprob;
	private final double XPRSconstants_MINUSINFINITY , XPRSconstants_PLUSINFINITY;
	private final int XPRSconstants_BESTBOUND , XPRSconstants_ERRORCODE , XPRSconstants_MAXTIME;
	private final int XPRSconstants_MIPSTATUS , XPRSconstants_MIP_OPTIMAL , XPRSconstants_MIP_SOLUTION;
	private final int XPRSconstants_MIP_INFEAS , XPRSconstants_MIP_UNBOUNDED , XPRSconstants_MIPBESTOBJVAL;
	private final int XPRSconstants_ROWS , XPRSconstants_COLS;
	private final int XPRSconstants_LPSTATUS , XPRSconstants_LP_OPTIMAL , XPRSconstants_LPOBJVAL;
	private final int XPRSconstants_PRIMALINFEAS , XPRSconstants_LP_UNFINISHED, XPRSconstants_LP_INFEAS , XPRSconstants_LP_UNBOUNDED;

	
	
	_SOLVERWRAPPER_XPRESS(_INTERNAL_SolverIO s, HashMap<String, Object> param)
	{
		this.s = s;
		this.solverLibraryName = (String) param.get("solverLibraryName");
		this.param = param;
		try
		{
			this.c_XPRSconstants = Class.forName("com.dashoptimization.XPRSconstants");
			this.c_XPRS = Class.forName("com.dashoptimization.XPRS");
			this.c_XPRSprob = Class.forName("com.dashoptimization.XPRSprob");
			this.c_XPRSenumerations_ObjSense = Class.forName("com.dashoptimization.XPRSenumerations$ObjSense");
			this.XPRSconstants_MINUSINFINITY = c_XPRSconstants.getField("MINUSINFINITY").getDouble(null);
			this.XPRSconstants_PLUSINFINITY = c_XPRSconstants.getField("PLUSINFINITY").getDouble(null);
			this.XPRSconstants_MAXTIME = c_XPRSconstants.getField("MAXTIME").getInt(null);
			this.XPRSconstants_BESTBOUND = c_XPRSconstants.getField("BESTBOUND").getInt(null);
			this.XPRSconstants_ERRORCODE = c_XPRSconstants.getField("ERRORCODE").getInt(null);
			this.XPRSconstants_MIPSTATUS = c_XPRSconstants.getField("MIPSTATUS").getInt(null);
			this.XPRSconstants_MIP_OPTIMAL = c_XPRSconstants.getField("MIP_OPTIMAL").getInt(null);
			this.XPRSconstants_MIP_SOLUTION = c_XPRSconstants.getField("MIP_SOLUTION").getInt(null);
			this.XPRSconstants_MIP_INFEAS = c_XPRSconstants.getField("MIP_INFEAS").getInt(null);
			this.XPRSconstants_MIPBESTOBJVAL = c_XPRSconstants.getField("MIPBESTOBJVAL").getInt(null);
			this.XPRSconstants_MIP_UNBOUNDED = c_XPRSconstants.getField("MIP_UNBOUNDED").getInt(null);
			this.XPRSconstants_ROWS = c_XPRSconstants.getField("ROWS").getInt(null);
			this.XPRSconstants_COLS = c_XPRSconstants.getField("COLS").getInt(null);

			this.XPRSconstants_LPSTATUS = c_XPRSconstants.getField("LPSTATUS").getInt(null);
			this.XPRSconstants_LP_OPTIMAL = c_XPRSconstants.getField("LP_OPTIMAL").getInt(null);
			this.XPRSconstants_LPOBJVAL = c_XPRSconstants.getField("LPOBJVAL").getInt(null);
			this.XPRSconstants_PRIMALINFEAS = c_XPRSconstants.getField("PRIMALINFEAS").getInt(null);
			this.XPRSconstants_LP_UNFINISHED = c_XPRSconstants.getField("LP_UNFINISHED").getInt(null);
			this.XPRSconstants_LP_INFEAS = c_XPRSconstants.getField("LP_INFEAS").getInt(null);
			this.XPRSconstants_LP_UNBOUNDED = c_XPRSconstants.getField("LP_UNBOUNDED").getInt(null);
		} catch (Exception e) { e.printStackTrace();throw new JOMException ("The classes of the EXPRESS solver cannot be found. Note that xprs.jar should be in the classpath, or accessible to the class loader."); }
	}

	int solve()
	{
		Object p_XPRSprob_obj = null;
		try
		{
			c_XPRS.getMethod("init", String.class).invoke(null , solverLibraryName); 
			//p = new XPRSprob();
			
			p_XPRSprob_obj = c_XPRSprob.newInstance();
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
				_dlb [col] = s.in.primalSolutionLowerBound.get(col) == -Double.MAX_VALUE? XPRSconstants_MINUSINFINITY : s.in.primalSolutionLowerBound.get(col);
				_dub [col] = s.in.primalSolutionUpperBound.get(col) == Double.MAX_VALUE? XPRSconstants_PLUSINFINITY : s.in.primalSolutionUpperBound.get(col);
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
			for (Entry<String, Object> entry : this.param.entrySet())
			{
				String keyStr = entry.getKey();
				if (keyStr.equalsIgnoreCase("solverLibraryName")) continue;
				if (keyStr.equals("maxSolverTimeInSeconds")) 
				{
					final Double val_maxSolverTimeInSeconds = ((Number) entry.getValue()).doubleValue();
					if (val_maxSolverTimeInSeconds > 0)
						c_XPRSprob.getMethod("setIntControl" , int.class , int.class).invoke(p_XPRSprob_obj , XPRSconstants_MAXTIME , (int) -Math.ceil(val_maxSolverTimeInSeconds));
					continue;
				}
				int key = new Integer(entry.getKey());
				Object value = entry.getValue();
				if (value instanceof String)
					c_XPRSprob.getMethod("setStrControl" , String.class , String.class).invoke(p_XPRSprob_obj , key , (String) value);
				else if (value instanceof Integer)
					c_XPRSprob.getMethod("setIntControl" , String.class , int.class).invoke(p_XPRSprob_obj , key , (Integer) value);
				else if (value instanceof Double)
					c_XPRSprob.getMethod("setDblControl" , String.class , double.class).invoke(p_XPRSprob_obj , key , (double) value);
				else
					throw new JOMException("JOM - XPRESS interface. Unknown value type in parameters");
			}
			
			
			if (s.in.hasIntegerVariables)
			{
				/* load the problem */
				c_XPRSprob.getMethod("loadGlobal",String.class,int.class,int.class,byte[].class,double[].class, double[].class,double[].class, int[].class, int[].class, int[].class, double[].class, double[].class, double[].class, int.class, int.class, byte[].class, int[].class, double[].class, byte[].class, int[].class, int[].class, double[].class).invoke(p_XPRSprob_obj , "" , ncols, nrows, _srowtypes, _drhs, _drange, _dobj, return_mstart [0],
						return_mnel [0], return_mrwind [0], return_dmatval [0], _dlb, _dub, ngents, nsets, _qgtype, _mgcols, _dlim,
						_stype, _msstart, _mscols, _dref);
				/* Minimize or maximize */
				c_XPRSprob.getMethod("chgObjSense" , c_XPRSenumerations_ObjSense).invoke(p_XPRSprob_obj , s.in.toMinimize? c_XPRSenumerations_ObjSense.getField("MINIMIZE").get(null) : c_XPRSenumerations_ObjSense.getField("MAXIMIZE").get(null));

				/* Call the solver */
				c_XPRSprob.getMethod("mipOptimize").invoke(p_XPRSprob_obj);

				s.problemAlreadyAttemptedTobeSolved = true;
				s.out.bestOptimalityBound = ((double) c_XPRSprob.getMethod("getDblAttrib" , int.class).invoke(p_XPRSprob_obj , XPRSconstants_BESTBOUND));
				s.out.statusCode = ((int) c_XPRSprob.getMethod("getIntAttrib" , int.class).invoke(p_XPRSprob_obj , XPRSconstants_ERRORCODE));
				s.out.statusMessage = (String) c_XPRSprob.getMethod("getLastError").invoke(p_XPRSprob_obj);

				final int mipStatus = (int) c_XPRSprob.getMethod("getIntAttrib",int.class).invoke(p_XPRSprob_obj , XPRSconstants_MIPSTATUS);
				s.out.solutionIsOptimal = (mipStatus == XPRSconstants_MIP_OPTIMAL);
				s.out.solutionIsFeasible = s.out.solutionIsOptimal || (mipStatus == XPRSconstants_MIP_SOLUTION);
				s.out.feasibleSolutionDoesNotExist = (mipStatus == XPRSconstants_MIP_INFEAS);
				s.out.foundUnboundedSolution = (mipStatus == XPRSconstants_MIP_UNBOUNDED);

				/* I may have a bound even if a feasible solution was not found */
				s.out.primalCost = (double) c_XPRSprob.getMethod("getDblAttrib",int.class).invoke(p_XPRSprob_obj , XPRSconstants_MIPBESTOBJVAL);

				if (!s.out.solutionIsFeasible)
				{
					try { c_XPRS.getMethod("free").invoke(null); } catch (Exception ee) {} // frees the license, can fail if already done
					return s.out.statusCode;
				}
				/* Check the number of constraitns and variables */
				if (((int) c_XPRSprob.getMethod("getIntAttrib" , int.class).invoke(p_XPRSprob_obj , XPRSconstants_ROWS)) != s.in.numConstraints) throw new JOMException("JOM - XPRESS interface. Unexpected error");
				if (((int) c_XPRSprob.getMethod("getIntAttrib" , int.class).invoke(p_XPRSprob_obj , XPRSconstants_COLS))   != s.in.numDecVariables) throw new JOMException("JOM - XPRESS interface. Unexpected error");

				/* Retrieve the optimal primal solution */
				double [] primalSolution = new double [s.in.numDecVariables];
				double [] slackSolution = new double [s.in.numConstraints];
				c_XPRSprob.getMethod("getMipSol" , double[].class,double[].class).invoke(p_XPRSprob_obj , primalSolution , slackSolution);
				s.out.primalSolution = DoubleFactory1D.dense.make(primalSolution);

				/* Retrieve the values of the constraints in the solution */
				double[] rhsCplex = (s.in.numConstraints == 0) ? new double[0] :
						Arrays.copyOf(s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray(), s.in
								.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray().length);
				for (int cont = 0; cont < rhsCplex.length; cont++) rhsCplex[cont] = -rhsCplex[cont];
				
				double[] slack = new double[s.in.numConstraints];
				c_XPRSprob.getMethod("calcSlacks" , double[].class,double[].class).invoke(p_XPRSprob_obj , primalSolution , slack);
				s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(rhsCplex).assign(DoubleFactory1D.dense.make(slack), DoubleFunctions.minus);
				s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(s.in.numConstraints);
				s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
				s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);

				try { c_XPRS.getMethod("free").invoke(null); } catch (Exception ee) {} // frees the license, can fail if already done
				return s.out.statusCode;
			}
			else
			{
				/* load the problem */
				
				c_XPRSprob.getMethod("loadLp" , String.class, int.class, int.class, byte[].class, double[].class, double[].class, double[].class, int[].class, int[].class, int[].class, 
						double[].class, double[].class, double[].class).invoke(p_XPRSprob_obj , "" , 
								ncols, nrows, _srowtypes, _drhs, _drange, _dobj, return_mstart [0], return_mnel [0], 
								return_mrwind [0],return_dmatval [0], _dlb, _dub);
				/* Minimize or maximize */
				c_XPRSprob.getMethod("chgObjSense" , c_XPRSenumerations_ObjSense).invoke(p_XPRSprob_obj , s.in.toMinimize? c_XPRSenumerations_ObjSense.getField("MINIMIZE").get(null) : c_XPRSenumerations_ObjSense.getField("MAXIMIZE").get(null));

				/* Call the solver */
				c_XPRSprob.getMethod("lpOptimize").invoke(p_XPRSprob_obj);

				s.problemAlreadyAttemptedTobeSolved = true;
				s.out.bestOptimalityBound = s.in.toMinimize? -Double.MAX_VALUE : Double.MAX_VALUE; //p.getDblAttrib(XPRSconstants.BESTBOUND); // pablo: this may not be the one for this
				s.out.statusCode = (int) c_XPRSprob.getMethod("getIntAttrib",int.class).invoke(p_XPRSprob_obj , XPRSconstants_ERRORCODE);
				s.out.statusMessage = (String) c_XPRSprob.getMethod("getLastError").invoke(p_XPRSprob_obj); 

				final int lpStatus = (int) c_XPRSprob.getMethod("getIntAttrib",int.class).invoke(p_XPRSprob_obj , XPRSconstants_LPSTATUS);
				s.out.solutionIsOptimal = (lpStatus == XPRSconstants_LP_OPTIMAL);
				final int a = (int) c_XPRSprob.getMethod("getIntAttrib",int.class).invoke(p_XPRSprob_obj , XPRSconstants_PRIMALINFEAS); 
				s.out.solutionIsFeasible = s.out.solutionIsOptimal || ((lpStatus == XPRSconstants_LP_UNFINISHED) && (a == 0));  
				s.out.feasibleSolutionDoesNotExist = (lpStatus == XPRSconstants_LP_INFEAS);
				s.out.foundUnboundedSolution = (lpStatus == XPRSconstants_LP_UNBOUNDED);

				if (!s.out.solutionIsFeasible)
				{
					try { c_XPRS.getMethod("free").invoke(null); } catch (Exception ee) {} // frees the license, can fail if already done
					return s.out.statusCode;
				}

				/* Check the number of constraitns and variables */
				if (((int) c_XPRSprob.getMethod("getIntAttrib" , int.class).invoke(p_XPRSprob_obj , XPRSconstants_ROWS)) != s.in.numConstraints) throw new JOMException("JOM - XPRESS interface. Unexpected error");
				if (((int) c_XPRSprob.getMethod("getIntAttrib" , int.class).invoke(p_XPRSprob_obj , XPRSconstants_COLS))   != s.in.numDecVariables) throw new JOMException("JOM - XPRESS interface. Unexpected error");

				/* Retrieve the optimal primal solution */
				final double [] primalSolution = new double [s.in.numDecVariables];
				double [] slackSolution = new double [s.in.numConstraints];
				double [] mulipliersSolution = new double [s.in.numConstraints];
				double [] reducedCosts = new double [s.in.numDecVariables];
				c_XPRSprob.getMethod("getLpSol",double[].class,double[].class,double[].class,double[].class).invoke(p_XPRSprob_obj,primalSolution , slackSolution , mulipliersSolution , reducedCosts);
				s.out.primalSolution = DoubleFactory1D.dense.make(primalSolution);
				//DoubleHolder obj = new DoubleHolder(); p.calcObjective(primalSolution , obj);
				s.out.primalCost = (double) c_XPRSprob.getMethod("getDblAttrib",int.class).invoke(p_XPRSprob_obj, XPRSconstants_LPOBJVAL);//s.in.objectiveFunction.evaluate_internal(primalSolution).toValue();
				//s.out.bestOptimalityBound;
				
				
				/* Retrieve the values of the constraints in the solution */
				double[] rhsCplex = (s.in.numConstraints == 0) ? new double[0] :
						Arrays.copyOf(s.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray(), s.in
								.lhsMinusRhsAccumulatedConstraint.getAffineExpression().getConstantCoefArray().length);
				for (int cont = 0; cont < rhsCplex.length; cont++) rhsCplex[cont] = -rhsCplex[cont];
				
				double[] slack = new double[s.in.numConstraints];
				c_XPRSprob.getMethod("calcSlacks" , double[].class,double[].class).invoke(p_XPRSprob_obj , primalSolution , slack);
				s.out.primalValuePerConstraint = DoubleFactory1D.dense.make(rhsCplex).assign(DoubleFactory1D.dense.make(slack), DoubleFunctions.minus);
				s.out.multiplierOfConstraint = DoubleFactory1D.dense.make(mulipliersSolution);
				s.out.multiplierOfLowerBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
				s.out.multiplierOfUpperBoundConstraintToPrimalVariables = DoubleFactory1D.dense.make(s.in.numDecVariables);
				for (int dv = 0; dv < s.in.numDecVariables; dv++)
				{
					double lb = s.in.primalSolutionLowerBound.get(dv);
					double ub = s.in.primalSolutionUpperBound.get(dv);
					double val = s.out.primalSolution.get(dv);
					if (lb == ub)
					{
						s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(dv, reducedCosts[dv]);
						s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(dv, -reducedCosts[dv]);
					} else if (Math.abs(val - lb) > Math.abs(val - ub))
						s.out.multiplierOfUpperBoundConstraintToPrimalVariables.set(dv, reducedCosts[dv]);
					else
						s.out.multiplierOfLowerBoundConstraintToPrimalVariables.set(dv, reducedCosts[dv]);
				}

				try { c_XPRS.getMethod("free").invoke(null); } catch (Exception ee) {} // frees the license, can fail if already done
				return s.out.statusCode;
			}
	
		} catch (Exception e)
		{
			try { c_XPRS.getMethod("free").invoke(null); } catch (Exception ee) {} // frees the license, can fail if already done
			e.printStackTrace();
			System.out.println(e.getLocalizedMessage());
			try { if (p_XPRSprob_obj != null) c_XPRSprob.getMethod("destroy").invoke(p_XPRSprob_obj); } catch (Exception ee) {} // frees any memory associated to the object
			throw new JOMException (e.getLocalizedMessage());
		} 
		
	}

}
