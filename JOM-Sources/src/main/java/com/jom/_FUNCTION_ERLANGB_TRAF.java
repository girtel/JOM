/*******************************************************************************
 * Copyright (c) 2015 Pablo Pavon Mari�o.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * <p>
 * Contributors:
 * Pablo Pavon Mari�o - initial API and implementation
 ******************************************************************************/

package com.jom;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;

class _FUNCTION_ERLANGB_TRAF extends Expression
{
	static double[][] precomputedBlocking = new double[100][101]; // [numServers-1 , percentageUtilization]
	private final Expression loads;
	private final int[]      numServers;
	private final boolean    isConstant;
	private final boolean    isLinear;

	_FUNCTION_ERLANGB_TRAF(OptimizationProblem model, Expression loads, Expression servers)
	{
		super(model);

		if (!Arrays.equals(servers.getSize(), loads.getSize()))
			throw new JOMException("Function 'erlangB': The array of number of servers and the array of number of loads must have the same size");
		if (!servers.isConstant()) throw new JOMException("Function 'erlangB': The number of servers must be a constant");
		this.resize(loads.getSize());

		double[] aux_numServers = servers.evaluateConstant().elements().toArray();
		this.numServers = new int[aux_numServers.length];
		for (int cont = 0; cont < aux_numServers.length; cont++) this.numServers[cont] = (int) aux_numServers[cont];

		this.isLinear = loads.isConstant();
		this.isConstant = loads.isConstant();
		if (isLinear)
		{
			/* change (destroy) the constant data in Expression a. This object keeps no reference to Expression a, but the data remains since we 
			 * keep reference to the affine expression object of a, now modified by this function */
			double[] data = loads.getAffineExpression().getConstantCoefArray();
			for (int cont = 0; cont < data.length; cont++) data[cont] = erlangBLossProbability(data[cont], this.numServers[cont]);
			this.affineExp = new _INTERNAL_AffineExpressionCoefs(model, servers.getSize(), data); //a.getAffineExpression();
		}
		this.loads = (isLinear()) ? null : loads; // if linear do not store the expressions. This info is already in the linear coefs matrix

		/* precompute the blockings for these amount of servers */
		for (int contNumServer = 0; contNumServer < numServers.length; contNumServer++)
		{
			int c = numServers[contNumServer];
			if (c > precomputedBlocking.length) continue; // too many servers
			for (int contUtil = 0; contUtil < precomputedBlocking[0].length; contUtil++)
				precomputedBlocking[c - 1][contUtil] = erlangBLossProbability(c * contUtil * 0.01, c);
		}

	}

	private static double erlangBLossProbability(double load, int numberOfServers)
	{
		if (load < 0) throw new JOMException("Function 'erlangB': Erlang function is undefined for negative loads: " + load);

		// Initialize output variable
		double gradeOfService = 0;

		// Compute the grade of service (or Erlang loss probability)
		if (load > 1E-10)
		{
			double s = 0;
			for (int i = 1; i <= numberOfServers; i++)
				s = (1 + s) * (i / load);
			gradeOfService = 1 / (1 + s);
		}

		return gradeOfService;
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrix1D loadValues = this.loads.evaluate_internal(valuesDVs).elements();
		return new DoubleMatrixND(this.getSize(), erlangBLossProbability_bulk(loadValues));
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix1D loads_value = this.loads.evaluate_internal(valuesDVs).elements();
		DoubleMatrix2D loads_jacobian = this.loads.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix1D gradientErlangBRespectToLoads = DoubleFactory1D.dense.make(this.getNumScalarExpressions());
		for (int contCell = 0; contCell < this.getNumScalarExpressions(); contCell++)
			gradientErlangBRespectToLoads.set(contCell, partialDerivative(numServers[contCell], loads_value.get(contCell)));

		return DoubleFactory2D.sparse.diagonal(gradientErlangBRespectToLoads).zMult(loads_jacobian, null);
	}

	/* (non-Javadoc)
	 * @see com.jom.Expression#isLinear()
	 */
	@Override
	boolean isLinear()
	{
		return isLinear;
	}

	/* (non-Javadoc)
	 * @see com.jom.Expression#isConstant()
	 */
	@Override
	boolean isConstant()
	{
		return isConstant;
	}

	/* (non-Javadoc)
	 * @see com.jom.Expression#getActiveVarIds()
	 */
	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		return this.loads.getActiveVarIds();
	}

	private DoubleMatrix1D erlangBLossProbability_bulk(DoubleMatrix1D loads)
	{
		DoubleMatrix1D blockings = DoubleFactory1D.dense.make(this.getNumScalarExpressions());
		for (int contCell = 0; contCell < this.getNumScalarExpressions(); contCell++)
			blockings.set(contCell, erlangBLossProbability(loads.get(contCell), this.numServers[contCell]));
		return blockings;
	}

	private double partialDerivative(int numServers, double load)
	{
		double loadVariation = numServers * 0.01;
		double prevLoad = (load - loadVariation < 0) ? 0 : load - loadVariation;
		double postLoad = load + loadVariation;
		double blocking_prev = erlangBLossProbability(prevLoad, numServers);
		double blocking_post = erlangBLossProbability(postLoad, numServers);
		return (blocking_post - blocking_prev) / (postLoad - prevLoad);
	}

}
