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

import cern.colt.matrix.tdouble.DoubleMatrix2D;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.LinkedHashMap;

/** @author Pablo */
class _FUNCTION_LINEAREXPRESSION extends Expression
{
	private static ArrayList<Long> times = new ArrayList<Long>();
	private final boolean isConstant;

	_FUNCTION_LINEAREXPRESSION(OptimizationProblem model, double constant)
	{
		super(model, new int[]{1, 1}, new _INTERNAL_AffineExpressionCoefs(model, constant));
		this.isConstant = this.affineExp.isConstant();
	}

	_FUNCTION_LINEAREXPRESSION(OptimizationProblem model, int[] size, _INTERNAL_AffineExpressionCoefs affineExpressionObject) // double
	{
		super(model, size, affineExpressionObject);
		this.isConstant = this.affineExp.isConstant();
	}

	_FUNCTION_LINEAREXPRESSION(OptimizationProblem model, _INTERNAL_DecisionVariableArray dv)
	{
		this(model, dv.getVarIds());
	}

	_FUNCTION_LINEAREXPRESSION(OptimizationProblem model, DoubleMatrixND constant)
	{
		super(model, constant.getSize().toArray(), new _INTERNAL_AffineExpressionCoefs(model, constant.getSize().toArray(), constant.elements()
				.toArray()));
		this.isConstant = this.affineExp.isConstant();
	}

	_FUNCTION_LINEAREXPRESSION(OptimizationProblem model, IntMatrixND varIds)
	{
		super(model, varIds.getSize().toArray(), new _INTERNAL_AffineExpressionCoefs(model, varIds));
		this.isConstant = this.affineExp.isConstant();
	}

	private static void tic(String msg)
	{
		times.add(System.nanoTime());
		if (times.size() > 1)
			System.out.println("" + times.size() + ": " + msg + ": " + (times.get(times.size() - 1) - times.get(times.size() - 2)) / 1E9);
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		throw new JOMException("Function 'linear expression': Unexpected error");
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		throw new JOMException("Function 'linear expression': Unexpected error");
	}

	/* (non-Javadoc)
	 * @see Expression#isLinear()
	 */
	@Override
	boolean isLinear()
	{
		return true;
	}

	/* (non-Javadoc)
	 * @see Expression#isConstant()
	 */
	@Override
	boolean isConstant()
	{
		return this.isConstant;
	}

	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		throw new JOMException("JOM. Function 'linear expression': Unexpected error");
	}

}
