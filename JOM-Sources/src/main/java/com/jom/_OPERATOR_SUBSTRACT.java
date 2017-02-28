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

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

/** @author Pablo */
@SuppressWarnings("unchecked")
class _OPERATOR_SUBSTRACT extends Expression
{
	private final Expression a;
	private final Expression b;
	private final boolean    isLinear;
	private final boolean    isConstant;

	_OPERATOR_SUBSTRACT(OptimizationProblem model, Expression a, Expression b)
	{
		super(model);
		
		/* compute size */
		if (!Arrays.equals(a.getSize(), b.getSize()))
			throw new JOMException("Operator '-'. Size mismatch: a and b should be of the same size, or one of them should be a scalar expression");
		this.resize(a.getSize());

		this.isLinear = (a.isLinear() && b.isLinear());
		this.affineExp = (isLinear) ? a.getAffineExpression().operator_substract(b.getAffineExpression()) : null;
		this.isConstant = (isLinear) ? this.affineExp.isConstant() : false;

		this.a = (this.isLinear()) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.b = (this.isLinear()) ? null : b; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrix1D a_value = a.evaluate_internal(valuesDVs).elements();
		DoubleMatrix1D b_value = b.evaluate_internal(valuesDVs).elements();
		a_value.assign(b_value, DoubleFunctions.minus);
		return new DoubleMatrixND(this.getSize(), a_value);
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix2D a_gradient = a.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D b_gradient = b.evaluateJacobian_internal(valuesDVs);
		return a_gradient.assign(b_gradient, DoubleFunctions.minus);
	}

	@Override
	boolean isLinear()
	{
		return this.isLinear;
	}

	/* (non-Javadoc)
	 * @see com.jom.Expression#isConstant()
	 */
	@Override
	boolean isConstant()
	{
		return this.isConstant;
	}

	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		/* Make a copy since I am modifying it */
		LinkedHashMap<Integer, HashSet<Integer>> newActiveVarIdsMatrix = new LinkedHashMap<Integer, HashSet<Integer>>();
		for (Entry<Integer, HashSet<Integer>> e : this.a.getActiveVarIds().entrySet())
			newActiveVarIdsMatrix.put(e.getKey(), (HashSet<Integer>) e.getValue().clone());
		for (Entry<Integer, HashSet<Integer>> e : this.b.getActiveVarIds().entrySet())
		{
			HashSet<Integer> hsRes = newActiveVarIdsMatrix.get(e.getKey());
			if (hsRes == null)
				newActiveVarIdsMatrix.put(e.getKey(), (HashSet<Integer>) e.getValue().clone());
			else
				hsRes.addAll(e.getValue());
		}
		return newActiveVarIdsMatrix;
	}

}
