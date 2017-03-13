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

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;

import java.util.HashSet;
import java.util.LinkedHashMap;

class _FUNCTION_ATAN extends Expression
{
	private final Expression a;
	private final boolean    isConstant;
	private final boolean    isLinear;

	_FUNCTION_ATAN(OptimizationProblem model, Expression a)
	{
		super(model);

		this.resize(a.getSize());

		this.isLinear = a.isConstant();
		this.isConstant = isLinear;
		if (isLinear)
		{
			/* change (destroy) the constant data in Expression a. This object keeps no reference to Expression a, but the data remains since we 
			 * keep reference to the affine expression object of a, now modified by this function */
			double[] data = a.getAffineExpression().getConstantCoefArray();
			for (int cont = 0; cont < data.length; cont++) data[cont] = Math.atan(data[cont]);
			this.affineExp = new _INTERNAL_AffineExpressionCoefs(model, a.getSize(), data); //a.getAffineExpression();
		}
		this.a = (this.isLinear()) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrix1D aValue = a.evaluate_internal(valuesDVs).elements();
		return new DoubleMatrixND(this.getSize(), aValue.assign(DoubleFunctions.atan));
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix1D a_value = a.evaluate_internal(valuesDVs).elements();
		DoubleMatrix2D a_jacobian = a.evaluateJacobian_internal(valuesDVs);
		return DoubleFactory2D.sparse.diagonal(a_value.assign(DoubleFunctions.square).assign(DoubleFunctions.plus(1)).assign(DoubleFunctions.inv))
				.zMult(a_jacobian, null);
	}

	/* (non-Javadoc)
	 * @see Expression#isLinear()
	 */
	@Override
	boolean isLinear()
	{
		return isLinear;
	}

	/* (non-Javadoc)
	 * @see Expression#isConstant()
	 */
	@Override
	boolean isConstant()
	{
		return isConstant;
	}

	/* (non-Javadoc)
	 * @see Expression#getActiveVarIds()
	 */
	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		return a.getActiveVarIds();
	}

}
