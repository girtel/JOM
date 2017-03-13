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

import cern.colt.matrix.tdouble.DoubleMatrix2D;

import java.util.HashSet;
import java.util.LinkedHashMap;

class _FUNCTION_EYE extends Expression
{

	_FUNCTION_EYE(OptimizationProblem model, Expression a)
	{
		this(model, a, a);
	}

	_FUNCTION_EYE(OptimizationProblem model, Expression a, Expression b)
	{
		super(model);
		
		/* compute size */
		if (a.getNumDim() != 2)
			throw new JOMException("Function 'eye': Eye function should receive one or two constant numbers stating the size of the diagonal "
					+ "matrix");
		if (b.getNumDim() != 2)
			throw new JOMException("Function 'eye': Eye function should receive one or two constant numbers stating the size of the diagonal "
					+ "matrix");
		if ((a.getNumScalarExpressions() != 1) || (b.getNumScalarExpressions() != 1))
			throw new JOMException("Function 'eye': Eye function should receive one or two constant numbers stating the size of the diagonal "
					+ "matrix");
		if (!a.isConstant() || !b.isConstant())
			throw new JOMException("Function 'eye': Eye function should receive one or two constant numbers stating the size of the diagonal "
					+ "matrix");

		int aValue = (int) a.evaluateConstant(0);
		int bValue = (int) b.evaluateConstant(0);
		this.resize(new int[]{aValue, bValue});
		this.affineExp = a.getAffineExpression().function_eye(aValue, bValue);
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		throw new JOMException("Function 'eye': Unexpected error");
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		throw new JOMException("Function 'eye': Unexpected error");
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
		return true;
	}

	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		throw new JOMException("Function 'eye': Unexpected error");
	}

}
