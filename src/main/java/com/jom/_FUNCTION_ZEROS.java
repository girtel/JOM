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

class _FUNCTION_ZEROS extends Expression
{

	_FUNCTION_ZEROS(OptimizationProblem model, Expression a)
	{
		super(model);
		
		/* compute size */
		if (a.getNumDim() != 2)
			throw new JOMException("Function 'zeros': Zeros function should receive a row or column vector of constants indicanting the array size");
		if ((a.getSize()[0] != 1) && (a.getSize()[1] != 1))
			throw new JOMException("Function 'zeros': Zeros function should receive a row or column vector of constants indicanting the array size");
		if (!a.isConstant())
			throw new JOMException("Function 'zeros': Zeros function should receive a row or column vector of constants indicanting the array size");

		double[] size_d = a.evaluateConstant().to1DArray();
		int[] size_i = new int[size_d.length];
		for (int cont = 0; cont < size_i.length; cont++) size_i[cont] = (int) size_d[cont];

		this.resize(size_i);
		this.affineExp = a.getAffineExpression().function_onesOrZeros(0);
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		throw new JOMException("Function 'zeros': Unexpected error");
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		throw new JOMException("Function 'zeros': Unexpected error");
	}

	/* (non-Javadoc)
	 * @see com.jom.Expression#isLinear()
	 */
	@Override
	boolean isLinear()
	{
		return true;
	}

	/* (non-Javadoc)
	 * @see com.jom.Expression#isConstant()
	 */
	@Override
	boolean isConstant()
	{
		return true;
	}

	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		throw new JOMException("Function 'zeros': Unexpected error");
	}

}
