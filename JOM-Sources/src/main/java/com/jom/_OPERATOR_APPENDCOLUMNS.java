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
import cern.colt.matrix.tdouble.DoubleMatrix2D;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

@SuppressWarnings("unchecked")
class _OPERATOR_APPENDCOLUMNS extends Expression
{

	private final Expression a;
	private final Expression b;
	private final boolean    isLinear;
	private final boolean    isConstant;

	_OPERATOR_APPENDCOLUMNS(OptimizationProblem model, Expression a, Expression b)
	{
		super(model);

		/* compute size */
		if ((a.getNumDim() != 2) || (b.getNumDim() != 2))
			throw new JOMException("Operator ';' (append columns): This operator only applies to 2-D expressions");
		if (a.getSize()[0] != b.getSize()[0])
			throw new JOMException("Operator ';' (append columns): This operator only applies to 2-D expressions, with the same number of rows");
		int[] newSize = new int[]{a.getSize()[0], a.getSize()[1] + b.getSize()[1]};
		this.resize(newSize);

		this.isLinear = ((a.isLinear()) && (b.isLinear()));
		this.isConstant = ((a.isConstant()) && (b.isConstant()));
		this.affineExp = (isLinear) ? a.getAffineExpression().operator_appendColumns(b.getAffineExpression()) : null;

		this.a = (this.isLinear()) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.b = (this.isLinear()) ? null : b; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrixND aValue = a.evaluate_internal(valuesDVs);
		DoubleMatrixND bValue = b.evaluate_internal(valuesDVs);
		DoubleMatrixND res = new DoubleMatrixND(this.getSize(), DoubleFactory1D.sparse.append(aValue.elements(), bValue.elements()));
		return res;
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix2D aValue = a.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D bValue = b.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D res = DoubleFactory2D.sparse.appendRows(aValue, bValue);
		return res;
	}

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
		return this.isConstant;
	}

	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		/* Make a copy since I am modifying it */
		final int aNumCells = a.getNumScalarExpressions();
		LinkedHashMap<Integer, HashSet<Integer>> newActiveVarIdsMatrix = new LinkedHashMap<Integer, HashSet<Integer>>();
		for (Entry<Integer, HashSet<Integer>> aCell : a.getActiveVarIds().entrySet())
			newActiveVarIdsMatrix.put(aCell.getKey(), (HashSet<Integer>) aCell.getValue().clone());
		for (Entry<Integer, HashSet<Integer>> bCell : b.getActiveVarIds().entrySet())
			newActiveVarIdsMatrix.put(bCell.getKey() + aNumCells, (HashSet<Integer>) bCell.getValue().clone());
		return newActiveVarIdsMatrix;
	}

}
