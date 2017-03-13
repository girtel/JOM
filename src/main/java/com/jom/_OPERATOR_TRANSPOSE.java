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

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;
import cern.colt.matrix.tint.IntMatrix2D;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

class _OPERATOR_TRANSPOSE extends Expression
{
	private final Expression     a;
	private final DoubleMatrix2D permutationMatrix;
	private final boolean        isLinear;
	private final boolean        isConstant;

	_OPERATOR_TRANSPOSE(OptimizationProblem model, Expression a)
	{
		super(model);

		/* compute size */
		if (a.getNumDim() != 2) throw new JOMException("Operator ' (transpose). Transpose operation is only valid in 2D arrays");
		this.resize(new int[]{a.getSize()[1], a.getSize()[0]});

		this.isLinear = a.isLinear();
		this.isConstant = a.isConstant();
		this.affineExp = (isLinear) ? a.getAffineExpression().operator_transpose() : null;

		this.a = (this.isLinear()) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.permutationMatrix = (this.isLinear()) ? null : computePermutationMatrix(a.getSize());
		; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	private static DoubleMatrix2D computePermutationMatrix(int[] oldSize)
	{
		IntMatrix1D newSize = IntFactory1D.dense.make(new int[]{oldSize[1], oldSize[0]});
		int numScalarExpressions = oldSize[0] * oldSize[1];
		IntMatrix1D indexesOld = new IntMatrixND(new int[]{1, numScalarExpressions}).ascending(0, 1).elements();
		IntMatrix2D subindexesOld = IntMatrixND.ind2sub(indexesOld, IntFactory1D.dense.make(oldSize));
		IntMatrix2D subindexesNew = subindexesOld.copy();
		subindexesNew.viewColumn(0).assign(subindexesOld.viewColumn(1));
		subindexesNew.viewColumn(1).assign(subindexesOld.viewColumn(0));
		IntMatrix1D indexesNew = IntMatrixND.sub2ind(subindexesNew, newSize);
		DoubleMatrix2D mat = DoubleFactory2D.sparse.make(numScalarExpressions, numScalarExpressions, 0.0);
		for (int cont = 0; cont < numScalarExpressions; cont++)
			mat.set(indexesNew.get(cont), indexesOld.get(cont), 1.0);
		return mat;
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		return new DoubleMatrixND(a.evaluate_internal(valuesDVs).view2D().viewDice());
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		return this.permutationMatrix.zMult(this.a.evaluateJacobian_internal(valuesDVs), null);
	}

	@Override
	boolean isLinear()
	{
		return this.isLinear;
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
		/* Make a copy since I am modifying it */
		LinkedHashMap<Integer, HashSet<Integer>> newActiveVarIdsMatrix = new LinkedHashMap<Integer, HashSet<Integer>>();
		final int oldRows = this.size[1];
		final int oldCols = this.size[0];
		for (Entry<Integer, HashSet<Integer>> e : a.getActiveVarIds().entrySet())
		{
			final int cellId = e.getKey();
			final HashSet<Integer> hs = e.getValue();
			final int oldCol = cellId / oldRows;
			final int oldRow = cellId % oldRows;
			newActiveVarIdsMatrix.put(oldCol + oldCols * oldRow, hs);
		}
		return newActiveVarIdsMatrix;
	}

}
