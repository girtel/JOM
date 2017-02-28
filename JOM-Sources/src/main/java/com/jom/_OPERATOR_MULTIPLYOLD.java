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

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.jet.math.tdouble.DoubleFunctions;

import java.util.*;
import java.util.Map.Entry;

@SuppressWarnings("unchecked")
class _OPERATOR_MULTIPLYOLD extends Expression
{
	private final Expression a;
	private final Expression b;
	private final boolean    isLinear;
	private final boolean    isConstant;

	_OPERATOR_MULTIPLYOLD(OptimizationProblem model, Expression a, Expression b)
	{
		super(model);

		/* compute size */
		if ((a.getNumDim() != 2) || (b.getNumDim() != 2)) throw new JOMException("Operator '*'  (matrix multiplication): Wrong size of the arrays");
		if ((a.getNumScalarExpressions() == 1) || (b.getNumScalarExpressions() == 1))
			throw new JOMException("Operator '*' (matrix multiplication): Unexpected error. For scalar multiplication use .* operator");
		if (a.getSize()[1] != b.getSize()[0])
			throw new JOMException("Operator '*'  (matrix multiplication): Wrong size of the arrays for multiplication. Size left matrix: " + Arrays
					.toString(a.getSize()) + ". Size right matrix: " + Arrays.toString(b.getSize()));
		this.resize(new int[]{a.getSize()[0], b.getSize()[1]});

		this.isLinear = ((a.isConstant() && b.isLinear()) || (b.isConstant() && a.isLinear()));
		this.affineExp = (isLinear) ? _INTERNAL_AffineExpressionCoefs.operator_multiply(a.getAffineExpression(), b.getAffineExpression()) : null;
		this.isConstant = (isLinear) ? this.affineExp.isConstant() : false;

		this.a = (isLinear) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.b = (isLinear) ? null : b; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrixND a_value = a.evaluate_internal(valuesDVs);
		DoubleMatrixND b_value = b.evaluate_internal(valuesDVs);
		return new DoubleMatrixND(this.getSize(), a_value.view2D().zMult(b_value.view2D(), null).vectorize());
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrixND a_value = a.evaluate_internal(valuesDVs);
		DoubleMatrixND b_value = b.evaluate_internal(valuesDVs);
		DoubleMatrix2D a_gradient = a.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D b_gradient = b.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D res = DoubleFactory2D.dense.make(this.getSize()[0] * this.getSize()[1], this.getModel().getNumScalarDecisionVariables());
		for (int newRow = 0; newRow < this.getSize()[0]; newRow++)
		{
			for (int newCol = 0; newCol < this.getSize()[1]; newCol++)
			{
				int cellIndex_result = DoubleMatrixND.sub2ind(IntFactory1D.dense.make(new int[]{newRow, newCol}), this.getSize());
				DoubleMatrix1D gradient = DoubleFactory1D.dense.make(this.getModel().getNumScalarDecisionVariables());
				for (int contSummand = 0; contSummand < a.getSize()[1]; contSummand++)
				{
					int cellIndex_a = DoubleMatrixND.sub2ind(IntFactory1D.dense.make(new int[]{newRow, contSummand}), a.getSize());
					int cellIndex_b = DoubleMatrixND.sub2ind(IntFactory1D.dense.make(new int[]{contSummand, newCol}), b.getSize());
					DoubleMatrix1D aux = a_gradient.viewRow(cellIndex_a).copy().assign(DoubleFunctions.mult(b_value.get(cellIndex_b)));
					gradient.assign(aux, DoubleFunctions.plus);
					aux = b_gradient.viewRow(cellIndex_b).copy().assign(DoubleFunctions.mult(a_value.get(cellIndex_a)));
					gradient.assign(aux, DoubleFunctions.plus);
				}
				res.viewRow(cellIndex_result).assign(gradient);
			}
		}
		return res;
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
		final int aRows = a.size[0];
		final int aCols = a.size[1];
		final int bRows = b.size[0];
		final int bCols = b.size[1];
		final int resRows = aRows;
		final int resCols = bCols;
		final int resNumCells = resRows * resCols;
		final int[] resSize = new int[]{aRows, bCols};
		LinkedHashMap<Integer, HashSet<Integer>> newActiveVarIdsMatrix = new LinkedHashMap<Integer, HashSet<Integer>>();
		int previousColumn_b = 0;
		ArrayList<HashSet<Integer>> aColumn = new ArrayList<HashSet<Integer>>(bRows);
		Iterator<Entry<Integer, HashSet<Integer>>> b_it = b.getActiveVarIds().entrySet().iterator();
		while (b_it.hasNext())
		{
			Entry<Integer, HashSet<Integer>> bCellEntry = b_it.next();
			final int bCellId = bCellEntry.getKey();
			final HashSet<Integer> bHs = bCellEntry.getValue();
			final int bCol = bCellId / bRows;
			final int bRow = bCellId % bRows;
			if ((!b_it.hasNext()) || (bCol != previousColumn_b))
			{
				/* the last time we execute this (possibly in the last column of b)  */
				if (!b_it.hasNext()) aColumn.set(bRow, bHs);

				/* multiply LM with the column of RM. Accumulate the result in res_linearCoefs */
				for (Entry<Integer, HashSet<Integer>> aCellEntry : a.getActiveVarIds().entrySet())
				{
					final int aCellId = aCellEntry.getKey();
					final HashSet<Integer> aHs = aCellEntry.getValue();
					final int aRow = aCellId % aRows;
					final int aCol = aCellId / aRows;
					final int resRow = aRow;
					final int resCol = previousColumn_b;
					final int resCellId = resRow + resRows * resCol;
					final HashSet<Integer> aCellColumnVector = aColumn.get(aCol);
					HashSet<Integer> summand = ((HashSet<Integer>) aHs.clone());
					summand.addAll(aCellColumnVector);
					HashSet<Integer> accumHs = newActiveVarIdsMatrix.get(resCellId);
					if (accumHs == null)
						newActiveVarIdsMatrix.put(resCellId, summand);
					else
						accumHs.addAll(summand);
				}
				previousColumn_b = bCol;
				aColumn = new ArrayList<HashSet<Integer>>(bRows);
			}
			aColumn.set(bRow, bHs);
		}

		return newActiveVarIdsMatrix;
	}

}
