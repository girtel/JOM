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

import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntMatrix1D;
import cern.jet.math.tdouble.DoubleFunctions;

class _OPERATOR_MULTIPLY extends Expression
{
	private final Expression a;
	private final Expression b;
	private final boolean    isLinear;
	private final boolean    isConstant;

	_OPERATOR_MULTIPLY(OptimizationProblem model, Expression a, Expression b)
	{
		super(model);

		/* compute size */
		if ((a.getNumDim() != 2) || (b.getNumDim() < 2))
			throw new JOMException("Operator '*'  (matrix multiplication): Wrong size of the arrays. First array must be 2D, second array can have 2"
					+ " or more dimensions");
		if ((a.getNumScalarExpressions() == 1) || (b.getNumScalarExpressions() == 1))
			throw new JOMException("Operator '*' (matrix multiplication): Unexpected error. For scalar multiplication use .* operator");
		if (a.getSize()[1] != b.getSize()[0])
			throw new JOMException("Operator '*'  (matrix multiplication): Wrong size of the arrays for multiplication. Size left matrix: " + Arrays
					.toString(a.getSize()) + ". Size right matrix: " + Arrays.toString(b.getSize()));

		final int[] newSize = Arrays.copyOf(b.size, b.size.length);
		newSize[0] = a.size[0];
		this.resize(newSize);

		this.isLinear = ((a.isConstant() && b.isLinear()) || (b.isConstant() && a.isLinear()));
		this.affineExp = (isLinear) ? _INTERNAL_AffineExpressionCoefs.operator_multiply(a.getAffineExpression(), b.getAffineExpression()) : null;
		this.isConstant = (isLinear) ? this.affineExp.isConstant() : false;

		this.a = (isLinear) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.b = (isLinear) ? null : b; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrix2D a_value = a.evaluate_internal(valuesDVs).view2D();
		DoubleMatrix1D b_value_elements = b.evaluate_internal(valuesDVs).elements();
		final int numberMatrixMultiplications = this.numScalarExpressions / (size[0] * size[1]);
		final int sizeEachBSubMatrix = b.size[0] * b.size[1];
		final int sizeEachResSubMatrix = a.size[0] * b.size[1];
		DoubleMatrix1D res_values = DoubleFactory1D.sparse.make(this.numScalarExpressions);
		for (int cont = 0; cont < numberMatrixMultiplications; cont++)
		{
			DoubleMatrix2D bMatrix = b_value_elements.viewPart(cont * sizeEachBSubMatrix, sizeEachBSubMatrix).reshape(b.size[0], b.size[1]);
			res_values.viewPart(cont * sizeEachResSubMatrix, sizeEachResSubMatrix).assign(a_value.zMult(bMatrix, null).vectorize());
		}
		return new DoubleMatrixND(this.getSize(), res_values);
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrixND a_value = a.evaluate_internal(valuesDVs);
		DoubleMatrixND b_value = b.evaluate_internal(valuesDVs);
		DoubleMatrix2D a_gradient = a.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D b_gradient = b.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D res = DoubleFactory2D.dense.make(this.numScalarExpressions, this.getModel().getNumScalarDecisionVariables());
		for (int newIndex = 0; newIndex < this.numScalarExpressions; newIndex++)
		{
			final IntMatrix1D subindexesInRes = DoubleMatrixND.ind2sub(newIndex, this.size);
			final int newRow = subindexesInRes.get(0);
			if (newRow != newIndex % a.size[0]) throw new RuntimeException("Bad");
			final int newCol = subindexesInRes.get(1);
			IntMatrix1D subIndexesInB = subindexesInRes.copy();
			subIndexesInB.set(1, newCol);
			DoubleMatrix1D gradient = DoubleFactory1D.dense.make(this.getModel().getNumScalarDecisionVariables());
			for (int contSummand = 0; contSummand < a.getSize()[1]; contSummand++)
			{
				final int cellIndex_a = newRow + contSummand * a.size[0];
				subIndexesInB.set(0, contSummand);
				final int cellIndex_b = DoubleMatrixND.sub2ind(subIndexesInB, b.getSize());
				DoubleMatrix1D aux = a_gradient.viewRow(cellIndex_a).copy().assign(DoubleFunctions.mult(b_value.get(cellIndex_b)));
				gradient.assign(aux, DoubleFunctions.plus);
				aux = b_gradient.viewRow(cellIndex_b).copy().assign(DoubleFunctions.mult(a_value.get(cellIndex_a)));
				gradient.assign(aux, DoubleFunctions.plus);
			}
			res.viewRow(newIndex).assign(gradient);
		}
		return res;
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
		HashMap<Integer,HashSet<Integer>> aColumn = new HashMap<Integer,HashSet<Integer>>();
		Iterator<Entry<Integer, HashSet<Integer>>> b_it = b.getActiveVarIds().entrySet().iterator();
		while (b_it.hasNext())
		{
			Entry<Integer, HashSet<Integer>> bCellEntry = b_it.next();
			final int bCellId = bCellEntry.getKey();
			final HashSet<Integer> bHs = bCellEntry.getValue();
			final IntMatrix1D subindexes_b = DoubleMatrixND.ind2sub(bCellId, b.size);
			final int bCol = subindexes_b.get(1);//bCellId / bRows;
			final int bRow = subindexes_b.get(0);//bCellId % bRows;
			if ((!b_it.hasNext()) || (bCol != previousColumn_b))
			{
				/* the last time we execute this (possibly in the last column of b)  */
				if (!b_it.hasNext()) aColumn.put(bRow, bHs);

				IntMatrix1D subindexes_res = subindexes_b.copy(); 
				/* multiply LM with the column of RM. Accumulate the result in res_linearCoefs */
				for (Entry<Integer, HashSet<Integer>> aCellEntry : a.getActiveVarIds().entrySet())
				{
					final int aCellId = aCellEntry.getKey();
					final HashSet<Integer> aHs = aCellEntry.getValue();
					final int aRow = aCellId % aRows;
					final int aCol = aCellId / aRows;
					final int resRow = aRow;
					final int resCol = previousColumn_b;
					subindexes_res.set(0, resRow);
					subindexes_res.set(1, resCol);
					final int resCellId = DoubleMatrixND.sub2ind(subindexes_res, this.size);//resRow + resRows * resCol;
					final HashSet<Integer> aCellColumnVector = aColumn.get(aCol);
					HashSet<Integer> summand = ((HashSet<Integer>) aHs.clone());
					if (aCellColumnVector != null) summand.addAll(aCellColumnVector);
					HashSet<Integer> accumHs = newActiveVarIdsMatrix.get(resCellId);
					if (accumHs == null)
						newActiveVarIdsMatrix.put(resCellId, summand);
					else
						accumHs.addAll(summand);
				}
				previousColumn_b = bCol;
				aColumn = new HashMap<Integer,HashSet<Integer>>();
			}
			aColumn.put(bRow, bHs);
		}

		return newActiveVarIdsMatrix;
	}

}
