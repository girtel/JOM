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

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

class _FUNCTION_SUM extends Expression
{
	private final boolean isLinear;
	private final boolean isConstant;

	private final Expression a;
	private final int dimension;
	private       DoubleMatrix2D sumOperationMatrix; // rows = num cells this (new) , columns = num cells a (old). 1 if new cell includes old cell
	// as summand

	_FUNCTION_SUM(OptimizationProblem model, Expression a)
	{
		this(model, a, -1);
	}

	_FUNCTION_SUM(OptimizationProblem model, Expression a, Expression dim)
	{
		/* For the user, first dimension starts in 1. In the whole JOM, first dimension is 0. This is changed in just next line !!! */
		//this (model , a , (dim == null)? computeFirstNonSingletonDime_0Index (a)  :  (int) dim.evaluate().get(0) ); 
		this(model, a, ((int) dim.evaluateConstant().get(0) - 1));
		if (!dim.isScalarConstant()) throw new JOMException("Function 'sum': Dimension argument must be a scalar constant");
	}

	private _FUNCTION_SUM(OptimizationProblem model, Expression a, int dimension)
	{
		super(model);
		
		/* compute size */
		this.dimension = dimension;
		/* If sum all the elements of the array => we produce a 1x1 matrix */
		if (dimension >= a.getNumDim()) throw new JOMException("Function 'sum': Dimension argument out of range");
		if (dimension == -1) this.resize(new int[]{1, 1});
		else if (a.getNumDim() == 2)
		{
			int[] newSize = Arrays.copyOf(a.getSize(), a.getSize().length);
			newSize[dimension] = 1;
			this.resize(newSize);
		} else
		{
			/* reduce the number of dimensions */
			int[] newSize = new int[a.getNumDim() - 1];
			for (int dim = 0; dim < dimension; dim++) newSize[dim] = a.getSize()[dim];
			for (int dim = dimension; dim < newSize.length; dim++) newSize[dim] = a.getSize()[dim + 1];
			this.resize(newSize);
		}

		this.isLinear = a.isLinear();
		this.affineExp = (isLinear) ? a.getAffineExpression().function_sum(dimension, this.size) : null;
		this.isConstant = (isLinear) ? this.affineExp.isConstant() : a.isConstant();
		
		/* The sum operation matrix precomputed, only for non-lienar expressions */
		this.sumOperationMatrix = (this.isLinear()) ? null : computeSumOperationMatrix(a, dimension);
		this.a = (this.isLinear()) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	private static DoubleMatrix2D computeSumOperationMatrix(Expression a, int dimension)
	{
		/* You sum all */
		if (dimension == -1)
			return DoubleFactory2D.dense.make(1, a.getNumScalarExpressions(), 1.0);

		/* You sum a dimension */
		int numCellsPreviousDimensions = (dimension == 0) ? 0 : IntMatrixND.prod(Arrays.copyOf(a.getSize(), dimension));
		int numCellsSum = a.getSize()[dimension];
		int numCellsPosteriorDimensions = (dimension == a.getNumDim() - 1) ?
				0 :
				IntMatrixND.prod(a.getSize(), dimension + 1, a.getNumDim() - dimension - 1);

		/* You sum nothing */
		if (numCellsSum < 1) throw new JOMException("Function 'sum': Unexpected error");
		if (numCellsSum == 1)
			return DoubleFactory2D.sparse.identity(a.getNumScalarExpressions());

		DoubleMatrix2D firstMatrix = DoubleFactory2D.sparse.identity(numCellsPreviousDimensions);
		DoubleMatrix2D secondMatrix = DoubleFactory2D.sparse.make(1, numCellsSum, 1.0);
		DoubleMatrix2D thirdMatrix = DoubleFactory2D.sparse.identity(numCellsPosteriorDimensions);

		DoubleMatrix2D auxSumOperationMatrix;
		if (numCellsPreviousDimensions > 0)
			auxSumOperationMatrix = DoubleMatrixND.kron(secondMatrix, firstMatrix);
		else
			auxSumOperationMatrix = secondMatrix;

		if (numCellsPosteriorDimensions > 0)
			return DoubleMatrixND.kron(thirdMatrix, auxSumOperationMatrix);
		else
			return auxSumOperationMatrix;
	}

	private static int[] computeBlockSizes(int[] size)
	{
		int[] blockSizes = new int[size.length];
		blockSizes[0] = 1;
		for (int dim = 1; dim < size.length; dim++)
			blockSizes[dim] = blockSizes[dim - 1] * size[dim - 1];
		return blockSizes;
	}

	private static int[] ind2sub(int index, int[] size, int[] blockSizes)
	{
		final int numDim = size.length;
		int[] res = new int[size.length];
		int indexToReduce = index;
		for (int contDim = numDim - 1; contDim >= 0; contDim--)
		{
			res[contDim] = indexToReduce / blockSizes[contDim];
			indexToReduce = indexToReduce % blockSizes[contDim];
		}
		return res;
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrix1D a_value = a.evaluate_internal(valuesDVs).elements();
		DoubleMatrix1D res = this.sumOperationMatrix.zMult(a_value, null);
		return new DoubleMatrixND(this.size(), res);
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix2D a_value = a.evaluateJacobian_internal(valuesDVs);
		return this.sumOperationMatrix.zMult(a_value, null);
	}

	@Override
	boolean isLinear()
	{
		return this.isLinear;
	}

	@Override
	boolean isConstant()
	{
		return this.isConstant;
	}

	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		LinkedHashMap<Integer, HashSet<Integer>> newActiveVarIds = new LinkedHashMap<Integer, HashSet<Integer>>();
		LinkedHashMap<Integer, HashSet<Integer>> a_activeVarIds = a.getActiveVarIds();

			/* sum all dimensions producing a scalar  */
		if (dimension == -1)
		{
			HashSet<Integer> newHs = new HashSet<Integer>();
			for (HashSet<Integer> hs : a_activeVarIds.values())
				for (Integer dv : hs)
					newHs.add(dv);
			newActiveVarIds.put(0, newHs);
			return newActiveVarIds;
		}

			/* Sum in a dimension */
		final int[] oldSize = a.getSize();
		final int[] newSize = this.size;
		final int[] blocksOldSize = computeBlockSizes(oldSize);
		final int[] blocksNewSize = computeBlockSizes(newSize);
		for (Entry<Integer, HashSet<Integer>> oldCellEntry : a.getActiveVarIds().entrySet())
		{
			final int oldCellId = oldCellEntry.getKey();
			final HashSet<Integer> varIdsToSum = oldCellEntry.getValue();
			final int[] oldCoordinates = ind2sub(oldCellId, oldSize, blocksOldSize);
			/* compute the index of the new cell */
			int newCellIndex = 0;
			for (int cont = 0; cont <= dimension - 1; cont++) newCellIndex += blocksNewSize[cont] + oldCoordinates[cont];
			for (int cont = dimension + 1; cont < oldSize.length; cont++) newCellIndex += blocksNewSize[cont - 1] * oldCoordinates[cont];
			/* Accumulate in that cell */
			HashSet<Integer> accumCell = newActiveVarIds.get(newCellIndex);
			if (accumCell == null)
				newActiveVarIds.put(newCellIndex, varIdsToSum);
			else
				accumCell.addAll(varIdsToSum);
		}
		return newActiveVarIds;
	}

}
