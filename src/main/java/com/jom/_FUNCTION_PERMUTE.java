/*******************************************************************************
 * Copyright (c) 2015 Pablo Pavon Mariño.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * 
 * Contributors:
 *     Pablo Pavon Mariño - initial API and implementation
 ******************************************************************************/



 




/**
 * 
 */
package com.jom;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;
import java.util.Set;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;

class _FUNCTION_PERMUTE extends Expression
{
	private final Expression a;
	private final int [] permutationOfDimensionsIndex0;
	private final DoubleMatrix2D permutationMatrix;
	private final boolean isLinear;
	private final boolean isConstant;
	
	_FUNCTION_PERMUTE(OptimizationProblem model, Expression a , Expression dimPermutation)
	{
		super(model);

		if (!dimPermutation.isConstant()) throw new JOMException ("Function 'permute': permutation argument must be a constant");
		if (!dimPermutation.isRowOrColumn2DExpression()) throw new JOMException ("Function 'permute': the permutation must be a row or column vector");
		if (dimPermutation.getNumScalarExpressions() != a.getNumDim()) throw new JOMException ("Function 'permute': the number of elements in the permutation argument must be equal to the number of dimensions");
		this.permutationOfDimensionsIndex0 = new int [a.getNumDim()];
		final int [] newSize = new int [a.size.length];
		Set<Integer> dimensionsUsed = new HashSet<Integer> ();
		for (int oldDimensionIndex0 = 0 ; oldDimensionIndex0 < a.numDim ; oldDimensionIndex0 ++)
		{ 
			final int newDimensionIndex0 = ((int) dimPermutation.evaluateConstant(oldDimensionIndex0)) - 1;
			if ((newDimensionIndex0 < 0) || (newDimensionIndex0 >= a.numDim)) throw new JOMException ("Function 'permute': wrong dimension identifier");
			permutationOfDimensionsIndex0 [oldDimensionIndex0] = newDimensionIndex0;
			newSize [oldDimensionIndex0] = a.size [newDimensionIndex0];
			dimensionsUsed.add (newDimensionIndex0);
		}
		if (dimensionsUsed.size() != a.numDim) throw new JOMException ("Function 'permute': dimensions cannot be repeated");
		
			/* compute size */
		this.resize(newSize);
		this.isLinear = a.isLinear();
		this.isConstant = a.isConstant();
		this.affineExp = (isLinear)? a.getAffineExpression().function_permute (permutationOfDimensionsIndex0) : null;

		this.a = (this.isLinear())? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.permutationMatrix = (this.isLinear())? null : computePermutationMatrix (a.getSize() , newSize , permutationOfDimensionsIndex0); // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	private static DoubleMatrix2D computePermutationMatrix (int [] oldSize , int [] newSize , int [] permutationOfDimensionsIndex0)
	{
		final int numScalarExpressions = IntMatrixND.prod(newSize);
		DoubleMatrix2D mat = DoubleFactory2D.sparse.make(numScalarExpressions , numScalarExpressions , 0.0);
		for (int oldCellIndex = 0 ; oldCellIndex < numScalarExpressions ; oldCellIndex ++)
		{
			final IntMatrix1D oldSubindexes = DoubleMatrixND.ind2sub(oldCellIndex, oldSize);
			final IntMatrix1D newSubindexes = IntFactory1D.dense.make(newSize.length); for (int dim = 0; dim < newSize.length ; dim ++) newSubindexes.set(dim , oldSubindexes.get(permutationOfDimensionsIndex0 [dim]));
			final int newIndex = DoubleMatrixND.sub2ind(newSubindexes, newSize);
			mat.set(newIndex , oldCellIndex , 1.0);
		}
		return mat;
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		return new DoubleMatrixND (a.evaluate_internal(valuesDVs).view2D().viewDice());
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		return this.permutationMatrix.zMult (this.a.evaluateJacobian_internal(valuesDVs) , null);
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
	LinkedHashMap<Integer,HashSet<Integer>> nl_getActiveVarIds()
	{
		/* Make a copy since I am modifying it */
		LinkedHashMap<Integer,HashSet<Integer>> newActiveVarIdsMatrix = new LinkedHashMap<Integer,HashSet<Integer>> ();
		for (Entry<Integer,HashSet<Integer>> e : a.getActiveVarIds().entrySet())
		{
			final int cellId = e.getKey();
			final HashSet<Integer> hs = e.getValue();
			final IntMatrix1D oldSubindexes = DoubleMatrixND.ind2sub(cellId, a.size);
			final IntMatrix1D newSubindexes = IntFactory1D.dense.make(this.numDim); for (int dim = 0; dim < this.size.length ; dim ++) newSubindexes.set(dim , oldSubindexes.get(this.permutationOfDimensionsIndex0 [dim]));
			final int newIndex = DoubleMatrixND.sub2ind(newSubindexes, this.size);
			newActiveVarIdsMatrix.put(newIndex, hs);
		}
		return newActiveVarIdsMatrix;
	}


	
}
