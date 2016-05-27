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



 




package com.jom;

import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

class _OPERATOR_APPENDROWS extends Expression
{
	private final Expression a;
	private final Expression b;
	private final boolean isLinear;
	private final boolean isConstant;
	private final DoubleMatrix2D permutationMatrix;
	
	_OPERATOR_APPENDROWS(OptimizationProblem model, Expression a, Expression b)
	{
		super(model);

		/* compute size */
		if ((a.getNumDim() != 2) || (b.getNumDim() != 2)) throw new JOMException("Operator ';;' (append rows): This operator only applies to 2-D expressions");
		if (a.getSize() [1] != b.getSize() [1]) throw new JOMException("Operator ';;' (append rows): This operator only applies to 2-D expressions, with the same number of columns");
		int [] newSize = new int[] { a.getSize() [0] + b.getSize() [0] , a.getSize() [1] };
		this.resize(newSize);

		this.isLinear = ((a.isLinear()) && (b.isLinear()));
		this.isConstant = ((a.isConstant()) && (b.isConstant()));
		this.affineExp = (isLinear)? a.getAffineExpression().operator_appendRows(b.getAffineExpression()) : null;
		
		this.a = (this.isLinear())? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.b = (this.isLinear())? null : b; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.permutationMatrix = (this.isLinear())? null : computePermutationMatrixDouble(a, b); // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	private static DoubleMatrix2D computePermutationMatrixDouble(Expression a, Expression b)
	{
		int numRows_a = a.getSize() [0];
		int numRows_b = b.getSize() [0];
		int numCols = a.getSize() [1];
		int newNumRows = numRows_a + numRows_b;
		int newNumElem = newNumRows * numCols;
		int contCellNew = 0;
		DoubleMatrix2D res = DoubleFactory2D.sparse.make(newNumElem , newNumElem);
		for (int contCol = 0 ; contCol  < numCols ; contCol ++)
		{
			/* column of first expression */
			int positionInBadAppendedMatrix = 0 + contCol*numRows_a;
			res.viewPart(contCellNew, positionInBadAppendedMatrix, numRows_a, numRows_a).assign(DoubleFactory2D.sparse.identity(numRows_a));
			contCellNew += numRows_a;
			
			/* column of first expression */
			positionInBadAppendedMatrix = a.getNumScalarExpressions() + contCol*numRows_b;
			res.viewPart(contCellNew, positionInBadAppendedMatrix, numRows_b, numRows_b).assign(DoubleFactory2D.sparse.identity(numRows_b));
			contCellNew += numRows_b;
		}
		if (contCellNew != newNumElem) throw new JOMException("Operator ';;' (append rows): Unexpected error");
		return res;
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrixND aValue = a.evaluate_internal(valuesDVs);
		DoubleMatrixND bValue = b.evaluate_internal(valuesDVs);
		DoubleMatrixND res = new DoubleMatrixND(this.getSize(), this.permutationMatrix.zMult(DoubleFactory1D.sparse.append(aValue.elements(), bValue.elements()), null));
		return res;
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix2D aValue = a.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D bValue = b.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D res = this.permutationMatrix.zMult(DoubleFactory2D.sparse.appendRows(aValue, bValue) , null);
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
	LinkedHashMap<Integer,HashSet<Integer>> nl_getActiveVarIds()
	{
		LinkedHashMap<Integer,HashSet<Integer>> newActiveVarIdsMatrix = new LinkedHashMap<Integer,HashSet<Integer>> ();
		final int resNumRows = this.size [0];  
		final int numRows_a = a.getSize() [0];
		final int numRows_b = b.getSize() [0];
		for (Entry<Integer,HashSet<Integer>> cellEntry : this.a.getActiveVarIds().entrySet())
		{
			final int a_cellId = cellEntry.getKey();
			final HashSet<Integer> hs = cellEntry.getValue();
			final int a_Col = a_cellId / numRows_a;
			final int a_Row = a_cellId % numRows_a;
			newActiveVarIdsMatrix.put(a_Row + a_Col * resNumRows , (HashSet<Integer>) hs.clone());
		}
		for (Entry<Integer,HashSet<Integer>> cellEntry : this.b.getActiveVarIds().entrySet())
		{
			final int b_cellId = cellEntry.getKey();
			final HashSet<Integer> hs = cellEntry.getValue();
			final int b_Col = b_cellId / numRows_b;
			final int b_Row = b_cellId % numRows_b;
			newActiveVarIdsMatrix.put((b_Row + numRows_a) + b_Col * resNumRows , (HashSet<Integer>) hs.clone());
		}
		return newActiveVarIdsMatrix;
	}

	
}
