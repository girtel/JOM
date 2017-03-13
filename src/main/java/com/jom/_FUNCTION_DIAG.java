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

import java.util.HashSet;
import java.util.LinkedHashMap;

@SuppressWarnings("unchecked")
class _FUNCTION_DIAG extends Expression
{
	private final Expression a;
	private final int        numElemDiag;
	private final boolean    isLinear;
	private final boolean    isConstant;

	_FUNCTION_DIAG(OptimizationProblem model, Expression a)
	{
		super(model);
		
		/* compute size */
		if (a.getNumDim() != 2) throw new JOMException("Function 'diag': Diag function is for 2-D expressions in the form of row or column vector");
		if ((a.getSize()[0] != 1) && (a.getSize()[1] != 1))
			throw new JOMException("Function 'diag': Diag function is for 2-D expressions in the form of row or column vector");

		this.numElemDiag = a.getNumScalarExpressions();
		this.resize(new int[]{numElemDiag, numElemDiag});

		this.isLinear = (a.isLinear());
		this.isConstant = a.isConstant();
		this.affineExp = (isLinear) ? a.getAffineExpression().function_diag() : null;

		this.a = (isLinear) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrix1D aValues = a.evaluate_internal(valuesDVs).elements();
		return new DoubleMatrixND(DoubleFactory2D.sparse.diagonal(aValues));
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix2D aGradient = a.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D res = DoubleFactory2D.dense.make(this.numScalarExpressions, aGradient.columns());
		for (int row = 0; row < res.rows(); row++)
			res.viewRow(row * (this.numElemDiag + 1)).assign(aGradient.viewRow(row).copy());
		return res;
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
		return this.isConstant;
	}

	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		/* Make a copy since I am modifying it */
		LinkedHashMap<Integer, HashSet<Integer>> activeVarIdsMatrix = new LinkedHashMap<Integer, HashSet<Integer>>();
		HashSet<Integer> copyOfOriginalActiveVarIds = a.getActiveVarIds().get(0);
		for (int cont = 0; cont < this.size[0]; cont++)
			activeVarIdsMatrix.put(cont * (size[0] + 1), (HashSet<Integer>) copyOfOriginalActiveVarIds.clone());
		return activeVarIdsMatrix;
	}

}
