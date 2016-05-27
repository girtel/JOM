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

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;

class _FUNCTION_REPEATEXPRESSION extends Expression
{
	private final Expression a;
	private final boolean isLinear;
	private final boolean isConstant;
	
	_FUNCTION_REPEATEXPRESSION(OptimizationProblem model, Expression a, int [] size)
	{
		super(model);
		
		/* compute size */
		if (a.getNumScalarExpressions() != 1) throw new JOMException ("Function 'repeat expression': Repeat function is for scalar expressions");
		this.resize(size);

		this.isLinear = (a.isLinear()); 
		this.isConstant = a.isConstant();
		this.affineExp = (isLinear)? a.getAffineExpression().function_repeatExpression(size) : null;
		
		this.a = (isLinear)? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		double aValue = a.evaluate_internal(valuesDVs).elements().get(0);
		return new DoubleMatrixND(this.getSize() , aValue);
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix2D aGradient = a.evaluateJacobian_internal(valuesDVs);
		return DoubleFactory2D.sparse.repeat(aGradient, this.getNumScalarExpressions() , 1);
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
		return this.isConstant;
	}

	@Override
	LinkedHashMap<Integer,HashSet<Integer>> nl_getActiveVarIds()
	{
		/* Make a copy since I am modifying it */
		LinkedHashMap<Integer,HashSet<Integer>> activeVarIdsMatrix = new LinkedHashMap<Integer,HashSet<Integer>> ();
		HashSet<Integer> copyOfOriginalActiveVarIds = (HashSet<Integer>) a.getActiveVarIds().get(0);
		for (int cont = 0 ;  cont < this.numScalarExpressions ; cont ++) activeVarIdsMatrix.put (cont , (HashSet <Integer>) copyOfOriginalActiveVarIds.clone());
		return activeVarIdsMatrix;
	}

}
