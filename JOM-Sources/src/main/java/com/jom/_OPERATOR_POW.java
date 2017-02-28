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
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.jet.math.tdouble.DoubleFunctions;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

class _OPERATOR_POW extends Expression
{
	private final Expression     a;
	private final double         aValue;
	private final DoubleMatrix1D aValues;
	private final Expression     b;
	private final double         bValue;
	private final DoubleMatrix1D bValues;
	private final boolean        constantAllBases;
	private final boolean        constantAllExponents;
	private final boolean        constantBase;
	private final boolean        constantExponent;

	private final boolean isConstant;
	private final boolean isLinear;

	//private final double bValue;

	_OPERATOR_POW(OptimizationProblem model, Expression a, Expression b)
	{
		super(model); //, _OPERATOR_POW.computeSize(a, b), computeExpressionString(a, b), computeCurvature(a, b), computeAffineExprObject(a, b));

		if (!Arrays.equals(a.getSize(), b.getSize()))
			throw new JOMException("Operator '^'. Size mismatch: a and b should be of the same size, or one of them should be a scalar expression");
		this.resize(a.getSize());

		this.isLinear = a.isConstant() && b.isConstant();
		this.isConstant = a.isConstant() && b.isConstant();

		if (isLinear)
		{
			/* change (destroy) the constant data in Expression a. This object keeps no reference to Expression a, but the data remains since we 
			 * keep reference to the affine expression object of a, now modified by this function */
			double[] data_a = a.getAffineExpression().getConstantCoefArray();
			double[] data_b = b.getAffineExpression().getConstantCoefArray();
			for (int cont = 0; cont < data_a.length; cont++) data_a[cont] = Math.pow(data_a[cont], data_b[cont]);
			this.affineExp = a.getAffineExpression();
		}
		
		/* Store some partial info to be used later */
		this.constantBase = (a.isConstant());
		this.aValues = (!this.isLinear() && constantBase) ? a.evaluateConstant().elements() : null;
		this.constantExponent = (b.isConstant());
		this.bValues = (!this.isLinear() && constantExponent) ? b.evaluateConstant().elements() : null;
		if (constantBase)
			this.constantAllBases = (aValues.equals(aValues.get(0)));
		else
			this.constantAllBases = false;
		if (constantExponent)
			this.constantAllExponents = (bValues.equals(bValues.get(0)));
		else
			this.constantAllExponents = false;
		this.aValue = (this.constantAllBases) ? aValues.get(0) : -1;
		this.bValue = (this.constantAllExponents) ? bValues.get(0) : -1;

		this.a = (this.isLinear()) ? null : a; // if linear do not store the expressions. This info is already in the linear coefs matrix
		this.b = (this.isLinear()) ? null : b; // if linear do not store the expressions. This info is already in the linear coefs matrix

	}

	@Override
	DoubleMatrixND nl_evaluate(double[] valuesDVs)
	{
		DoubleMatrix1D a_value = (constantBase) ? aValues : a.evaluate_internal(valuesDVs).elements();
		DoubleMatrix1D b_value = (constantExponent) ? bValues : b.evaluate_internal(valuesDVs).elements();

		a_value.assign(b_value, DoubleFunctions.pow);
		return new DoubleMatrixND(this.getSize(), a_value);
	}

	@Override
	DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs)
	{
		DoubleMatrix1D a_value = (constantBase) ? aValues : a.evaluate_internal(valuesDVs).elements();
		DoubleMatrix1D b_value = (constantExponent) ? bValues : b.evaluate_internal(valuesDVs).elements();
		DoubleMatrix2D a_jacobian = a.evaluateJacobian_internal(valuesDVs);
		DoubleMatrix2D b_jacobian = b.evaluateJacobian_internal(valuesDVs);

		if (constantAllExponents)
		{
			a_value.assign(DoubleFunctions.pow(bValue - 1));
			a_value.assign(DoubleFunctions.mult(bValue));
			return DoubleFactory2D.sparse.diagonal(a_value).zMult(a_jacobian, null);
		}
		if (constantExponent)
		{
			DoubleMatrix1D exponentsMinus1 = bValues.copy().assign(DoubleFunctions.minus(1));
			a_value.assign(exponentsMinus1, DoubleFunctions.pow);
			a_value.assign(bValues, DoubleFunctions.mult);
			return DoubleFactory2D.sparse.diagonal(a_value).zMult(a_jacobian, null);
		}
		if (constantBase)
		{
			DoubleMatrix1D lnBase = this.aValues.copy().assign(DoubleFunctions.log); // ln (a)
			a_value.assign(b_value, DoubleFunctions.pow); // a^b
			a_value.assign(this.aValues, DoubleFunctions.mult); // ln(a) * a^b
			return DoubleFactory2D.sparse.diagonal(a_value).zMult(b_jacobian, null); // // ln(a) * a^b
		}
		/* General case: neither base nor exponenets are constant */
		DoubleMatrix1D f_pow_g = a_value.copy().assign(b_value, DoubleFunctions.pow); // f^g
		DoubleMatrix1D lnBase = a_value.copy().assign(DoubleFunctions.log); // ln (f)
		DoubleMatrix1D lnBase_by_f_pow_g = lnBase.assign(f_pow_g, DoubleFunctions.mult); // ln (f) * f^g
		DoubleMatrix1D g_div_f = b_value.copy().assign(a_value, DoubleFunctions.div); // g/f
		DoubleMatrix1D g_div_f_by_f_pow_g = g_div_f.assign(f_pow_g, DoubleFunctions.mult); // (g/f) * f^g

		DoubleMatrix2D sum1 = DoubleFactory2D.sparse.diagonal(lnBase_by_f_pow_g).zMult(b_jacobian, null); // ln(f)*(f^g)*g'
		DoubleMatrix2D sum2 = DoubleFactory2D.sparse.diagonal(g_div_f_by_f_pow_g).zMult(a_jacobian, null); // (g/f)*(f^g)*f'
		return sum1.assign(sum2, DoubleFunctions.plus);
	}

	/* (non-Javadoc)
	 * @see Expression#isLinear()
	 */
	@Override
	boolean isLinear()
	{
		return isLinear;
	}

	/* (non-Javadoc)
	 * @see Expression#isConstant()
	 */
	@Override
	boolean isConstant()
	{
		return isConstant;
	}

	/* (non-Javadoc)
	 * @see Expression#getActiveVarIds()
	 */
	@Override
	LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds()
	{
		/* Make a copy since I am modifying it */
		LinkedHashMap<Integer, HashSet<Integer>> newActiveVarIdsMatrix = new LinkedHashMap<Integer, HashSet<Integer>>();
		for (Entry<Integer, HashSet<Integer>> e : this.a.getActiveVarIds().entrySet())
			newActiveVarIdsMatrix.put(e.getKey(), (HashSet<Integer>) e.getValue().clone());
		for (Entry<Integer, HashSet<Integer>> e : this.b.getActiveVarIds().entrySet())
		{
			HashSet<Integer> hsRes = newActiveVarIdsMatrix.get(e.getKey());
			if (hsRes == null)
				newActiveVarIdsMatrix.put(e.getKey(), (HashSet<Integer>) e.getValue().clone());
			else
				hsRes.addAll(e.getValue());
		}
		return newActiveVarIdsMatrix;
	}
}
