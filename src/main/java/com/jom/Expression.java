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

import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Map.Entry;

/** Expressions objects represent a function of the decision variables defined in its related OptimizationProblem object. JOM creates Expression
 * objects parsing Strings in the JOM syntax.
 * @author Pablo Pavon Mariño
 * @see http://www.net2plan.com/jom
 * */
public abstract class Expression
{
	protected final OptimizationProblem model;
	protected final int                 numVarIds;
	//protected ExpressionCurvatureInfo curvature;
	protected _INTERNAL_AffineExpressionCoefs affineExp;
	protected int   numScalarExpressions;
	protected int   numDim; // reshaping can change this
	protected int[] size; // reshaping can change this
	private   DoubleMatrixND                  cachedValueIfConstant;

	Expression(OptimizationProblem model)
	{
		this.model = model;
		this.numVarIds = model.getNumScalarDecisionVariables();
		this.cachedValueIfConstant = null;
		this.affineExp = null;
		this.numDim = -1;
		this.numScalarExpressions = -1;
		this.size = null;
	}

	/* A constant scalar */
	Expression(OptimizationProblem model, double constant)
	{
		this(model, new DoubleMatrixND(constant));
	}

	/* A linear expression */
	Expression(OptimizationProblem model, _INTERNAL_AffineExpressionCoefs affineExpression)
	{
		this.model = model;
		this.numVarIds = model.getNumScalarDecisionVariables();
		this.size = affineExpression.getSize();
		this.numScalarExpressions = IntMatrixND.prod(size);
		this.numDim = this.size.length;
		this.cachedValueIfConstant = null;

		if (this.numVarIds != affineExpression.getNumVarIds()) throw new JOMException("Wrong number of columns in coefs matrix");
		if (this.numScalarExpressions != affineExpression.getNumScalarExpressions()) throw new JOMException("Wrong number of rows in coefs matrix");

		this.affineExp = affineExpression;
		if (size.length < 2) throw new JOMException("1-D vectors are not allowed. Try row vectors or column vectors");
	}

	/* General constructor to be called by overriding classes */
	Expression(OptimizationProblem model, int[] size, _INTERNAL_AffineExpressionCoefs affineExpression)
	{
		this.model = model;
		this.numVarIds = model.getNumScalarDecisionVariables();
		this.size = Arrays.copyOf(size, size.length);
		this.numScalarExpressions = IntMatrixND.prod(size);
		this.numDim = size.length;
		this.cachedValueIfConstant = null;
		this.affineExp = affineExpression;

		if (size.length < 2) throw new JOMException("1-D vectors are not allowed. Try row vectors or column vectors");
	}

	/* A decision variable (and thus linear) */
	Expression(OptimizationProblem model, _INTERNAL_DecisionVariableArray dv)
	{
		this(model, dv.getVarIds());
	}

	/* A constant expression */
	Expression(OptimizationProblem model, DoubleMatrixND constant)
	{
		this.model = model;
		this.numVarIds = model.getNumScalarDecisionVariables();
		this.size = constant.getSize().toArray();
		this.numScalarExpressions = IntMatrixND.prod(size);
		this.numDim = this.size.length;
		this.cachedValueIfConstant = constant.copy();
		this.affineExp = new _INTERNAL_AffineExpressionCoefs(model, size, constant.elements().toArray());
		if (size.length < 2) throw new JOMException("1-D vectors are not allowed. Try row vectors or column vectors");
	}

	/* An array or subarray of decision variables */
	Expression(OptimizationProblem model, IntMatrixND varIds)
	{
		this.model = model;
		this.numVarIds = model.getNumScalarDecisionVariables();
		this.size = varIds.getSize().toArray();
		this.numScalarExpressions = IntMatrixND.prod(size);
		this.numDim = this.size.length;
		this.cachedValueIfConstant = null;
		this.affineExp = new _INTERNAL_AffineExpressionCoefs(model, varIds);
		if (size.length < 2) throw new JOMException("1-D vectors are not allowed. Try row vectors or column vectors");
	}

	static String toStringLinearScalarExpression(OptimizationProblem model, DoubleMatrix1D row)
	{
		double constant = row.get((int) row.size() - 1);
		String s = "" + constant;
		DoubleMatrix1D lCoefs = row.viewPart(0, (int) row.size() - 1);
		IntArrayList indexes = new IntArrayList();
		DoubleArrayList values = new DoubleArrayList();
		lCoefs.getNonZeros(indexes, values);
		indexes.trimToSize();
		values.trimToSize();
		for (int cont = 0; cont < indexes.size(); cont++)
		{
			int varId = indexes.get(cont);
			double coef = values.get(cont);
			_INTERNAL_DecisionVariableArray varInfo = model.getVarInfo(varId);
			String varName = varInfo.getName();
			int indexInVariable = varId - varInfo.getFirstVarId();
			IntMatrix1D subindexesInVariable = IntMatrixND.ind2sub(indexInVariable, varInfo.getSize());
			s = s + " + " + coef + "�" + varName + Arrays.toString(subindexesInVariable.toArray());
		}
		return s;
	}

	/** Evaluates this expression. The expression should be a constant: do not depend from the decision variables
	 * @return an array with the same size as the expression, with the values evaluated in each cell */
	public final DoubleMatrixND evaluateConstant()
	{
		if (!this.isConstant()) throw new JOMException("The expression is not constant");
		if (cachedValueIfConstant != null) return cachedValueIfConstant.copy();
		this.cachedValueIfConstant = this.evaluate_internal(new double[this.numVarIds]);
		return cachedValueIfConstant.copy();
	}

	/** Evaluates this expression, in the point provided by the decision variables. For those variables for which value is not provided, zeros are
	 * assumed.
	 * @param varNameVarValuePairs An even number of parameters, first in each pair is a String with the variable name, second a DoubleMatrixND
	 *                                object with its value
	 * @return an array with the same size as the expression, with the values evaluated in each cell */
	public final DoubleMatrixND evaluate(Object... varNameVarValuePairs)
	{
		if (varNameVarValuePairs.length % 2 != 0) throw new JOMException("Wrong number of parameters");
		double[] x = new double[this.model.getNumScalarDecisionVariables()];
		for (int cont = 0; cont < varNameVarValuePairs.length / 2; cont++)
		{
			if (!(varNameVarValuePairs[cont * 2] instanceof String)) throw new JOMException("Variable name must be a String");
			if (!(varNameVarValuePairs[cont * 2 + 1] instanceof DoubleMatrixND)) throw new JOMException("Value name must be a DoubleMatrixND "
					+ "object");
			String varName = (String) varNameVarValuePairs[cont * 2];
			DoubleMatrixND val = (DoubleMatrixND) varNameVarValuePairs[cont * 2 + 1];
			_INTERNAL_DecisionVariableArray dv = this.model.getDecisionVariable(varName);
			if (dv == null) throw new JOMException("Unknown decision variable name");
			if (!dv.getSize().equals(val.getSize())) throw new JOMException("Wrong size of variable: '" + varName + "'");
			double[] data = val.elements().toArray();
			System.arraycopy(data, 0, x, dv.getFirstVarId(), data.length);
		}
		return this.evaluate_internal(x);
	}

	/** Evaluates the given cell in this expression. The expression should be a constant: do not depend from the decision variables
	 * @return an array with the same size as the expression, with the values evaluated in each cell */
	public final double evaluateConstant(int cellIndex)
	{
		if (!this.isConstant()) throw new JOMException("The expression is not constant");
		return this.evaluate_internal(new double[this.numVarIds]).get(cellIndex);
	}

	final DoubleMatrix2D evaluateJacobian_internal(double[] valuesDVs)
	{
		if (!this.isLinear()) return this.nl_evaluateJacobian(valuesDVs);
		else return this.affineExp.getJacobian();
	}

	/** Evaluates the jacobian of this expression, in the point provided by the decision variables. For those variables for which value is not
	 * provided, zeros are assumed.
	 * @param varNameVarValuePairs An even number of parameters, first in each pair is a String with the variable name, second a DoubleMatrixND
	 *                                object with its value
	 * @return a DoubleMatrix2D array representing the jacobian matrix, with one row for each cell in the expression (in the position of its linear
	 * index), and one column per decision variable */
	public final DoubleMatrix2D evaluateJacobian(Object... varNameVarValuePairs)
	{
		if (varNameVarValuePairs.length % 2 != 0) throw new JOMException("Wrong number of parameters");
		double[] x = new double[this.model.getNumScalarDecisionVariables()];
		for (int cont = 0; cont < varNameVarValuePairs.length / 2; cont++)
		{
			if (!(varNameVarValuePairs[cont * 2] instanceof String)) throw new JOMException("Variable name must be a String");
			if (!(varNameVarValuePairs[cont * 2 + 1] instanceof DoubleMatrixND)) throw new JOMException("Value name must be a DoubleMatrixND "
					+ "object");
			String varName = (String) varNameVarValuePairs[cont * 2];
			DoubleMatrixND val = (DoubleMatrixND) varNameVarValuePairs[cont * 2 + 1];
			_INTERNAL_DecisionVariableArray dv = this.model.getDecisionVariable(varName);
			if (dv == null) throw new JOMException("Unknown decision variable name");
			if (!dv.getSize().equals(val.getSize())) throw new JOMException("Wrong size of variable: " + varName);
			double[] data = val.elements().toArray();
			System.arraycopy(data, 0, x, dv.getFirstVarId(), data.length);
		}
		return this.evaluateJacobian_internal(x);
	}

	/** Returns the number of dimensions in this expression
	 * @return the number of dimensions */
	public final int getNumDim()
	{
		return this.numDim;
	}

	/** Returns the number of cells (scalar expressions) inside this arrayed expression
	 * @return the number of cells */
	public final int getNumScalarExpressions()
	{
		return this.numScalarExpressions;
	}

	/** Returns true if the expression has two dimensions, and is a row (1xS) or column (Sx1) array. Returns false otherwise
	 * @return true or false as stated above */
	final boolean isRowOrColumn2DExpression()
	{
		if (this.numDim != 2) return false;
		return ((this.size[0] == 1) || (this.size[1] == 1));
	}

	/** Returns true if the expression has two dimensions, and is a 1x1 expression (a scalar)
	 * @return true or false as stated above */
	public final boolean isScalar()
	{
		return ((this.numScalarExpressions == 1) && (this.numDim == 2));
	}

	/** Returns true if the expression has two dimensions, is a 1x1 expression (a scalar), and also is constant (does not depend on the decision
	 * variables)
	 * @return true or false as stated above */
	final boolean isScalarConstant()
	{
		return (isScalar() && isConstant());
	}

	/** Resizes the array. The new size must have the same number of cells as the old one. The cell in index i in the old array has the same index
	 * in the new. Cell values are NOT copied.
	 * @param newSize the new size of the array */
	final void reshape(int[] newSize)
	{
		if (IntMatrixND.prod(this.size) != this.numScalarExpressions)
			throw new JOMException("Reshaping an expression cannot change its number of elements");
		this.size = newSize;
		this.numDim = this.size.length;
		if (this.affineExp != null)
			this.affineExp.reshape(newSize);
	}

	/** Returns the size of the array. The size array produced is a copy.
	 * @return Array size. The i-th element in size indicates the size of the array in the corresponding dimension. */
	public final int[] size()
	{
		return this.size;
	}

	/* A formatted String with information of the expression
	 * 
	 * @see java.lang.Object#toString() */
	@Override
	public String toString()
	{
		String s = "Expression " + ((this.isLinear()) ?
				"(linear)" :
				"(not linear)") + " size: " + Arrays.toString(size) + " (" + this.numScalarExpressions + " scalar expressions) , number variables "
				+ "model: " + numVarIds + "\n";

		if (!this.isLinear()) return s;

		for (int contCell = 0; contCell < this.numScalarExpressions; contCell++)
		{
			IntMatrix1D coord = IntMatrixND.ind2sub(contCell, IntFactory1D.dense.make(size));
			s = s + "  Cell " + Arrays.toString(coord.toArray()) + ": ";
			s = s + this.affineExp.getCellConstantCoef(contCell);
			LinkedHashMap<Integer, Double> c = this.affineExp.getScalarExpressionLinearCoefs(contCell);
			if (c != null)
				for (Entry<Integer, Double> e : c.entrySet())
				{
					int thisActiveVarId = e.getKey();
					double coef = e.getValue();
					_INTERNAL_DecisionVariableArray vSet = model.getVarInfo(thisActiveVarId);
					int indexInVariable = thisActiveVarId - vSet.getFirstVarId();
					IntMatrix1D coordInVariableThisCoef = IntMatrixND.ind2sub(indexInVariable, vSet.getSize());
					s = s + " + " + coef + " * " + vSet.getName() + Arrays.toString(coordInVariableThisCoef.toArray());
				}
			s = s + "\n";
		}
		return s;
	}

	/* tries to evaluate the expression to a number if out of domain somewhere, or not numbers => exception */
	/* This implementation works for linear function. To be overriden for other expressions */
	final DoubleMatrixND evaluate_internal(double[] valuesDVs)
	{
		if (!this.isLinear()) return this.nl_evaluate(valuesDVs); // not lineal
		return this.affineExp.evaluate(valuesDVs); // linear
	}

	/** Returns the optimization problem object that is the framework where this expression was created
	 * @return the optimization problem object */
	public final OptimizationProblem getModel()
	{
		return model;
	}

	final int[] getSize()
	{
		return this.size;
	}

	final _INTERNAL_AffineExpressionCoefs getAffineExpression()
	{
		return this.affineExp;
	}

	void resize(int[] size)
	{
		this.size = size;
		this.numDim = size.length;
		this.numScalarExpressions = IntMatrixND.prod(size);
		if (this.numScalarExpressions < 0)
			throw new JOMException("Wrong size of Expression object"); // maybe because of a big number greater than Integer.MAX_VALUE creates this
	}

	/* TO BE IMPLEMENTED BY THE SUBCLASSES. CALLED WHEN THE EXPRESSION IS NOT LINEAR */
	abstract DoubleMatrixND nl_evaluate(double[] valuesDVs);

	abstract DoubleMatrix2D nl_evaluateJacobian(double[] valuesDVs);

	abstract boolean isLinear();

	abstract boolean isConstant();

	abstract LinkedHashMap<Integer, HashSet<Integer>> nl_getActiveVarIds();

	final LinkedHashMap<Integer, HashSet<Integer>> getActiveVarIds()
	{
		if (this.isLinear()) return this.affineExp.getNonZeroLinearCoefsCols();
		else return this.nl_getActiveVarIds();
	}

}
