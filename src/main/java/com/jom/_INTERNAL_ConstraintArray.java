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



class _INTERNAL_ConstraintArray 
{
	private final String connector;
	private final _INTERNAL_ExpressionParser evaluator;
	private final int index_0_FirstConstraint;
	private final String lhs;
	private Expression lhsMinusRhs;
	private final String name;

	private final int numScalarConstraints;

	private final String originalExpression;
	private final String rhs;
	private final int [] size;
	private final String stringLhsMinusRhs;
	public static String[] ACCEPTEDSYMBOLS_EQUAL = { "==" };
	public static String[] ACCEPTEDSYMBOLS_GREATEREQUAL = { ">=", "=>" };
	public static String[] ACCEPTEDSYMBOLS_LESSEQUAL = { "<=", "=<" };
	public static String SYMBOL_EQUAL = "==";
	public static String SYMBOL_GREATEREQUAL = ">=";

	public static String SYMBOL_LESSEQUAL = "<=";

	_INTERNAL_ConstraintArray(_INTERNAL_ExpressionParser evaluator, String name, String expression, int index_0_FirstConstraint)
	{
		this.evaluator = evaluator;
		this.name = name;
		this.index_0_FirstConstraint = index_0_FirstConstraint;
		this.originalExpression = expression;
		String[] parsedConstraint = parseConstraint(expression); // some syntax
																															// errors are
																															// checked here
		this.lhs = parsedConstraint[0];
		this.connector = parsedConstraint[1];
		this.rhs = parsedConstraint[2];
		this.stringLhsMinusRhs = "(" + lhs + ") - (" + rhs + ") ";
		this.lhsMinusRhs = this.evaluator.evaluate(stringLhsMinusRhs, null);

		this.size = lhsMinusRhs.getSize();
		this.numScalarConstraints = lhsMinusRhs.getNumScalarExpressions();
	}

	String getConnector()
	{
		return connector;
	}

	int getIndex_0_FirstConstraint()
	{
		return index_0_FirstConstraint;
	}

	int getIndex_0_LastConstraint()
	{
		return index_0_FirstConstraint + numScalarConstraints - 1;
	}

	String getLhs()
	{
		return lhs;
	}

	Expression getLhsMinusRhs()
	{
		return this.lhsMinusRhs;
	}

	String getName()
	{
		return name;
	}

	int getNumScalarConstraints()
	{
		return numScalarConstraints;
	}

	String getOriginalExpression()
	{
		return originalExpression;
	}

	String getRhs()
	{
		return rhs;
	}

	int [] getSize()
	{
		return size;
	}

	String getStringLhsMinusRhs()
	{
		return this.stringLhsMinusRhs;
	}

	boolean isLinear()
	{
		return this.lhsMinusRhs.isLinear();
	}

	private String[] parseConstraint(String constraintExpression)
	{
		String lhs = null;
		String rhs = null;
		String connector = null;

		for (int cont = 0; cont < ACCEPTEDSYMBOLS_LESSEQUAL.length; cont++)
		{
			String separator = ACCEPTEDSYMBOLS_LESSEQUAL[cont];
			String[] result = constraintExpression.split(separator);
			if (result.length == 1) continue;
			if ((result.length == 0) || (result.length > 2)) throw new JOMException ("Syntax error in the JOM expression");
			/* there are lhs and rhs */
			if (lhs != null) throw new JOMException ("Syntax error in the JOM expression. More than one delimiter");
			lhs = result[0];
			rhs = result[1];
			connector = SYMBOL_LESSEQUAL;
		}

		for (int cont = 0; cont < ACCEPTEDSYMBOLS_GREATEREQUAL.length; cont++)
		{
			String separator = ACCEPTEDSYMBOLS_GREATEREQUAL[cont];
			String[] result = constraintExpression.split(separator);
			if (result.length == 1) continue;
			if ((result.length == 0) || (result.length > 2)) throw new JOMException ("Syntax error in the JOM expression");
			/* there is lhs and rhs */
			if (lhs != null) throw new JOMException ("Syntax error in the JOM expression. More than one delimiter");
			lhs = result[0];
			rhs = result[1];
			connector = SYMBOL_GREATEREQUAL;
		}

		for (int cont = 0; cont < ACCEPTEDSYMBOLS_EQUAL.length; cont++)
		{
			String separator = ACCEPTEDSYMBOLS_EQUAL[cont];
			String[] result = constraintExpression.split(separator);
			if (result.length == 1) continue;
			if ((result.length == 0) || (result.length > 2)) throw new JOMException ("Syntax error in the JOM expression");
			/* there is lhs and rhs */
			if (lhs != null) throw new JOMException ("Syntax error in the JOM expression. More than one delimiter");
			lhs = result[0];
			rhs = result[1];
			connector = SYMBOL_EQUAL;
		}

		if (lhs == null) throw new JOMException ("Syntax error in the JOM expression");

		return new String[] { lhs, connector, rhs };

	}

	void releaseMemoryJOMExpression()
	{
		this.lhsMinusRhs = null;
	}

}
