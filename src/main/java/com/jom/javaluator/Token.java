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

package com.jom.javaluator;

/** A token. <br> When evaluating an expression, it is first split into tokens. These tokens can be operators, constants, etc ...
 *
 * @author Jean-Marc Astesana
 * @see <a href="../../../license.html">License information</a> */
public class Token
{
	static final Token FUNCTION_ARG_SEPARATOR = new Token(Kind.FUNCTION_SEPARATOR, null);
	private Object content;
	private Kind kind;

	private Token(Kind kind, Object content)
	{
		super();
		if ((kind.equals(Kind.OPERATOR) && !(content instanceof Operator)) || (kind.equals(Kind.FUNCTION) && !(content instanceof Function)) ||
				(kind.equals(Kind.LITERAL) && !(content instanceof String)))
			throw new IllegalArgumentException();
		this.kind = kind;
		this.content = content;
	}

	static Token buildCloseToken(BracketPair pair)
	{
		return new Token(Kind.CLOSE_BRACKET, pair);
	}

	static Token buildFunction(Function function)
	{
		return new Token(Kind.FUNCTION, function);
	}

	static Token buildLiteral(String literal)
	{
		return new Token(Kind.LITERAL, literal);
	}

	static Token buildOpenToken(BracketPair pair)
	{
		return new Token(Kind.OPEN_BRACKET, pair);
	}

	static Token buildOperator(Operator ope)
	{
		return new Token(Kind.OPERATOR, ope);
	}

	/** Tests whether the token is a close bracket.
	 *
	 * @return true if the token is a close bracket */
	public boolean isCloseBracket()
	{
		return kind.equals(Kind.CLOSE_BRACKET);
	}

	/** Tests whether the token is a function.
	 *
	 * @return true if the token is a function */
	public boolean isFunction()
	{
		return kind.equals(Kind.FUNCTION);
	}

	/** Tests whether the token is a function argument separator.
	 *
	 * @return true if the token is a function argument separator */
	public boolean isFunctionArgumentSeparator()
	{
		return kind.equals(Kind.FUNCTION_SEPARATOR);
	}

	/** Tests whether the token is a literal or a constant or a variable name.
	 *
	 * @return true if the token is a literal, a constant or a variable name */
	public boolean isLiteral()
	{
		return kind.equals(Kind.LITERAL);
	}

	/** Tests whether the token is an open bracket.
	 *
	 * @return true if the token is an open bracket */
	public boolean isOpenBracket()
	{
		return kind.equals(Kind.OPEN_BRACKET);
	}

	/** Tests whether the token is an operator.
	 *
	 * @return true if the token is an operator */
	public boolean isOperator()
	{
		return kind.equals(Kind.OPERATOR);
	}

	Operator.Associativity getAssociativity()
	{
		return getOperator().getAssociativity();
	}

	BracketPair getBrackets()
	{
		return (BracketPair) this.content;
	}

	Function getFunction()
	{
		return (Function) this.content;
	}

	Kind getKind()
	{
		return kind;
	}

	String getLiteral()
	{
		if (!this.kind.equals(Kind.LITERAL)) throw new IllegalArgumentException();
		return (String) this.content;
	}

	Operator getOperator()
	{
		return (Operator) this.content;
	}

	int getPrecedence()
	{
		return getOperator().getPrecedence();
	}

	public enum Kind
	{ // CHANGE: Kind was private (created some compilation errors)
		CLOSE_BRACKET, FUNCTION, FUNCTION_SEPARATOR, LITERAL, OPEN_BRACKET, OPERATOR
	}
}
