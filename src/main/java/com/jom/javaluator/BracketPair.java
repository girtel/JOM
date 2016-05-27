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



package com.jom.javaluator;

/** A <a href="http://en.wikipedia.org/wiki/Bracket_(mathematics)">bracket pair</a>.
 * 
 * @author Jean-Marc Astesana
 * @see <a href="../../../license.html">License information</a> */
public class BracketPair
{
	private String close;
	private String open;
	/** The angle brackets pair: <>. */
	public static final BracketPair ANGLES = new BracketPair('<', '>');
	/** The braces pair: {}. */
	public static final BracketPair BRACES = new BracketPair('{', '}');

	/** The square brackets pair: []. */
	public static final BracketPair BRACKETS = new BracketPair('[', ']');

	/** The parentheses pair: (). */
	public static final BracketPair PARENTHESES = new BracketPair('(', ')');

	/** Constructor.
	 * 
	 * @param open The character used to open the brackets.
	 * @param close The character used to close the brackets. */
	public BracketPair(char open, char close)
	{
		super();
		this.open = new String(new char[] { open });
		this.close = new String(new char[] { close });
	}

	/** Gets the close bracket character.
	 * 
	 * @return a char */
	public String getClose()
	{
		return close;
	}

	/** Gets the open bracket character.
	 * 
	 * @return a char */
	public String getOpen()
	{
		return open;
	}

	@Override
	public String toString()
	{
		return open + close;
	}
}
