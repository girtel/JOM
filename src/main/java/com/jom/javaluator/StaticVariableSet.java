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

import java.util.HashMap;
import java.util.Map;

/** A static variable set. <br> Here, static means that the values of variables are set before starting to evaluate the expressions.
 * 
 * @param <T> The type of the values of the variable (the one handled by the evaluator).
 * @author Jean-Marc Astesana
 * @see <a href="../../../license.html">License information</a> */
public class StaticVariableSet<T> implements AbstractVariableSet<T>
{
	private final Map<String, T> varToValue;

	/** Constructor. <br> Builds a new empty variable set. */
	public StaticVariableSet()
	{
		this.varToValue = new HashMap<String, T>();
	}

	@Override
	public T get(String variableName)
	{
		return this.varToValue.get(variableName);
	}

	/** Sets a variable value.
	 * 
	 * @param variableName The variable name
	 * @param value The variable value (null to remove a variable from the set). */
	public void set(String variableName, T value)
	{
		this.varToValue.put(variableName, value);
	}
}
