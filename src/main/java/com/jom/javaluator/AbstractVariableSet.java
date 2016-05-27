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

/** An abstract variable set. <br> Javaluator supports expression that contains variables (for example <i>sin(x)</i>). <br> An AbstractVariableSet converts, during the expression evaluation, each variable to its value.
 * 
 * @param <T> The type of the values of the variable (the one handled by the evaluator).
 * @author Jean-Marc Astesana
 * @see <a href="../../../license.html">License information</a> */
public interface AbstractVariableSet<T>
{
	/** Gets the value of a variable.
	 * 
	 * @param variableName The name of a variable
	 * @return the variable's value or null if the variable is unknown */
	public T get(String variableName);
}
