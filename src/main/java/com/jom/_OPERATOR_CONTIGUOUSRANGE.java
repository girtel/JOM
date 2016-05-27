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

import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;

class _OPERATOR_CONTIGUOUSRANGE extends _FUNCTION_LINEAREXPRESSION
{
	_OPERATOR_CONTIGUOUSRANGE(OptimizationProblem model, Expression init, Expression end)
	{
		super(model, computeConstant(init, end));
	}

	private static DoubleMatrixND computeConstant(Expression init, Expression end)
	{
		if ((!init.isScalarConstant()) || (!end.isScalarConstant())) throw new JOMException("Operator ':' (range definition): A constant expression is expected in a range definition");

		int initValue = (int) init.evaluateConstant().get(0);
		int endValue = (int) end.evaluateConstant().get(0);
		IntMatrix1D size = IntFactory1D.dense.make(new int[] { 1, endValue - initValue + 1 });

		return new DoubleMatrixND(size.toArray(), "dense").ascending(initValue, 1.0);
	}

}
