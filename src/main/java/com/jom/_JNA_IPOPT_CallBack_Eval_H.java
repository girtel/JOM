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



 




/**
 * 
 */
package com.jom;

import com.sun.jna.Callback;
import com.sun.jna.Pointer;

/** @author Pablo */
public interface _JNA_IPOPT_CallBack_Eval_H extends Callback
{
	public boolean callback(int n, Pointer x, boolean new_x, double obj_factor, int m, Pointer lambda, boolean new_lambda, int nele_hess, Pointer iRow, Pointer jCol, Pointer values, Pointer user_data);
}
