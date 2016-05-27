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

import com.sun.jna.Library;
import com.sun.jna.Pointer;

/** @author Pablo */
public interface _JNA_IPOPT extends Library
{

	public boolean AddIpoptIntOption(Pointer ipopt_problem, String keyword, int val);

	public boolean AddIpoptNumOption(Pointer ipopt_problem, String keyword, double val);

	public boolean AddIpoptStrOption(Pointer ipopt_problem, String keyword, String val);

	public Pointer CreateIpoptProblem(int n, double[] x_L, double[] x_U, int m, double[] g_L, double[] g_U, int nele_jac, int nele_hess, int index_style, _JNA_IPOPT_CallBack_Eval_F eval_f, _JNA_IPOPT_CallBack_Eval_G eval_g, _JNA_IPOPT_CallBack_Eval_Grad_F eval_grad_f, _JNA_IPOPT_CallBack_Eval_Jac_G eval_jac_g, _JNA_IPOPT_CallBack_Eval_H eval_h);

	public void FreeIpoptProblem(Pointer ipopt_problem);

	public int IpoptSolve(Pointer ipopt_problem, double[] x, double[] g, double[] obj_val, double[] mult_g, double[] mult_x_L, double[] mult_x_U, Pointer user_data);

	public boolean SetIpoptProblemScaling(Pointer ipopt_problem, double obj_scaling, double[] x_scaling, double[] g_scaling);

}
