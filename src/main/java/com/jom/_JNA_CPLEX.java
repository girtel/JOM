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

import com.sun.jna.Library;
import com.sun.jna.Pointer;

/** @author Pablo */
public interface _JNA_CPLEX extends Library
{
	public static int CPX_MAX                 = -1;
	public static int CPX_MIN                 = 1;
	public static int CPX_OFF                 = 0;
	public static int CPX_ON                  = 1;
	//
	public static int CPX_PARAM_SCRIND        = 1035;
	public static int CPX_PARAM_TUNINGDISPLAY = 1113;
	public static int CPX_PARAM_TILIM         = 1039; // maximum execution time in seconds
	public static int CPX_PARAM_EPGAP         = 2009; // optimality gap in MIPs (p.e. 5% => 0.05). stops when it is reached
	public static int CPX_PARAM_EPOPT         = 1014; // optimality gap for simplex (p.e. 5% => 0.05). stops when it is reached

	public int CPXsetstrparam(Pointer env, int whichparam, byte[] newvalue_str);

	int CPXchgcoeflist(Pointer env, Pointer lp, int numcoefs, int[] rowlist, int[] collist, double[] vallist);

	void CPXchgobjsen(Pointer env, Pointer lp, int maxormin); // set maximization or minimization

	Pointer CPXcreateprob(Pointer env, int[] status_p, byte[] probname_str);

	int CPXfreeprob(Pointer env, Pointer lp_p);

	int CPXgetcoef(Pointer env, Pointer lp, int i, int j, double[] coef_p);

	int CPXgetctype(Pointer env, Pointer lp, byte[] xctype, int begin, int end);

	int CPXgetdj(Pointer env, Pointer lp, double[] dj, int begin, int end);

	Pointer CPXgeterrorstring(Pointer env, int errcode, byte[] buffer_str); // CPXCCHARptr CPXgeterrorstring(CPXCENVptr env, int errcode, char *
	// buffer_str)

	int CPXgetnumcols(Pointer env, Pointer lp);

	int CPXgetnumrows(Pointer env, Pointer lp);

	int CPXgetobjval(Pointer env, Pointer lp, double[] objval);

	int CPXgetpi(Pointer env, Pointer lp, double[] pi, int begin, int end);

	int CPXgetrhs(Pointer env, Pointer lp, double[] rhs, int begin, int end);

	int CPXgetslack(Pointer env, Pointer net, double[] slack, int begin, int end);

	int CPXgetstat(Pointer env, Pointer lp);

	int CPXgetx(Pointer env, Pointer lp, double[] x, int begin, int end);

	int CPXlpopt(Pointer env, Pointer lp);

	int CPXmipopt(Pointer env, Pointer lp);

	int CPXnewcols(Pointer env, Pointer lp, int ccnt, double[] obj, double[] lb, double[] ub, byte[] xctype, byte[][] colname);

	int CPXnewrows(Pointer env, Pointer lp, int rcnt, double[] rhs, byte[] sense, double[] rngval, byte[][] rowname);

	Pointer CPXopenCPLEX(int[] status); // receives pointer int, returns pointer to CPXENVptr --> this struct no info about their fields

	int CPXsetdblparam(Pointer env, int whichparam, double newvalue);

	int CPXsetintparam(Pointer env, int whichparam, int newvalue);

	int CPXsolninfo(Pointer env, Pointer lp, int[] solnmethod_p, int[] solntype_p, int[] pfeasind_p, int[] dfeasind_p);

	int CPXsolution(Pointer env, Pointer lp, int[] lpstat_p, double[] objval_p, double[] x, double[] pi, double[] slack, double[] dj);

	int CPXwriteprob(Pointer env, Pointer lp, byte[] filename_str, byte[] filetype_str);

}
