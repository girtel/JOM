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
import com.sun.jna.Structure;

import java.util.Arrays;
import java.util.List;

/** @author Pablo */
public interface _JNA_GLPK extends Library
{
	public static int GLP_MAX = 2; /* maximization */

	public static int GLP_MIN = 1; /* minimization */

	public static int GLP_MSG_OFF = 0; /* disable something */

	public static int GLP_OFF = 0; /* disable something */

	public static int GLP_ON = 1; /* enable something */

	int glp_add_cols(Pointer lp, int ncs);

	int glp_add_rows(Pointer lp, int nrs);

	Pointer glp_create_prob();

	void glp_delete_prob(Pointer lp);

	double glp_get_col_dual(Pointer lp, int j);

	double glp_get_col_prim(Pointer lp, int j);

	int glp_get_num_cols(Pointer lp);

	int glp_get_num_rows(Pointer lp);

	double glp_get_obj_val(Pointer lp);

	double glp_get_row_dual(Pointer lp, int i);

	double glp_get_row_prim(Pointer lp, int i);

	int glp_get_status(Pointer lp);

	int glp_init_env();

	/* MIP solver */
	void glp_init_iocp(Pointer parm);

	/* INTERIOR POINT SOLVER */
	int glp_init_iptcp(Pointer parm);

	/* SIMPLEX SOLVER */
	int glp_init_smcp(Pointer parm);

	int glp_interior(Pointer P, Pointer parm);

	int glp_intopt(Pointer mip, Pointer parm);

	double glp_ipt_col_dual(Pointer lp, int j);

	double glp_ipt_col_prim(Pointer lp, int j);

	double glp_ipt_obj_val(Pointer lp);

	double glp_ipt_row_dual(Pointer lp, int i);

	double glp_ipt_row_prim(Pointer lp, int i);

	int glp_ipt_status(Pointer lp);

	void glp_load_matrix(Pointer lp, int ne, int[] ia, int[] ja, double[] ar);

	double glp_mip_col_val(Pointer mip, int j);

	double glp_mip_obj_val(Pointer mip);

	double glp_mip_row_val(Pointer mip, int i);

	int glp_mip_status(Pointer mip);

	void glp_set_col_bnds(Pointer lp, int j, int type, double lb, double ub);

	void glp_set_col_kind(Pointer mip, int j, int kind);

	void glp_set_obj_coef(Pointer lp, int j, double coef);

	void glp_set_obj_dir(Pointer lp, int dir);

	void glp_set_row_bnds(Pointer lp, int i, int type, double lb, double ub);

	int glp_simplex(Pointer lp, Pointer parm);

	public static class GLP_IOCP extends Structure
	{
		public int     alien; /* use alien solver */
		public int     binarize; /* try to binarize integer variables */
		public int     br_tech; /* branching technique: */
		public int     bt_tech; /* backtracking technique: */
		public Pointer cb_func; // void (*cb_func)(glp_tree *T, void *info);
		public Pointer cb_info; // void *cb_info; /* mip.cb_info */
		public int     cb_size; /* mip.cb_size */
		public int     clq_cuts; /* clique cuts (GLP_ON/GLP_OFF) */
		public int     cov_cuts; /* cover cuts (GLP_ON/GLP_OFF) */
		public double foo_bar[] = new double[29]; /* (reserved) */
		public int    fp_heur; /* feasibility pump heuristic */
		public int    gmi_cuts; /* Gomory's cuts (GLP_ON/GLP_OFF) */
		public double mip_gap; /* relative MIP gap tolerance */
		public int    mir_cuts; /* MIR cuts (GLP_ON/GLP_OFF) */
		public int    msg_lev; /* message level (see glp_smcp) */
		public int    out_dly; /* mip.out_dly (milliseconds) */
		public int    out_frq; /* mip.out_frq (milliseconds) */
		public int    pp_tech; /* preprocessing technique: */
		public int    presolve; /* enable/disable using MIP presolver */
		public int    tm_lim; /* mip.tm_lim (milliseconds) */
		public double tol_int; /* mip.tol_int */
		public double tol_obj; /* mip.tol_obj */

		@Override
		protected List getFieldOrder()
		{
			return Arrays.asList(new String[]{"msg_lev", "br_tech", "bt_tech", "tol_int", "tol_obj", "tm_lim", "out_frq", "out_dly", "cb_func",
					"cb_info", "cb_size", "pp_tech", "mip_gap", "mir_cuts", "gmi_cuts", "cov_cuts", "clq_cuts", "presolve", "binarize", "fp_heur",
					"alien", "foo_bar"});
		}
	}

	public static class GLP_IPTCP extends Structure
	{
		public double foo_bar[] = new double[48]; /* (reserved) */
		public int msg_lev; /* message level (see glp_smcp) */
		public int ord_alg; /* ordering algorithm: */

		@Override
		protected List getFieldOrder()
		{
			return Arrays.asList(new String[]{"msg_lev", "ord_alg", "foo_bar"});
		}
	}

	public static class GLP_PROB extends Structure
	{
		public double _opaque_prob[] = new double[100];

		@Override
		public List getFieldOrder()
		{
			return Arrays.asList(new String[]{"_opaque_prob"});
		}
	}

	public static class GLP_SMCP extends Structure
	{
		public double foo_bar[] = new double[36]; /* (reserved) */
		public int    it_lim; /* spx.it_lim */
		public int    meth; /* simplex method option: */
		public int    msg_lev; /* message level: */
		public double obj_ll; /* spx.obj_ll */
		public double obj_ul; /* spx.obj_ul */
		public int    out_dly; /* spx.out_dly (milliseconds) */
		public int    out_frq; /* spx.out_frq */
		public int    presolve; /* enable/disable using LP presolver */
		public int    pricing; /* pricing technique: */
		public int    r_test; /* ratio test technique: */
		public int    tm_lim; /* spx.tm_lim (milliseconds) */
		public double tol_bnd; /* spx.tol_bnd */
		public double tol_dj; /* spx.tol_dj */
		public double tol_piv; /* spx.tol_piv */

		@Override
		protected List getFieldOrder()
		{
			return Arrays.asList(new String[]{"msg_lev", "meth", "pricing", "r_test", "tol_bnd", "tol_dj", "tol_piv", "obj_ll", "obj_ul", "it_lim",
					"tm_lim", "out_frq", "out_dly", "presolve", "foo_bar"});
		}
	}

}
