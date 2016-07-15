
import com.dashoptimization.*;

import cern.colt.Arrays;

// 1. The name of the file where the license is (xpauth.xpr) should be given in the initialization
// 2. 

public class TestXPRESS 
{
 public static void main (String [] args)
 {
	 System.out.println("A");
	 
	 XPRS.init("C:\\xpressmp\\xpauth.xpr");
	 
	 XPRSprob p = new XPRSprob();
	 
//	 final int ncol = 2; // num variables?
//	 final int ncoeffs = 4; // num coeffs matrix?
	 final byte[] _srowtypes = new byte [] { 'R' , 'R'};
//	 final double[] _dobj = new double [] { 1,1}; // func obj?
//	 final int[] _mstart = new int [] { 0 , 0,1,1}; // rows of the matrix?
//	 final int[] _mrwind = new int [] { 0 , 1,0,1}; // cols of the matrix?
//	 final double[] _dmatval = new double [] { 1,-1,1,1}; // vals of the coeffs?
//	 final double[] _dbdl = new double [] { 0 , 0}; // lower bounds constraints
//	 final double[] _dbdu = new double [] { 1 , 2}; // upper bounds constraints
	 
	 //loadGlobal(int[] _mgcols, double[] _dlim, byte[] _stype, int[] _msstart, int[] _mscols, double[] _dref) 
	 
	 p.loadGlobal("Problem name", 2, 2, _srowtypes , new double [] {1,1.5} , new double [] {1,1.5},new double [] {1,1} , 
			 new int [] {0,2} , new int [] {2,2} , new int [] {0,1,0,1}, new double [] {1,1,-1,1}, new double [] {0,0} , 
			 new double [] {XPRSconstants.PLUSINFINITY , XPRSconstants.PLUSINFINITY } , 
			 2 , 0 , new byte [] {'I' , 'I'} , 
			 new int [] {0,1} , new double [2] , new byte [2] , null ,  null , null);
	 
	 IntHolder status = new IntHolder();
	 p.loadMipSol(new double [] {-1 ,-1} , status);
	 System.out.println("Setting the initial solution status: " + status);
	 
	 p.chgObjSense(XPRSenumerations.ObjSense.MAXIMIZE);
	 
	 p.mipOptimize();

	 
	 double [] primalSolution = new double [2];
	 double [] slackSolution = new double [2];
	 p.getMipSol(primalSolution , slackSolution);

	 System.out.println("Primal: " + Arrays.toString(primalSolution));
	 System.out.println("Slacks: " + Arrays.toString(slackSolution));
	 //	 p.initGlobal();
//	 p.addCols(ncols, ncoeffs, _dobj, _mstart, _mrwind, _dmatval, _dbdl, _dbdu);
	 
	 
	 /* */
	 XPRS.free(); // frees the license
	 
	 System.out.println("Ok");
 }
}
