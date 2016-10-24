import cern.colt.matrix.Transpose;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import com.jom.*;
import com.jom.javaluator.*;

import java.lang.Object;



public class Example {

	public static void main(String[] args) {
				
		double[] x_i = {0, 10, 0, 0}; 		
		double[] x_u = {400, 400, 400, 300}; 
		double[][] Q = {{0.009, 0, 0, 0}, {0, 0.008, 0, 0},{0, 0, 0.008, 0},{ 0, 0, 0, 0}};
		double[] b = {80, 90, 70, 80};
		double[][] Bloss= {{0.0002,0.0001,0.,0.0001}, {0.0001,0.0002,0.0001,0.0001},{0.,0.0001,0.0002,0.0001},{0.,0.0001,0.0002,0.0003}};
		double[] E_trans1 = {1,1,0,0};
		double[] E_trans2 = {0,0,1,1};
		double[] U = {1,1,1,1};		
			
		DoubleMatrixND x_ii = new DoubleMatrixND(new int[] {4, 1}, x_i);
		DoubleMatrixND x_uu = new DoubleMatrixND(new int[] {4, 1}, x_u);		
		DoubleMatrixND QQ = new DoubleMatrixND(Q);	
		DoubleMatrixND bb = new DoubleMatrixND(new int[] {4, 1}, b);	
		DoubleMatrixND Blosss = new DoubleMatrixND(Bloss);
		DoubleMatrixND E_tra1 = new DoubleMatrixND(new int[] {1, 4}, E_trans1);
		DoubleMatrixND E_tra2 = new DoubleMatrixND(new int[] {1, 4}, E_trans2);
		DoubleMatrixND UU = new DoubleMatrixND(new int[] {4, 1}, U);
		DoubleMatrixND PG0 = new DoubleMatrixND(new int[] {4, 1}, new double[] {200, 200, 200, 200});
		
		OptimizationProblem Opt = new OptimizationProblem();		
		Opt.addDecisionVariable("PG", false, new int[] {4, 1}, x_ii, x_uu);

		Opt.setInputParameter("QQ", QQ);		
		Opt.setInputParameter("bb", bb);
		Opt.setInputParameter("Blosss", Blosss);
		Opt.setInputParameter("E_trans1", E_tra1);
		Opt.setInputParameter("E_trans2", E_tra2);
		Opt.setInputParameter("UU", UU);
		Opt.setInputParameter("AA", UU);
		Opt.setInputParameter("x_ii", x_ii);
		Opt.setInputParameter("x_uu", x_uu);
		Opt.setInputParameter("PG0", PG0);
		
		Opt.setObjectiveFunction("minimize",  "PG'*QQ*PG+bb'*PG+430");
		
		Opt.addConstraint("sum(PG) - PG' * Blosss * PG - 400 == 0");  //Here is the problem line. If I take out this part:<- PG' * Blosss * PG > the code runs without problem!!!!!!!!!!?????????!!!!!!!!!!!!!!!
		Opt.setInitialSolution("PG", PG0);		
		
		//Opt.solve("ipopt", "solverLibraryName" , "ipopt.dll", "mu_init", 1.0, "print_level", 12);
		Opt.solve("ipopt", "solverLibraryName" , "ipopt.dll", "mu_init", 1.0);
		if (!Opt.solutionIsOptimal ()) throw new RuntimeException ("An optimal solution was not found");
		DoubleMatrixND XX = Opt.getPrimalSolution("PG");
		System.out.println(XX);
		
	}
}