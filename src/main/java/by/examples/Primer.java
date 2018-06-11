package by.examples;
import by.bsu.JVmipcl.*;
import java.io.IOException;

public class Primer
{
	static {
            System.loadLibrary("JVmipcl");
        }

   	public static void main(String[] args) {
		final double[] c = {100, 64};
		final double[] b = {250, 4};
		final double[][] A = {{50, 31}, {-3,2}};
		try {
			Primer pr = new Primer();
			pr.solveMIP(c,b, A);
		}
		catch(Throwable e) {
			System.out.println(e);
		}
	}

	void solveMIP(double[] c, double[] b, double[][] A) throws IOException
	{
		int m = b.length; // m - number of rows
		int n = c.length; // n - number of columns
		int[] ind = new int[n]; // ind - auxiliary array

		MIP prob = new MIP("MIPCLtest"); // create a new MIP problem

		prob.openMatrix(n,m,m*n); // open a dence mxn-matrix

		for (int j=0; j < n; ++j) { // add n variables:
			ind[j] = j;
			prob.addVar(j,MIP.VAR_INT, c[j],0.0, LP.VAR_INF);
		}

		for (int i=0; i < m; ++i) // add m rows (constraints):
			prob.addRow(i,0,-LP.INF,b[i],A[i],ind,true);

		prob.closeMatrix(); // close the matrix

		prob.optimize(); // solve the problem

		prob.printSolution("primer.sol"); // print the solution to the file
		prob.dispose(); // delete the CMIP object referenced by prob
	}
}

