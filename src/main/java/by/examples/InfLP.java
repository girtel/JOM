package by.examples;
import by.bsu.JVmipcl.*;

public class InfLP
{
	static {
            System.loadLibrary("JVmipcl");
        }

   	public static void main(String[] args) {
		try {
			LP lp = new LP("infLP");
			lp.openMatrix(2,3,5);
			lp.addVar(0,0,4.0,0.0,LP.VAR_INF);
			lp.addVar(1,0,1.0,0.0,LP.VAR_INF);
			lp.addVar(2,0,2.0,0.0,LP.VAR_INF);
			lp.addCtr(0,0,-LP.INF,5.0);
			lp.addCtr(1,0,3.0,LP.INF);
			lp.addEntry(2.0,0,0);
			lp.addEntry(5.0,0,1);
			lp.addEntry(1.0,0,2);
			lp.addEntry(1.0,1,0);
			lp.addEntry(2.0,1,1);
			lp.preprocOff();
			lp.closeMatrix();
			lp.optimize();
			lp.printSolution("infLP.sol");
		}
		catch(Throwable e) {
			System.out.println(e);
		}
	}
}

