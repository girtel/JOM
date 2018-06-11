package by.examples;

import by.bsu.JVmipshell.*;
import java.io.IOException;

import java.util.Collections;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Multlotsize extends MIPshell
{
	Task task;
	VarVector s, x, y;

	int l(int t) {return task.getPeriodLength(t);}
	double s0(int j) {return task.getinitStock(j);}
	double sf(int j) {return task.getfinalStock(j);}
	double u(int j) {return task.getcapacity(j);}
	double c(int j, int t) {return task.getcost(j,t);}
	double h(int j, int t) {return task.getholdingCost(j,t);}
	double d(int j, int t) {return task.getdemand(j,t);}
	double tau(int i, int j) {return task.getprocTime(j,i);}
	Set<Integer> M(int j) {return task.M(j);}
	double rho(int j, int k) {return task.getResUsed(j,k);}
	Set<Integer> R(int j) {return task.R(j);}
	int m(int i, int t) {return task.getMachineNum(i,t);}
	double f(int i, int t) {return task.getMachineCost(i,t);}

	public Multlotsize(String name)
	{
		super(name);
		task = null;
		s = null; x = null; y = null;
	}

	public Multlotsize(String name, Task task)
	{
		super(name);
		model(task);
	} // end of Multlotsize(String name, Task task)

	void model(Task task)
	{
		this.task =task;
		int T = task.getHorizonLength();
		int m = task.getMachineNum();
		int n = task.getProductNum();

		VarVector s = new VarVector(this,"s",0,n,T);
		VarVector x = new VarVector(this,"x",VAR_INT,n,T);
		VarVector y = new VarVector(this,"y",VAR_INT,m,T);
		this.s = s; this.x = x; this.y =y;

		LinSum lsum = new LinSum();
		for (int t=0; t < T; ++t) {
			for (int j=0; j < n; ++j)
				lsum.add(c(j,t),x.get(j,t))
				    .add(h(j,t),s.get(j,t)); 
			for (int i=0; i < m; ++i)
				lsum.add(f(i,t),y.get(i,t));
		}
		minimize(lsum);
	  
		for (int j=0; j < n; ++j) {
			lsum = new LinSum();
			lsum.add(1.0,x.get(j,0))
			    .add(-1.0,s.get(j,0));
			for (int k=0; k < n; ++k) {
				if (R(k).contains(j))
					lsum.add(-rho(k,j),x.get(k,0));
			}
			addCtr(lsum,d(j,0)-s0(j),d(j,0)-s0(j));
		} // for (int j=0; j < n; ++j)
		
		for (int t=1; t < T; ++t) {
			for (int j=0; j < n; ++j) {
				lsum = new LinSum();
				lsum.add(1.0,s.get(j,t-1))
				    .add(1.0,x.get(j,t))
				    .add(-1.0,s.get(j,t)); // s(j,t-1) + x(j,t) - s(j,t)
				for (int k=0; k < n; ++k) {
					if (R(k).contains(j))
						lsum.add(-rho(k,j),x.get(k,t));
				}
				addCtr(lsum,d(j,t),d(j,t));
			}
		} // for (int t=1; t < T; ++t)
	
		for (int i=0; i < m; ++i) {
			for (int t=0; t < T; ++t) {
				y.get(i,t).setub(m(i,t));
				lsum = new LinSum();
				for (int j=0; j < n; ++j) {
					if (M(j).contains(i))
						lsum.add(tau(i,j),x.get(j,t));
				}
				addCtr(lsum.add(-l(t),y.get(i,t)),-VAR_INF,0.0);
			}
		} // for (int i=0; i < m; ++i)
		
		for (int j=0; j < n; ++j) {
			for (int t=0; t < T-1; ++t) {
				s.get(j,t).setub(u(j));
			}
			s.get(j,T-1).setbounds(sf(j),sf(j));
		}
	} // end of model()

     /**
      * @param[in] c string to be repeated;
      * @param[in] count number of time input string be repeated.
      * @return string that is composed of `count` copies of input string `c`.
      */
	static String repeatChar(String c, int count)
	{
		return String.join("", Collections.nCopies(count, c));
	}

     /**
      * The procedure prints a table filling it with data stored in a given `VarVector` object.
      * @param[in] tblName table name;
      * @param[in] rowTitle title printed in head row;
      * @param[in] colTitle title of first column;
      * @param[in] T horizon length;
      * @param[in] n (n,T) is size of 2-dimensional vector `x`;
      * @param[in] x reference to 2-dimensional `VarVector` object from `{this.x,this.s,this.y}`.      
      */
	void printTbl(String tblName, String rowTitle, String colTitle, int T, int n, VarVector x)
	{
		System.out.println("    " + tblName + ":");
		System.out.println(" " + repeatChar("_",10*n+8));
		int spl=(10*n-1 - rowTitle.length())/2;
		int spr=(10*n-1) - spl -rowTitle.length();
		System.out.println("|        |" + repeatChar(" ",spl) + rowTitle + repeatChar(" ",spr) + "|");
                System.out.println("|--------+" + repeatChar("-",10*n-1) + "|");
		spl=(6-colTitle.length())/2;
		spr=6-spl-colTitle.length();
		System.out.print("| " + repeatChar(" ",spl) + colTitle + repeatChar(" ",spr) + " |");
		for (int j=0; j < n; ++j) {
			System.out.printf("%8d |",j);
		}
        	System.out.printf("%n|--------+");
		for (int j=0; j < n-1; ++j) {
			System.out.printf("---------+");
		}
		System.out.printf("---------|");
		for (int t=0; t < T; ++t) {
			System.out.printf("%n| %6d |",(t+1));
			for (int j=0; j < n; ++j) {
				System.out.printf("%8.2f |",x.get(j,t).getval());
			}
		}
                System.out.printf("%n " + repeatChar("-",10*n+8) + "%n");
	}
 
	@Override
	public void printSolution()
	{
		if (isSolution()) {
			String optSol = (isSolutionOptimal())? "optimality proven": "optimality not proven";
			int T = task.getHorizonLength();
			int m = task.getMachineNum();
			int n = task.getProductNum();
        		System.out.printf("Objective value = %.4f - %s%n",getobjVal(),optSol);
			printTbl("Production","Products","Period",T,n,x);
			printTbl("Stock","Products","Period",T,n,s);
			printTbl("Machines","Machines","Period",T,m,y);
		}
		else {
			if (isInfeasible())
				System.out.println("Problem is infeasible");
			else				
				System.out.println("Time limit reached, no solution has been found");
		}
	}	
} // end of Multlotsize

