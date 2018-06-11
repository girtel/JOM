package by.examples;

import by.bsu.JVmipshell.*;
import java.io.IOException;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

public class Schedule extends MIPshell
{
	static {
	        System.loadLibrary("JVmipcl");
	}

	static int min(int a, int b)
	{return (a < b)? a: b;}

	static int max(int a, int b)
	{return (a > b)? a: b;}
	
	/** Array of jobs to be scheduled */
	Job[] jobs;
	
	/** Number of jobs */
	int n;
	
	/** Number of machines */
	int m;

	/** Project start time */
	int Rmin;
	
	/** Project finish time */
	int Dmax;
	
	/** Variables that define the schedule */
	Map<Integer, Var> x;
	
	int r(int j) {return jobs[j].getrelease();}
	
	int d(int j) {return jobs[j].getdue();}
    
	int w(int j) {return jobs[j].getweight();}
    
	int p(int j, int i) {return jobs[j].getmachines().get(i);}
    		
	Set<Integer> M(int j) {return jobs[j].getmachines().keySet();}
    
	int[] prec(int j) {return jobs[j].getprec();}
    
	void processData()
	{
		int n = jobs.length;
		int m = 0, Rmin = 10000000, Dmax = 0;
		for (int j=0; j < n; ++j) {
			if (Rmin > r(j))
				Rmin = r(j);
			if (Dmax < d(j))
				Dmax = d(j);
			
			for (int i : M(j))
				if (m < i)
					m = i;
		}

		this.n = n;
		this.m = m+1;
		this.Rmin = Rmin;
		this.Dmax = Dmax;
	}
    
	int ind(int j, int i, int t)
	{return j + (i + t*m)*n;}

	public Schedule(String name)
	{
		super(name);
		jobs = null;
		x = null;
	}

	public Schedule(String name, Job[] jobs)
	{
		super(name);
		model(jobs);
	} //

	public void model(Job[] jobs)
	{
		this.jobs = jobs;
		processData();
		
		VarVector s = new VarVector(this,"s",0,n);
		VarVector pt = new VarVector(this,"pt",0,n);
		VarVector y = new VarVector(this,"y",VAR_BIN,n);
        
		Map<Integer, Var> x = new HashMap<Integer,Var>();
		this.x = x; // remember x to display the result later
		for (int j=0; j < n; ++j) {
			for (int i : M(j)) {
				for (int t=r(j); t <= d(j)-p(j,i); ++t) {
					x.put(ind(j,i,t), new Var(this,"x("+j+","+i+","+t+")",VAR_BIN));
				}
			}
		} // for (int j=0; j < n; ++j)
		
		LinSum lsum = new LinSum();
		
		for (int j=0; j < n; ++j) {
			lsum.add(w(j),y.get(j));
		}
		maximize(lsum);
		
		for (int i=0; i < m; ++i) {
			for (int t=Rmin; t < Dmax; ++t) {
				lsum = new LinSum();
				for (int j=0; j < n; ++j) {
					if (M(j).contains(i)) {
						int t1 = max(t-p(j,i),r(j));
						int t2 = min(t,d(j)-p(j,i));
						for (int tau=t1; tau <= t2; ++tau)
							lsum.add(1.0,x.get(ind(j,i,tau)));
					}
				}
				addCtr(lsum,-VAR_INF,1.0);
			}
		} // for (int i=0; i < m; ++i)
		
		for (int j=0; j < n; ++j) {
			lsum = new LinSum();
			LinSum lsum1 = new LinSum();
			LinSum lsum2 = new LinSum();
			for (int i : M(j)) {
				for (int t=r(j); t <= d(j)-p(j,i); ++t) {
					Var x0 = x.get(ind(j,i,t));
					lsum.add(1.0,x0);
					lsum1.add(t,x0);
					lsum2.add(p(j,i),x0);
				}
			}
			addCtr(lsum.add(-1.0,y.get(j)),0.0,0.0);
			addCtr(lsum1.add(-1.0,s.get(j)),0.0,0.0);
			addCtr(lsum2.add(-1.0,pt.get(j)),0.0,0.0);
            		
			if (prec(j) != null) {
				for (int j1 : prec(j)) { // inequalities for precedence relations
					lsum = new LinSum();
					addCtr(lsum.add(1.0,y.get(j)).add(-1.0,y.get(j1)),-VAR_INF,0.0);
					lsum = new LinSum();
					addCtr(lsum.add(1.0,y.get(j)).add(-1.0,y.get(j1))
							.add(-1.0,pt.get(j1)),0.0,VAR_INF);
				}
			}
		} // for (int j=0; j < n; ++j)

	} // end of model()
	
	@Override
	public void printSolution()
	{
		if (isSolution()) {
			String optSol = (isSolutionOptimal())? "optimality proven": "optimality not proven";
        		System.out.printf("Cost of scheduled jobs = %.2f - %s%n",getobjVal(),optSol);
			System.out.println(" ------------------------------------");
			System.out.println("|  Job | Processor | Starts |   Ends |");
			System.out.println("+------------------------------------+");
			for (int j=0; j < n; ++j) {
				for (int i : M(j)) {
					for (int t=r(j); t <= d(j)-p(j,i); ++t) {
						if (x.get(ind(j,i,t)).getval() > 0.5) {
							System.out.format("| %4d | %9d | %6d | %6d |%n",j,i,t,t+p(j,i));
        	                break;
						}
					}
				}
			}
			System.out.println(" ------------------------------------");
		}
		else {
			if (isInfeasible())
				System.out.println("Problem is infeasible");
			else				
				System.out.println("Time limit reached, no solution has been found");
		}
	}
} // end of Schedule

