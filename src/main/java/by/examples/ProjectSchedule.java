package by.examples;

import by.bsu.JVmipshell.*;
import java.util.List;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Collections;

public class ProjectSchedule extends MIPshell
{
	static {
	        System.loadLibrary("JVmipcl");
	}
	
    /** Project to be scheduled */
	Project prj;

    /**
     * Length of the planning horizon
     */
	int H;
	
    /** Number of jobs */
	int n;
	
    /** Number of renewable resources */
	int renResNum;

    /** Number of non-renewable resources */
	int nonRenResNum;

    /** Optimal value of makespan */
	int makeSpan;

    /**
     * In optimization phase, `ef[j]` and `lf[j]` are, resp., early and late finish time of job `j`.
     *
     * When a  solution has been found, `getSchedule()` sets `ef[j]` to be the finish time of job `j`.
     */
	int[] ef, lf;

    /**
     * When a solution has been found, `getSchedule()` sets `mode[j]` to be the processing mode  of job `j`.
     * It should be noted that `mode` and `lf` are references to the same object. 
     */
	int[] mode;

    /**
     * This object represents _x_-variables in the MIP formulation of the problem.
     */
	List<List<Map<Integer,Var>>> x;
    
    /**
     * The constructor.
     *
     * @param[in] prj 
     */
	public ProjectSchedule(Project prj)
	{
		super(prj.getName());
		this.prj = prj;
		H = prj.getHorizon();
		n = prj.getJobNum();
		renResNum = prj.getRenResNum();
		nonRenResNum = prj.getNonRenResNum();
		ef = lf = mode = null;
		x = null;
	} // end of ProjectSchedule()

    /**
     * @param[in] a,b two integers.
     * @return minimum of `a` and `b`. 
     */
	static int min(int a, int b)
	{return (a < b)? a: b;}

    /**
     * @param[in] a,b two integers.
     * @return maximum of `a` and `b`. 
     */
	static int max(int a, int b)
	{return (a > b)? a: b;}

    /**
     * @param[in] j job index;
     * @param[in] m mode index.
     * @return processing time of job `j` in mode `m`.
     */
	int p(int j, int m) {return prj.getJobDuration(j,m);}
    
    /**
     * @param[in] i index of renewable resource.
     * @return number of units of renewable resource `i` available.
     */
	int R_r(int i) {return prj.getRenResQuantity(i);}

    /**
     * @param[in] j job index;
     * @param[in] m mode index;
     * @param[in] i index of renewable resource.
     * @return number of units of renewable resource `i` used when processing job `j` in mode `m`.
     */
	int rho_r(int j, int m, int i) {return prj.getRenResForJob(j,m,i);}

    /**
     * @param[in] i index of non-renewable resource.
     * @return number of units of non-renewable resource `i` available.
     */
	int R_n(int i) {return prj.getNonRenResQuantity(i);}

    /**
     * @param[in] j job index;
     * @param[in] m mode index;
     * @param[in] i index of non-renewable resource.
     * @return number of units of non-renewable resource `i` used when processing job `j` in mode `m`.
     */
	int rho_n(int j, int m, int i) {return prj.getNonRenResForJob(j,m,i);}

    /**
     * @param[in] j job index. 
     * @return number of processing modes for job `j`.
     */
	int M(int j) {return prj.getJobModeNum(j);}

    /**
     * @param[in] j job index.
     * @return array (list) of preceding jobs for job `j`.
     */
	int[] prec(int j) {return prj.getPrecJobs(j);}

    /** This procedure implements the critical path method.
     *
     * @param[in] H duration of time horizon;
     * @param[out] ef array of size `n`, `ef[j]` is early finish time of job `j`;
     * @param[out] lf array of size `n`, `ej[j]` is late finish time of job `j`.
     */
	private void cpm(int H)
	{
		ef = new int[n];
		lf = new int[n];

		for (int j=0; j < n; ++j) { // compute job early start times
            		int es = 0;
			if (prec(j) != null) {
            			for (int j1 : prec(j)) {
                			if (es < ef[j1])
                    				es = ef[j1];
				}
			}
			int minDur = H;
			for (int m=0; m < M(j); ++m) {
				int d = p(j,m);
				if (minDur > d)
					minDur = d;
			}
			ef[j] = es + minDur;
			lf[j] = H;
		}

		for (int j=n; --j >= 0; ) {
			int minDur = H;
			for (int m=0; m < M(j); ++m) {
				int d = p(j,m);
				if (minDur > d)
					minDur = d;
			}
			int ls = lf[j] - minDur;
			if (prec(j) != null) {
            			for (int j1 : prec(j)) {
                			if (lf[j1] > ls)
                    				lf[j1] = ls;
				}
			}
		}
	} // end of cpm()

	public void model()
	{
		cpm(H);
 
		Var T = new Var(this,"T",VAR_INT);
		VarVector e = new VarVector(this,"e",VAR_INT,n);
		VarVector d = new VarVector(this,"d",VAR_INT,n);
		
		x = new ArrayList<List<Map<Integer,Var>>>(n);
		for (int j=0; j < n; ++j) {
			x.add(new ArrayList<Map<Integer,Var>>(M(j)));
            		for (int m=0; m < M(j); ++m) {
				x.get(j).add(new HashMap<Integer,Var>());
				for (int t=ef[j]; t <= lf[j]; ++t) {
					x.get(j).get(m).put(t,new Var(this,"x("+j+","+m+","+t+")",VAR_BIN));
				}
			}
		}

		LinSum lsum, lsum1, lsum2;

		minimize(T);

		for (int j=0; j < n; ++j) {
			lsum = new LinSum();
			lsum1 = new LinSum();
			lsum2 = new LinSum();
			for (int m=0; m < M(j); ++m) {
				for (int t=ef[j]; t <= lf[j]; ++t) {
					lsum.add(1.0,x.get(j).get(m).get(t));
					lsum1.add(p(j,m),x.get(j).get(m).get(t));
					lsum2.add(t,x.get(j).get(m).get(t));
				}
			}
			addCtr(lsum,1.0,1.0); // each job ends only once and is processd in only one mode
			addCtr(lsum1.add(-1.0,d.get(j)),0.0,0.0); // sets duration of job j 
			addCtr(lsum2.add(-1.0,e.get(j)),0.0,0.0); // sets completion time of job j
			
			lsum = new LinSum();
			addCtr(lsum.add(1.0,e.get(j)).add(-1.0,T),-INF,0.0); // sets makespan: T >= e[j]
		} // for (int j=0; j < n; ++j)

		// renewable resource limits
		for (int tau=1; tau <= H; ++tau) {
			for (int i=0; i < renResNum; ++i) {
				lsum = new LinSum();
				for (int j=0; j < n; ++j) {
					for (int m=0; m < M(j); ++m) {
						for (int t=max(tau,ef[j]); t <= min(tau+p(j,m),lf[j]); ++t) {
							lsum.add(rho_r(j,m,i),x.get(j).get(m).get(t));
						}
					}
				}
				addCtr(lsum,-INF,R_r(i));
			}
		} // for (int tau=1; tau <= H; ++tau)

		// non-renewable resource limits
		for (int i=0; i < nonRenResNum; ++i) {
			lsum = new LinSum();
			for (int j=0; j < n; ++j) {
				for (int m=0; m < M(j); ++m) {
					for (int t=ef[j]; t <= lf[j]; ++t) {
						lsum.add(rho_n(j,m,i),x.get(j).get(m).get(t));
					}
				}
			}
			addCtr(lsum,-INF,R_n(i));
		} // for (int i=0; i < nonRenResNum; ++i)

		// precedence relations
		for (int j2=0; j2 < n; ++j2) {
			if (prec(j2) != null) {
				for (int j1 : prec(j2)) {
					lsum = new LinSum();
					addCtr(lsum.add(1.0,e.get(j1)) // e[j1] <= e[j2] - d[j2]
					           .add(-1.0,e.get(j2))
						   .add(1.0,d.get(j2)),-INF,0.0);
				}
			}
			else { // j2 has not predecessors
				lsum = new LinSum();
				addCtr(lsum.add(1.0,e.get(j2)) // e[j2] >= d[j2]
				           .add(-1.0,d.get(j2)),0.0,INF);
			}
		}
	} // end of model()

    /**
     * This procedure "extracts" from the values of `x`-variables an optimal schedule,
     * and strores it in two arrays, `ef` and `mode`:
     * job `j` is processed in mode `mode[j]`, it starts at `ef[j]-p(j,mode(j))` and ends at `ef[j]`.  
     */
	public void getSchedule()
	{
		if (isSolution()) {
			mode = lf;
			makeSpan = (int)(getobjVal() + 0.5);
			for (int j=0; j < n; ++j) {
				boolean stop = false;
				for (int m=0; m < M(j); ++m) {
					for (int t=ef[j]; t <= lf[j]; ++t) {
						if (x.get(j).get(m).get(t).getval() > 0.5) {
							mode[j] = m;
							ef[j] = t;
							stop = true;
							break;
						}
					}
					if (stop == true)
						break;
				}
			}
		}
	} // end of getSchedule()

    /**
     * @param[in] i renewable resource index;
     * @param[in] t time period that starts at time `t-1` and ends at time `t`.
     * @return amount of resource `i` used in period `t`. 
     */
	int getNonRenResUse(int i, int t)
	{
		int q=0;
		for (int j=0; j < n; ++j) {
			if ( t > ef[j]-p(j,mode[j]) && t <= ef[j])
				q += rho_r(j,mode[j],i);
		}
		return q;
	} // end of getNonRenResUse()

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
      * The procedure prints the schedule found by the solver.
      */
	public void printSchedule()
	{
		if (isSolution()) {
			String optSol = (isSolutionOptimal())? "optimality proven": "optimality not proven";
        		System.out.printf("Makespan = %6d - %s%n",makeSpan,optSol);
			System.out.println("Job schedule:");
			System.out.println(" _______________________________");
			System.out.println("|  Job | Mode |  Start |    End |");
			System.out.println("|------+------+--------+--------|");
			for (int j=0; j < n; ++j) {
				System.out.printf("| %4d | %4d | %6d | %6d |%n",j,mode[j],ef[j]-p(j,mode[j]),ef[j]);
			}
			System.out.println(" -------------------------------");

			System.out.printf("%nResource usage:%n");
			System.out.println(" _____" + repeatChar("_",9*renResNum));
			System.out.print("| Period |");
			for (int i=0; i < renResNum; ++i) {
				System.out.printf(" %5d |",i);
			}
			System.out.printf("%n|--------");
			for (int i=0; i < renResNum; ++i) {
				System.out.print("+-------");
			}
			System.out.println("|");
			for (int t=1; t <= makeSpan; ++t) {
				System.out.printf("| %6d |",t);
				for (int i=0; i < renResNum; ++i) {
					System.out.printf(" %5d |",getNonRenResUse(i,t));
				}
				System.out.println("");
			}
			System.out.println(" -----" + repeatChar("-",9*renResNum));
		}
		else {
			if (isInfeasible())
				System.out.println("Problem is infeasible");
			else				
				System.out.println("Time limit reached, no solution has been found");
		}
	} // end of printSchedule() 

} // end of ProjectSchedule

