package by.examples;

import java.util.Map;

class Resource {
	/** Resource name */
	String name;

	/** Available quantity */
	int quantity;
} // end of Resource

class Mode {
	/** Processing time of Job in this mode */
	int duration;

    /**
     *  Keys are resources needed for this job;
     *  the value of `renRes.get(r)` is amount of resource `r` needed at any time during the execution of this job.
     */
	Map<Integer,Integer> renRes;

    /**
     *  Keys are resources needed for this job;
     *  the value of `nonRenRes.get(r)` is amount of resource `r` needed to fulfill this job.
     */
	Map<Integer,Integer> nonRenRes;
} // end of Mode

class Job {
    /** Job name */
	String name;

    /**
     * Array of this job processing modes
     */
	Mode[] modes;

    /** Array of preceding jobs */
	int[] prec;
} // end of Job

public class Project {
    /** Project name */
	String name;

    /**
     * @return project name.
     */
	String getName()
	{return name;}

    /**
     * Length of the planning horizon
     */
	int horizon;

    /**
     * @return length of the planning horizon.
     */
	int getHorizon()
	{return horizon;}

    /** Renewable resources.
     *
     * The length of this array is the number of renewable resources used by the project.
     */
	Resource[] renRes;

    /**
     * @return number of renewable resources.
     */
	int getRenResNum()
	{return (renRes != null)? renRes.length: 0;}

    /**
     * @param[in] i index of a renewable resource.
     * @return name of resource `i`.
     */
	String getRenResName(int i)
	{return renRes[i].name;}

    /**
     * @param[in] i index of a renewable resource.
     * @return available quantity of renewable resource `i`.
     */
	int getRenResQuantity(int i)
	{return renRes[i].quantity;}

    /** Non-renewable resources.
     *
     * The length of this array is the number of non-renewable resources used by the project.
     */
	Resource[] nonRenRes;

    /**
     * @return number of non-renewable resources.
     */
	int getNonRenResNum()
	{return (nonRenRes != null)? nonRenRes.length: 0;}

    /**
     * @param[in] i index of a non-renewable resource.
     * @return name of resource `i`.
     */
	String getNonRenResName(int i)
	{return nonRenRes[i].name;}

    /**
     * @param[in] i index of a non-renewable resource.
     * @return available quantity of non-renewable resource `i`.
     */
	int getNonRenResQuantity(int i)
	{return nonRenRes[i].quantity;}

    /** Array of jobs */
	Job[] jobs;

    /**
     * @return number of jobs in the project.
     */
	int getJobNum()
	{return (jobs != null)? jobs.length: 0;}

    /**
     * @param[in] j job index.
     * @return name of job indexed by `j`.
     */
	String getJobName(int j)
	{return jobs[j].name;}

    /**
     * @param[in] j job index.
     * @return number of modes available for processing job `j`.
     */
	int getJobModeNum(int j)
	{return jobs[j].modes.length;}

    /**
     * @param[in] j job index;
     * @param[in] m processing mode.
     * @return job `j` processing time in mode `m`.
     */
	int getJobDuration(int j, int m)
	{return jobs[j].modes[m].duration;}

    /**
     * @param[in] j job index;
     * @param[in] m processing mode;
     * @param[in] i renewable resource. 
     * @return number of units of renewable resource `i` used by job `j` in mode `m`.
     */
	int getRenResForJob(int j, int m, int i)
	{
		int q = 0;
		if (jobs[j].modes[m].renRes != null) {
			Integer Q = jobs[j].modes[m].renRes.get(i);
			q = (Q != null)? Q: 0;
		}
		return q;
	}

    /**
     * @param[in] j job index;
     * @param[in] m processing mode;
     * @param[in] i non-renewable resource.
     * @return number of units of non-renewable resource `i` used by job `j` in mode `m`.
     */
	int getNonRenResForJob(int j, int m, int i)
	{
		int q = 0;
		if (jobs[j].modes[m].nonRenRes != null) {
			Integer Q = jobs[j].modes[m].nonRenRes.get(i);
			q = (Q != null)? Q: 0;
		}
		return q;
	}

    /**
     * @param[in] j job index.
     * @return array (list) of preceding jobs.
     */
	int[] getPrecJobs(int j)
	{return jobs[j].prec;}
} // end of Project

