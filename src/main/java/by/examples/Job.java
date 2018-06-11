package by.examples;

import java.util.Map;

public class Job {
	/** Job name */
	String name;

	/** Release date */
	int release;

	int getrelease()
	{return release;}

	/** Due date */
	int due;

	int getdue()
	{return due;}

	/** Weight oh this job */
	int weight;

	int getweight()
	{return weight;}

	/**
	 *  Keys are machines capable of doing this job;
	 *  the value of `machines.get(i)` is processing time of this job on machine `i`.
	 */
	Map<Integer,Integer> machines;

	Map<Integer,Integer> getmachines()
	{return machines;}

	/** Array of preceding jobs */
	int[] prec;

	int[] getprec() {return prec;}
}
