package by.examples;

import java.io.IOException;
import java.io.Reader;
import java.io.FileReader;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

/**
 * This class is designed to solve scheduling problems
 * written in text files in json format.
 */
public class testSchedule {
	static {
		System.loadLibrary("JVmipcl");
	}

	public static void main(String[] args) {
		try {
			Job[] jobs = readJobs(args[0]);
			Schedule sh = new Schedule("test");
			sh.model(jobs);
			sh.optimize();
			sh.printSolution();		
		} catch(Throwable e) {
			e.printStackTrace();
		}
	}
	
	static Job[] readJobs(String fileName) throws IOException {
		try(Reader reader = new FileReader(fileName)){
			Gson gson = new GsonBuilder().create();
			Job[] jobs = gson.fromJson(reader, Job[].class);
			return jobs;
		}
	}

}
