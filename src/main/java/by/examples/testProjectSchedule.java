package by.examples;

import java.io.IOException;
import java.io.Reader;
import java.io.FileReader;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;

/**
 * This class is designed to solve project scheduling problems
 * written in text files in json format.
 */
public class testProjectSchedule {
	public static void main(String[] args) {
		try {
			Project prj = readProject(args[0]);
			ProjectSchedule sh = new ProjectSchedule(prj);
			sh.model();
			sh.optimize();
			sh.getSchedule();
			sh.printSchedule();		
		} catch(Throwable e) {
			e.printStackTrace();
		}
	}
	
	static Project readProject(String fileName) throws IOException {
		Project prj = null;
		Reader reader = null;
		try {
			reader = new FileReader(fileName);
			Gson gson = new GsonBuilder().create();
			prj = gson.fromJson(reader, Project.class);
		} finally {
			if (reader != null)
				reader.close();
		}
		return prj;
	}
}
