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
public class testMultlotsize {
	static {
		System.loadLibrary("JVmipcl");
	}

	public static void main(String[] args) {
		try {
			Task task = readTask(args[0]);
			Multlotsize ml = new Multlotsize("test");
			ml.model(task);
			ml.optimize();
			ml.printSolution();		
		} catch(Throwable e) {
			e.printStackTrace();
		}
	}
	
	static Task readTask(String fileName) throws IOException {
		try(Reader reader = new FileReader(fileName)){
			System.out.println(fileName);
			Gson gson = new GsonBuilder().create();
			Task task = gson.fromJson(reader, Task.class);
			return task;
		}
	}

} // end of testMultlotsize

