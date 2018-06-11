package by.examples;

public class testFcnf
{
	static {
            System.loadLibrary("JVmipcl");
        }

   	public static void main(String[] args) {
		if (args.length < 1) {
			System.out.println("File name is omitted!\n");
			return;
		}
		String path = args[0];
		int i1 = path.lastIndexOf('/') + 1;
		int i2 = path.lastIndexOf('.');
		if (i2 < 0)
			i2 = path.length();
		String name = path.substring(i1,i2);
		try {
			Fcnf prob = new Fcnf(name);
			prob.readNet(path);
			prob.buildProb();
			prob.optimize();
			prob.printSolution(name+".sol");
			prob.dispose();
		}
		catch(Throwable e) {
			System.out.println(e);
		}
	}
}
