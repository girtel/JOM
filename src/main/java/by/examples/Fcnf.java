package by.examples;

import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;

import by.bsu.JVmipcl.*;

public class Fcnf extends MIP
{
	int vertNum, edgeNum;
	int[] demand;
	int[] tail, head, capacity, fixedCost, cost;

	public Fcnf(String name)
	{
		super(name);
		vertNum = edgeNum = 0;
	} // end of Fcnf(String name)

    	public void readNet(String filename) throws FileNotFoundException {
		Scanner scan = new Scanner(new FileReader(filename));
		vertNum = scan.nextInt();
		edgeNum = scan.nextInt();
		demand = new int[vertNum];
		tail = new int[edgeNum];
		head = new int[edgeNum];
		capacity = new int[edgeNum];
		fixedCost = new int[edgeNum];
		cost = new int[edgeNum];
		for (int v=0; v < vertNum; ++v) {
			demand[v] = scan.nextInt();
		}
		for (int e=0; e < edgeNum; ++e) {
			tail[e] = scan.nextInt();
			head[e] = scan.nextInt();
			capacity[e] = scan.nextInt();
			fixedCost[e] = scan.nextInt();
			cost[e] = scan.nextInt();
		}
		scan.close();
	} // end of readNet()

	public void buildProb()
	{
		final double[] val = {-1.0, 1.0, 1.0};
		double[] w = new double[1];
		int[] row = new int[3];
		int[] r = new int[1];
		int ctrNum=vertNum+edgeNum;

		openMatrix(vertNum+edgeNum,2*edgeNum,4*edgeNum);

		for (int v=0; v < vertNum; ++v) {
			addCtr(v,0,demand[v],demand[v]);
		}

		for (int e=vertNum; e < ctrNum; ++e) {
			addCtr(e,0,-INF,0.0);
		}

		for (int e=0; e < edgeNum; ++e) {
			row[0] = tail[e]; row[1] = head[e];
			row[2] = r[0] = vertNum+e;
			addColumn(e,0,-cost[e],0.0,VAR_INF,val,row,true);
			w[0] = -capacity[e];
			addColumn(e+edgeNum,VAR_BIN,-fixedCost[e],0.0,1.0,w,r,true);
		}
		closeMatrix();
	} // end of buildProb()

	@Override
	public void printSolution(String filename) throws IOException
	{
		BufferedWriter bfWriter =
			new BufferedWriter(new FileWriter(filename));
		if (isSolution()) {
			int varNum = 2*edgeNum;
			getSolution();
			double[] X = getX();
			int[] hd = getVarHandles();
			bfWriter.write("Nonzero flows:\n");
			for (int i=0; i < varNum; ++i) {
				if (hd[i] < edgeNum && X[i] > 0.5) {
					bfWriter.write("f(" + tail[hd[i]] + "," + head[hd[i]]
						 + ")=" + X[i] + "\n");
				}
			}
		}
		else
			bfWriter.write("Problem has no solution!\n");
		bfWriter.close();
	} // end of printSolution()
}
