/*******************************************************************************
 * Copyright (c) 2015 Pablo Pavon Mari�o.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * <p>
 * Contributors:
 * Pablo Pavon Mari�o - initial API and implementation
 ******************************************************************************/

package com.jom;

import cern.colt.function.tdouble.DoubleDoubleFunction;
import cern.colt.function.tdouble.DoubleFunction;
import cern.colt.function.tdouble.DoubleProcedure;
import cern.colt.list.tdouble.DoubleArrayList;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.*;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntFactory2D;
import cern.colt.matrix.tint.IntMatrix1D;
import cern.colt.matrix.tint.IntMatrix2D;
import cern.jet.math.tdouble.DoublePlusMultSecond;
import cern.jet.math.tint.IntFunctions;

/** This class implements the functionality to handle an array of doubles of an arbitrary number of dimensions. Internally, data is stored in a 1D
 * vector of the type DoubleMatrix1D, in the COLT library. The n-dim array can be dense or sparse, with NO DIFFERENCE AT ALL in the interface with
 * the user. Internally, the data is stored in a sparse or dense 1D vector in each case, but all the operations have the same effect.
 *
 * The methods to manipulate and operate with the array follow a similar philosophy as the ones COLT library provides for handling 1D, 2D and 3D
 * arrays. See
 * http://sourceforge.net/projects/parallelcolt/ for more details
 *
 * @author Pablo Pavon Mariño
 * @see <a href="http://www.net2plan.com/jom">http://www.net2plan.com/jom</a>
 */
public class DoubleMatrixND
{
	private int            numDim;
	private int            numElem;
	private IntMatrix1D    size;
	private DoubleMatrix1D x;

	/** Creates a 1x1 dense array with the given value.
	 * @param value The initializing data */
	public DoubleMatrixND(double value)
	{
		this(new int[]{1, 1}, DoubleFactory1D.dense.make(1, value));
	}

	/** Creates a 2D dense array, with the same size as value, and initializes with the data in values
	 * @param values The initializing data */
	public DoubleMatrixND(double[][] values)
	{
		this(values, "dense");
	}

	/** Creates a 2D array of the appropriated type, with the same size as value, and initializes with the data in values
	 * @param values The initializing data
	 * @param type "sparse" or "dense" for creating sparse or dense arrays respectively */
	public DoubleMatrixND(double[][] values, String type)
	{
		this((type.equals("dense")) ? DoubleFactory2D.dense.make(values) : DoubleFactory2D.sparse.make(values));
	}

	/** Creates a 3D dense array, with the same size as value, and initializes with the data in values
	 * @param values The initializing data */
	public DoubleMatrixND(double[][][] values)
	{
		this(values, "dense");
	}

	/** Creates a 3D array of the appropriated type, with the same size as value, and initializes with the data in values
	 * @param values The initializing data
	 * @param type "sparse" or "dense" for creating sparse or dense arrays respectively */
	public DoubleMatrixND(double[][][] values, String type)
	{
		this(new int[]{values.length, values[0].length, values[0][0].length}, (type.equals("dense")) ?
				DoubleFactory1D.dense.make(flattenArray(values)) :
				DoubleFactory1D.sparse.make(flattenArray(values)));
	}

	/** Creates a 4D dense array, with the same size as value, and initializes with the data in values
	 * @param values The initializing data */
	public DoubleMatrixND(double[][][][] values)
	{
		this(values, "dense");
	}

	/** Creates a 4D array of the appropriated type, with the same size as value, and initializes with the data in values
	 * @param values The initializing data
	 * @param type "sparse" or "dense" for creating sparse or dense arrays respectively */
	public DoubleMatrixND(double[][][][] values, String type)
	{
		this(new int[]{values.length, values[0].length, values[0][0].length, values[0][0][0].length}, (type.equals("dense")) ?
				DoubleFactory1D.dense.make(flattenArray(values)) :
				DoubleFactory1D.sparse.make(flattenArray(values)));
	}

	/** Creates a 2D array of the same type (sparse, dense) as values, and also with the same size and data
	 * @param values The initializing data */
	public DoubleMatrixND(DoubleMatrix2D values)
	{
		this(new int[]{values.rows(), values.columns()}, values.copy().vectorize());
	}

	/** Creates a 3D array of the same type (sparse, dense) as values, and also with the same size and data
	 * @param values The initializing data */
	public DoubleMatrixND(DoubleMatrix3D values)
	{
		this(new int[]{values.slices(), values.rows(), values.columns()}, values.copy().vectorize());
	}

	/** Creates a dense empty array of the given size, initialized to zeros
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension. */
	public DoubleMatrixND(int[] size)
	{
		this(size, "dense");
	}

	/** Creates a dense array of the given size, initialized to value
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param value The initializing data */
	public DoubleMatrixND(int[] size, double value)
	{
		this(size, value, "dense");
	}

	/** Creates an array of the given size and type. All the elements of the array equal to value
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param value The initializing data
	 * @param type "sparse" or "dense" for creating sparse or dense arrays respectively */
	public DoubleMatrixND(int[] size, double value, String type)
	{
		this(size, (type.equals("dense")) ?
				DoubleFactory1D.dense.make(numElements(size), value) :
				DoubleFactory1D.sparse.make(numElements(size), value));
	}

	/** Creates a dense array of the given size, using values as initializing data. values must be a vector of as many elements as the number of
	 * cells in the array (product of the numbers in size parameter). i-th index in values vector is put in the array position of i-th linear index.
	 * See the JOM help to see how linear indexes correspond to positions in the array.
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param values The initializing data */
	public DoubleMatrixND(int[] size, double[] values)
	{
		this(size, values, "dense");
	}

	/** Creates an array of the given size and type, using values as initializing data. values must be a vector of as many elements as the number of
	 *  cells in the array (product of the numbers in size parameter). i-th index in values vector is put in the array position of i-th linear index
	 *  . See the JOM help to see how linear indexes correspond to positions in the array.
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param values The initializing data
	 * @param type "sparse" or "dense" for creating sparse or dense arrays respectively */
	public DoubleMatrixND(int[] size, double[] values, String type)
	{
		this(size, (type.equals("dense")) ? DoubleFactory1D.dense.make(values) : DoubleFactory1D.sparse.make(values));
	}

	/** Creates an array of the given size and the same type as values, using values as initializing data. values must have as many elements as the
	 * number of cells in the array (product of the numbers in size parameter). i-th index in values vector is put in the array position of i-th
	 * linear index. See the JOM help to see how linear indexes correspond to positions in the array.
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param values The initializing data */
	public DoubleMatrixND(int[] size, DoubleMatrix1D values)
	{
		long totalNumberCells = 1;
		for (int cont = 0; cont < size.length; cont++)
			totalNumberCells *= size[cont];
		if (totalNumberCells >= Integer.MAX_VALUE)
			throw new RuntimeException("The size of the array is limited to " + Integer.MAX_VALUE + " cells. This limitation comes from the COLT "
					+ "libraries");

		this.size = IntFactory1D.dense.make(size).copy();
		this.numDim = size.length;
		if (this.size.getMinLocation()[0] <= 0) throw new RuntimeException("Wrong size of n-array");
		this.numElem = IntMatrixND.prod(this.size);
		if (this.numElem != (int) values.size()) throw new RuntimeException("Wrong size of n-array");
		this.x = values;
		this.x.trimToSize();
	}

	/** Creates a dense array of the given size, using values as initializing data, previously casted to double. values must be a vector of as many
	 * elements as the number of cells in the array (product of the numbers in size parameter). i-th index in values vector is put in the array
	 * position of i-th linear index. See the JOM help to see how linear indexes correspond to positions in the array.
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param values The initializing data */
	public DoubleMatrixND(int[] size, int[] values)
	{
		this(size, IntFactory1D.dense.make(values));
	}

	/** Creates an array of the given size and type, using values as initializing data, previously casted to double. values must be a vector of as
	 * many elements as the number of cells in the array (product of the numbers in size parameter). i-th index in values vector is put in the array
	 * position of i-th linear index. See the JOM help to see how linear indexes correspond to positions in the array.
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param values The initializing data
	 * @param type "sparse" or "dense" for creating sparse or dense arrays respectively */
	public DoubleMatrixND(int[] size, int[] values, String type)
	{
		this(size, (type.equals("dense")) ? IntFactory1D.dense.make(values) : IntFactory1D.sparse.make(values));
	}

	/** Creates an empty array of the given size and type, initialized to zeros
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @param type "sparse" or "dense" for creating sparse or dense arrays respectively. "random" for creating a matrix with random numbers
	 * between 0 and 1 sampled from a uniform distribution  */
	public DoubleMatrixND(int[] size, String type)
	{
		this(size, (type.equals("random")) ?
				DoubleFactory1D.dense.random(numElements(size)) :
				(type.equals("dense")) ? DoubleFactory1D.dense.make(numElements(size), 0) : DoubleFactory1D.sparse.make(numElements(size), 0));
	}

	/* Just one internal accessible copy with IntMatrix1D in the size */
	DoubleMatrixND(IntMatrix1D size, DoubleMatrix1D values)
	{
		this(size.toArray(), values);
	}

	DoubleMatrixND(int[] size, IntMatrix1D values)
	{
		this(size, new IntMatrixND(size, values).cast2Double().elements());
	}

	/** For arrays of the given size, converts the position of a cell given as a linear index, into its associated coordinates in the arrays.
	 * @param index Linear index of the cell position
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @return The coordinates of the cell in the array */
	public static IntMatrix1D ind2sub(int index, int[] size)
	{
		IntMatrix2D res = DoubleMatrixND.ind2sub(new int[]{index}, size);
		return res.viewRow(0);
	}

	static IntMatrix1D ind2sub(int index, IntMatrix1D size)
	{
		return ind2sub(index, size.toArray());
	}

	/** For arrays of the given size, converts the positions of the cells given as a linear index in indexes, into their associated coordinates in
	 * the arrays.
	 * @param indexes Linear indexes of the cells to convert
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @return A 2D array, as many columns as dimensions in the array. i.th row of the array contains the coordinates of the i-th linear index to
	 * convert */
	public static IntMatrix2D ind2sub(int[] indexes, int[] size)
	{
		IntMatrix1D indexes1D = IntFactory1D.dense.make(indexes);
		int numCells = IntMatrixND.prod(size);
		int numDim = size.length;
		int numIndexes = indexes.length;
		if (indexes1D.size() > 0)
			if (indexes1D.getMaxLocation()[0] > numCells - 1){ throw new RuntimeException("Index out of bounds in n-dim array"); }
		if (indexes1D.size() > 0) if (indexes1D.getMinLocation()[0] < 0){ throw new RuntimeException("Index out of bounds in n-dim array"); }

		IntMatrix1D blockSizes = DoubleMatrixND.computeBlockSizes(size);
		IntMatrix2D res = IntFactory2D.dense.make((int) indexes1D.size(), numDim);
		IntMatrix1D indexesToReduce = indexes1D.copy();
		for (int contDim = numDim - 1; contDim >= 0; contDim--)
		{
			IntMatrix1D thisBlockSizeVector = IntFactory1D.dense.make(numIndexes, blockSizes.get(contDim));
			IntMatrix1D coord_n = indexesToReduce.copy().assign(thisBlockSizeVector, cern.jet.math.tint.IntFunctions.div);
			res.viewColumn(contDim).assign(coord_n);

			indexesToReduce.assign(thisBlockSizeVector, cern.jet.math.tint.IntFunctions.mod);
		}
		return res;
	}

	static IntMatrix2D ind2sub(int[] indexes, IntMatrix1D size)
	{
		return ind2sub(indexes, size.toArray());
	}

	/** For arrays of the given size, converts the position of the cell given as a their coordinates, into its associated linear index.
	 * @param coord Coordinates of the cell position
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @return The linear index associated to the cell position */
	public static int sub2ind(IntMatrix1D coord, int[] size)
	{
		IntMatrix1D blockSizes = DoubleMatrixND.computeBlockSizes(size);
		return coord.zDotProduct(blockSizes);
	}

	static int sub2ind(IntMatrix1D coord, IntMatrix1D size)
	{
		return sub2ind(coord, size.toArray());
	}

	/** For arrays of the given size, converts the positions of the cells given as a their coordinates, into their associated linear index.
	 * @param coord Coordinates of the cell position. Each row is a coordinate
	 * @param size Array size. The i-th element in size indicates the size of the array in the corresponding dimension.
	 * @return The set of linear indexes associated to the cell positions. i-th element corresponds to i-th row in coord */
	public static IntMatrix1D sub2ind(IntMatrix2D coord, int[] size)
	{
		IntMatrix1D blockSizes = DoubleMatrixND.computeBlockSizes(size);
		IntMatrix1D indexes = coord.zMult(blockSizes, null);
		return indexes;
	}

	static IntMatrix1D sub2ind(IntMatrix2D coord, IntMatrix1D size)
	{
		return sub2ind(coord, size.toArray());
	}

	static DoubleMatrix1D kron(DoubleMatrix1D a, DoubleMatrix1D b)
	{
		DoubleMatrix1D res = DoubleFactory1D.sparse.make((int) (a.size() * b.size()));
		IntArrayList indexList = new IntArrayList();
		DoubleArrayList valueList = new DoubleArrayList();
		a.getNonZeros(indexList, valueList);
		for (int cont = 0; cont < indexList.size(); cont++)
		{
			/* for each non-zero value in a, set a block in res = (value * b) */
			int index = indexList.get(cont);
			double value = valueList.get(cont);
			res.viewPart(index * (int) b.size(), (int) b.size()).assign(b.copy(), DoublePlusMultSecond.plusMult(value));
		}
		return res;
	}

	static DoubleMatrix2D kron(DoubleMatrix2D a, DoubleMatrix2D b)
	{
		long init = System.nanoTime();

		DoubleMatrix2D res = DoubleFactory2D.sparse.make(a.rows() * b.rows(), a.columns() * b.columns());
		IntArrayList rowList = new IntArrayList();
		IntArrayList columnList = new IntArrayList();
		DoubleArrayList valueList = new DoubleArrayList();
		a.getNonZeros(rowList, columnList, valueList);
		for (int cont = 0; cont < rowList.size(); cont++)
		{
			/* for each non-zero value in a, set a block in res = (value * b) */
			int row = rowList.get(cont);
			int column = columnList.get(cont);
			double value = valueList.get(cont);
			res.viewPart(row * b.rows(), column * b.columns(), b.rows(), b.columns()).assign(b.copy(), DoublePlusMultSecond.plusMult(value));
		}
		return res;
	}

	private static IntMatrix1D computeBlockSizes(int[] size)
	{
		IntMatrix1D blockSizes = IntFactory1D.dense.make(size.length);
		blockSizes.setQuick(0, 1);
		for (int dim = 1; dim < size.length; dim++)
			blockSizes.setQuick(dim, blockSizes.getQuick(dim - 1) * size[dim - 1]);
		return blockSizes;
	}

	private static double[] flattenArray(double[][][] a)
	{
		int s1 = a.length;
		int s2 = a[0].length;
		int s3 = a[0][0].length;
		double[] res = new double[s1 * s2 * s3];
		int cont = 0;
		for (int c3 = 0; c3 < s3; c3++)
			for (int c2 = 0; c2 < s2; c2++)
				for (int c1 = 0; c1 < s1; c1++)
					res[cont++] = a[c1][c2][c3];
		return res;
	}

	private static double[] flattenArray(double[][][][] a)
	{
		int s1 = a.length;
		int s2 = a[0].length;
		int s3 = a[0][0].length;
		int s4 = a[0][0][0].length;

		double[] res = new double[s1 * s2 * s3 * s4];
		int cont = 0;
		for (int c4 = 0; c4 < s4; c4++)
			for (int c3 = 0; c3 < s3; c3++)
				for (int c2 = 0; c2 < s2; c2++)
					for (int c1 = 0; c1 < s1; c1++)
						res[cont++] = a[c1][c2][c3][c4];
		return res;
	}

	private static int numElements(int[] size)
	{
		int res = 1;
		for (int cont = 0; cont < size.length; cont++)
			res = res * size[cont];
		return res;
	}

	/****************** AGGREGATE METHODS ********************************/

	/** Applies a function to each cell and aggregates the results. Applies a function to each cell and aggregates the results. Returns a value v
	 * such that v==a(size()) where a(i) == aggr( a(i-1), f(get(row,column)) ) and terminators are a(1) == f(get(0,0)), a(0)==Double.NaN.
	 * @param aggr an aggregation function taking as first argument the current aggregation and as second argument the transformed current cell
	 *                value.
	 * @param f a function transforming the current cell value.
	 * @return the aggregated measure */
	public double aggregate(DoubleDoubleFunction aggr, DoubleFunction f)
	{
		return x.aggregate(aggr, f);
	}

	/** Applies a function to all cells with the given indexes and aggregates the results
	 * @param aggr an aggregation function taking as first argument the current aggregation and as second argument the transformed current cell
	 *                value.
	 * @param f a function transforming the current cell value
	 * @param indexList indexes of the cells
	 * @return the aggregated measure */
	public double aggregate(DoubleDoubleFunction aggr, DoubleFunction f, IntArrayList indexList)
	{
		return x.aggregate(aggr, f, indexList);
	}

	/** Applies a function to each corresponding cell of two matrices and aggregates the results. Returns a value v such that v==a(size()) where a
	 * (i) == aggr( a(i-1), f(get(row,column),other.get(row,column)) ) and terminators are a(1) == f(get(0,0),other.get(0,0)), a(0)==Double.NaN.
	 * @param other the second array to make the aggregation
	 * @param aggr an aggregation function taking as first argument the current aggregation and as second argument the transformed current cell
	 *                values
	 * @param f a function transforming the current cell values
	 * @return the aggregated measure */
	public double aggregate(DoubleMatrixND other, DoubleDoubleFunction aggr, DoubleDoubleFunction f)
	{
		return x.aggregate(other.elements(), aggr, f);
	}

	/** Fills the current array with ascending numbers (in the order of the linear indexes)
	 * @param first First number, assigned to index 0 of the array
	 * @param step step between the values in two consecutive cells
	 * @return this */
	public DoubleMatrixND ascending(double first, double step)
	{
		for (int cont = 0; cont < this.numElem; cont++)
			this.x.set(cont, first + cont * step);
		return this;
	}

	/** Sets all cells to the state specified by value
	 * @param value the value to be filled into the cells
	 * @return this */
	public DoubleMatrixND assign(double value)
	{
		x.assign(value);
		return this;
	}

	/** Sets all cells to the state specified by values. values is required to have as many cells as this. The values are copied. So subsequent
	 * changes in values are not reflected in the matrix, and vice-versa.
	 * @param values the values to be filled into the cells.
	 * @return this */
	public DoubleMatrixND assign(double[] values)
	{
		this.x.assign(values);
		return this;
	}

	/** Assigns the result of a function to each cell. x[i] = function(x[i]).
	 * @param f a function object taking as argument the current cell's value.
	 * @return this */
	public DoubleMatrixND assign(DoubleFunction f)
	{
		this.x.assign(f);
		return this;
	}

	/** Replaces all cell values of the receiver with the values of another matrix. Both matrices must have the same size. If both matrices share
	 * the same cells (as is the case if they are views derived from the same matrix) and intersect in an ambiguous way, then replaces as if using
	 * an intermediate auxiliary deep copy of other.
	 * @param other the source matrix to copy from (may be identical to the receiver).
	 * @return this */
	public DoubleMatrixND assign(DoubleMatrixND other)
	{
		if (!this.size.equals(other.getSize())) throw new RuntimeException("Wrong size of n-array");
		this.x.assign(other.elements());
		return this;
	}

	/** Assigns the result of a function to each cell; x[i] = function(x[i],y[i]).
	 * @param y the secondary array to operate on.
	 * @param function a function object taking as first argument the current cell's value of this, and as second argument the current cell's value
	 *                    of y,
	 * @return this */
	public DoubleMatrixND assign(DoubleMatrixND y, DoubleDoubleFunction function)
	{
		this.x.assign(y.elements(), function);
		return this;
	}

	/** Assigns the result of a function to all cells with a given indexes function - a function object taking as first argument the current cell's
	 * value of this, and as second argument the current cell's value of y
	 * @param y the secondary array to operate on
	 * @param function a function object taking as first argument the current cell's value of this, and as second argument the current cell's value
	 *                    of y,
	 * @param indexList the list of indexes
	 * @return this */
	public DoubleMatrixND assign(DoubleMatrixND y, DoubleDoubleFunction function, IntArrayList indexList)
	{
		this.x.assign(y.elements(), function, indexList);
		return this;
	}

	/** Assigns a value to all cells that satisfy a condition.
	 * @param cond a condition
	 * @param value a value
	 * @return this */
	public DoubleMatrixND assign(DoubleProcedure cond, double value)
	{
		this.x.assign(cond, value);
		return this;
	}

	/** Assigns the result of a function to all cells that satisfy a condition.
	 * @param cond a condition
	 * @param f a function object
	 * @return this */
	public DoubleMatrixND assign(DoubleProcedure cond, DoubleFunction f)
	{
		this.x.assign(cond, f);
		return this;
	}

	/** Returns the number of cells having non-zero values; ignores tolerance
	 * @return cardinality */
	public int cardinality()
	{
		return x.cardinality();
	}

	/** Creates a new copy of the array, with all the elements transformed to integer
	 * @return the copy */
	IntMatrixND cast2Int()
	{
		IntMatrix1D res = IntFactory1D.dense.make(this.numElem);
		IntArrayList indexList = new IntArrayList();
		DoubleArrayList valueList = new DoubleArrayList();
		getNonZeros(indexList, valueList);
		indexList.trimToSize();
		int[] indexesNonZero = indexList.elements();
		for (int cont = 0; cont < indexesNonZero.length; cont++)
			res.set(indexesNonZero[cont], (int) this.x.get(indexesNonZero[cont]));
		return new IntMatrixND(this.size, res);
	}

	/** Constructs and returns a deep copy of the receiver. Note that the returned matrix is an independent deep copy. The returned matrix is not
	 * backed by this matrix, so changes in the returned matrix are not reflected in this matrix, and vice-versa.
	 * @return a deep copy of the receiver */
	public DoubleMatrixND copy()
	{
		return new DoubleMatrixND(this.size.toArray(), this.x.copy());
	}

	/** Fills the current array with descending numbers (in the order of the linear indexes)
	 * @param first First number, assigned to index 0 of the array
	 * @param step step between the values in two consecutive cells
	 * @return this */
	public DoubleMatrixND descending(double first, double step)
	{
		for (int cont = 0; cont < this.numElem; cont++)
			this.x.set(cont, first - cont * step);
		return this;
	}

	/** Returns the elements of this matrix
	 * @return the elements */
	public DoubleMatrix1D elements()
	{
		return x;
	}

	/** Returns whether all cells are equal to the given value
	 * @param value the value to test against.
	 * @return true if all cells are equal to the given value, false otherwise */
	public boolean equals(double value)
	{
		return x.equals(value);
	}

	/** Compares this object against the specified object. The result is true if and only if the argument is not null and is at least a
	 * DoubleMatrixND object that has the same size as the receiver and has exactly the same values at the same coordinates
	 * @param other the object to compare with
	 * @return true if the objects are the same; false otherwise */
	public boolean equals(DoubleMatrixND other)
	{
		return (this.x.equals(other.elements()) && (this.size.equals(other.getSize())));
	}

	/** Returns the matrix cell value at the given index.
	 * @param index the index of the cell
	 * @return the value of the specified cell */
	public double get(int index)
	{
		return x.get(index);
	}

	/** Returns the matrix cell value at the given subindexes
	 * @param subindexes the coordinates of the cell
	 * @return the value of the specified cell */
	public double get(int[] subindexes)
	{
		return x.get(DoubleMatrixND.sub2ind(IntFactory1D.dense.make(subindexes), this.size.toArray()));
	}

	/** Returns the matrix cell value at the given subindexes
	 * @param subindexes the coordinates of the cell
	 * @return the value of the specified cell */
	double get(IntMatrix1D subindexes)
	{
		return x.get(DoubleMatrixND.sub2ind(subindexes, this.size.toArray()));
	}

	/** Return maximum value of this matrix together with its location. If the array has no elements, an empty array is returned
	 * @return {maximum_value, cell index (cast to double)} */
	public double[] getMaxLocation()
	{
		return x.size() > 0 ? this.x.getMaxLocation() : new double[]{};
	}

	/** Return minimum value of this matrix together with its location. If the array has no elements, an empty array is returned
	 * @return {minimum_value, cell index (cast to double)} */
	public double[] getMinLocation()
	{
		return x.size() > 0 ? this.x.getMinLocation() : new double[]{};
	}

	/** Fills the indexes and values of cells having negative values into the specified lists. Fills into the lists, starting at index 0. After this
	 *  call returns the specified lists all have a new size, the number of non-zero values.
	 * @param indexList the list to be filled with the indexes, can have any size
	 * @param valueList the list to be filled with values, can have any size */
	public void getNegativeValues(IntArrayList indexList, DoubleArrayList valueList)
	{
		this.x.getNegativeValues(indexList, valueList);
	}

	/** Fills the indexes and values of cells having non zero values into the specified lists. Fills into the lists, starting at index 0. After this
	 *  call returns the specified lists all have a new size, the number of non-zero values. In general, fill order is unspecified.
	 * @param indexList the list to be filled with the indexes, can have any size
	 * @param valueList the list to be filled with values, can have any size */
	public void getNonZeros(IntArrayList indexList, DoubleArrayList valueList)
	{
		this.x.getNonZeros(indexList, valueList);
	}

	/** Gets the number of dimensions of the array
	 * @return the number of dimensions */
	public int getNumDim()
	{
		return this.numDim;
	}

	/** Gets the number of cells in the array (e.g for an array of size 2x3x5 the number of cells is 30)
	 * @return the number of cells */
	public int getNumElements()
	{
		return this.numElem;
	}

	/** Fills the indexes and values of cells having positive values into the specified lists. Fills into the lists, starting at index 0. After this
	 *  call returns the specified lists all have a new size, the number of non-zero values.
	 * @param indexList the list to be filled with the indexes, can have any size
	 * @param valueList the list to be filled with values, can have any size */
	public void getPositiveValues(IntArrayList indexList, DoubleArrayList valueList)
	{
		this.x.getPositiveValues(indexList, valueList);
	}

	/** Returns the matrix cell value at coordinate index. Provided with invalid parameters this method may return invalid objects without throwing
	 * any exception. You should only use this method when you are absolutely sure that the coordinate is within bounds. Precondition (unchecked):
	 * \(index {@literal <} 0 \parallel index \ge size()\).
	 * @param index the index of the cell
	 * @return the value of the specified cell */
	public double getQuick(int index)
	{
		return this.x.getQuick(index);
	}

	/** Gets the size of the array
	 * @return a copy of the vector with the array size. The i-th element in size indicates the size of the array in the corresponding dimension. */
	public IntMatrix1D getSize()
	{
		return this.size.copy();
	}

	/** Gets the size of the array in the given dimension
	 * @param dim Dimension from which the size will be taken.
	 * @return the size */
	public int getSize(int dim)
	{
		return this.size.get(dim);
	}

	/** Returns whether the receiver is a view or not. Being a view means that the data in the array is backed by the data in other array, so
	 * changes in one side are automatically reflected in the other
	 * @return true if the receiver is a view, false otherwise */
	public boolean isView()
	{
		return this.x.isView();
	}

	/** Normalizes this matrix, i.e. makes the sum of all elements equal to 1.0 If the matrix contains negative elements then all the values are
	 * shifted to ensure non-negativity. */
	public void normalize()
	{
		this.x.normalize();
	}

	/** Changes the shape of this array. The number of elements is the same, the elements are not copied.
	 * @param newSize new size of the array. The array must have the same number of cells
	 * @return new array matrix */
	public DoubleMatrixND reshape(int[] newSize)
	{
		if (IntMatrixND.prod(newSize) != this.numElem) throw new RuntimeException("Wrong size of n-array");
		this.size = IntFactory1D.dense.make(newSize);
		this.numDim = newSize.length;
		return this;
	}

	/** Sets the matrix cell at coordinate given by the subindexes to the specified value.
	 * @param index the index of the cell.
	 * @param value the value to be filled into the specified cell. */
	public void set(int index, double value)
	{
		x.set(index, value);
	}

	/** Sets the matrix cell at coordinate given by the subindexes to the specified value.
	 * @param subindexes the cell coordinates
	 * @param value the value to be filled into the specified cell. */
	public void set(int[] subindexes, double value)
	{
		set(IntFactory1D.dense.make(subindexes), value);
	}

	void set(IntMatrix1D subindexes, double value)
	{
		x.set(DoubleMatrixND.sub2ind(subindexes, this.size.toArray()), value);
	}

	/** Sets the matrix cell at coordinate index to the specified value. Provided with invalid parameters this method may access illegal indexes
	 * without throwing any exception. You should only use this method when you are absolutely sure that the coordinate is within bounds.
	 * @param index the index of the cell.
	 * @param value the value to be filled into the specified cell */
	public void setQuick(int index, double value)
	{
		x.setQuick(index, value);
	}

	/** Sets the matrix cell at coordinate given by the subindexes to the specified value. Provided with invalid parameters this method may access
	 * illegal indexes without throwing any exception. You should only use this method when you are absolutely sure that the coordinate is within
	 * bounds.
	 * @param subindexes the coordinates of the cell.
	 * @param value the value to be filled into the specified cell */
	public void setQuick(int[] subindexes, double value)
	{
		set(IntFactory1D.dense.make(subindexes), value);
	}

	void setQuick(IntMatrix1D subindexes, double value)
	{
		x.setQuick(DoubleMatrixND.sub2ind(subindexes, this.size.toArray()), value);
	}

	/** Constructs and returns a 1-dimensional array containing the array values. The values are copied.
	 * @return an array filled with the values of the cells */
	public double[] to1DArray()
	{
		return this.x.toArray();
	}

	/** Constructs and returns an n-dimensional double array containing the array values. The values are copied. If the array has 2 dimensions
	 * returns a double [][] object If the array has 3 dimensions returns a double [][][] object If the array has 4 dimensions returns a double
	 * [][][][] object If the array has 5 dimensions returns a double [][][][][] object If the array has 6 dimensions returns a double [][][][][][]
	 * object For arrays of more dimensions, an exception is raised
	 * @return an array filled with the values of the cells */
	public Object toArray()
	{
		switch (this.numDim)
		{
			case 1:
				return this.x.toArray();
			case 2:
			{
				double[][] res = new double[size.get(0)][size.get(1)];
				int cont = 0;
				for (int i1 = 0; i1 < size.get(1); i1++)
					for (int i0 = 0; i0 < size.get(0); i0++)
						res[i0][i1] = this.x.get(cont++);
				return res;
			}
			case 3:
			{
				double[][][] res = new double[size.get(0)][size.get(1)][size.get(2)];
				int cont = 0;
				for (int i2 = 0; i2 < size.get(2); i2++)
					for (int i1 = 0; i1 < size.get(1); i1++)
						for (int i0 = 0; i0 < size.get(0); i0++)
							res[i0][i1][i2] = this.x.get(cont++);
				return res;
			}
			case 4:
			{
				double[][][][] res = new double[size.get(0)][size.get(1)][size.get(2)][size.get(3)];
				int cont = 0;
				for (int i3 = 0; i3 < size.get(3); i3++)
					for (int i2 = 0; i2 < size.get(2); i2++)
						for (int i1 = 0; i1 < size.get(1); i1++)
							for (int i0 = 0; i0 < size.get(0); i0++)
								res[i0][i1][i2][i3] = this.x.get(cont++);
				return res;
			}
			case 5:
			{
				double[][][][][] res = new double[size.get(0)][size.get(1)][size.get(2)][size.get(3)][size.get(4)];
				int cont = 0;
				for (int i4 = 0; i4 < size.get(4); i4++)
					for (int i3 = 0; i3 < size.get(3); i3++)
						for (int i2 = 0; i2 < size.get(2); i2++)
							for (int i1 = 0; i1 < size.get(1); i1++)
								for (int i0 = 0; i0 < size.get(0); i0++)
									res[i0][i1][i2][i3][i4] = this.x.get(cont++);
				return res;
			}
			case 6:
			{
				double[][][][][][] res = new double[size.get(0)][size.get(1)][size.get(2)][size.get(3)][size.get(4)][size.get(5)];
				int cont = 0;
				for (int i5 = 0; i5 < size.get(5); i5++)
					for (int i4 = 0; i4 < size.get(4); i4++)
						for (int i3 = 0; i3 < size.get(3); i3++)
							for (int i2 = 0; i2 < size.get(2); i2++)
								for (int i1 = 0; i1 < size.get(1); i1++)
									for (int i0 = 0; i0 < size.get(0); i0++)
										res[i0][i1][i2][i3][i4][i5] = this.x.get(cont++);
				return res;
			}

		}
		throw new RuntimeException("Too many dimensions to use toArray ()");
	}

	@Override
	public String toString()
	{
		if (this.numDim == 1)
		{
			String s = "";
			for (int cont = 0; cont < x.size() - 1; cont++)
				s += this.x.get(cont) + " ; ";
			if (x.size() != 0) s += this.x.get((int) x.size() - 1);
			return s;
		}

		if (this.numDim == 2)
		{
			String s = "\n";
			for (int r = 0; r < size.get(0); r++)
			{
				for (int c = 0; c < size.get(1) - 1; c++)
					s += this.get(new int[]{r, c}) + " ; ";
				s = s + this.get(new int[]{r, size.get(1) - 1}) + ((r == size.get(0) - 1) ? "" : " ;;") + "\n";
			}
			return s;
		}

		//		String s = "Size: " + Arrays.toString(size.toArray()) + ", isView: " + this.x.isView() + ", values: ";
		//		s = s + "\n";
		String s = "";
		IntMatrix1D size2DmatrixToPrint = this.size.viewPart(0, 2).copy();
		IntMatrix1D[] ranges = new IntMatrix1D[this.numDim];
		for (int d = 2; d < this.numDim; d++)
			ranges[d] = IntFactory1D.dense.make(1);
		ranges[0] = null;
		ranges[1] = null;
		for (int cont = 0; cont < IntMatrixND.prod(this.size.viewPart(2, numDim - 2)); cont++)
		{
			IntMatrix1D subindexes = DoubleMatrixND.ind2sub(cont, this.size.viewPart(2, numDim - 2).toArray());
			s = s + "(all,all";
			for (int i : subindexes.toArray())
				s = s + "," + i;
			s += ")";
			for (int d = 2; d < this.numDim; d++)
				ranges[d].set(0, subindexes.get(d - 2));
			s = s + this.viewSelection(ranges).reshape(size2DmatrixToPrint.toArray()).toString() + "\n";
		}
		return s;
	}

	/** If the array has only one cell, returns its value. Otherwise raises an exception
	 * @return The value */
	public double toValue()
	{
		if (this.numElem != 1) throw new RuntimeException("This is not a scalar");
		return this.x.get(0);
	}

	/** If the array is 2D, returns a DoubleMatrix2D object copying the data in the array
	 * @return The DoubleMatrix2D object */
	public DoubleMatrix2D view2D()
	{
		if (this.numDim != 2) throw new RuntimeException("The N-DIM array must have two dimensions");
		return DoubleFactory2D.dense.make(this.x.toArray(), this.size.get(0));
	}

	/** If the array is 2D, returns a DoubleMatrix2D object copying the data in the array
	 * @return The DoubleMatrix2D object */
	public DoubleMatrix1D view1D()
	{
		if (this.numDim != 2) throw new RuntimeException("The N-DIM array must be a row or column vector");
		if ((this.size.get(0) != 1) && (this.size.get(1) != 1)) throw new RuntimeException("The N-DIM array must be a row or column vector");
		return x.copy();
	}

	/** If the array is 3D, returns a DoubleMatrix3D object copying the data in the array
	 * @param type Type of the matrix.
	 * @return The DoubleMatrix3D object */
	public DoubleMatrix3D view3D(String type)
	{
		if (this.numDim != 3) throw new RuntimeException("The N-DIM array must have three dimensions");
		if (!type.equalsIgnoreCase("sparse") && !type.equalsIgnoreCase("dense"))
			throw new RuntimeException("Wrong type of matrix. Please use \"dense\" or \"sparse\"");
		DoubleMatrix3D res = (type.equalsIgnoreCase("sparse")) ?
				DoubleFactory3D.sparse.make(size.get(0), size.get(1), size.get(2)) :
				DoubleFactory3D.dense.make(size.get(0), size.get(1), size.get(2));
		for (int c0 = 0; c0 < size.get(0); c0++)
			for (int c1 = 0; c1 < size.get(1); c1++)
				for (int c2 = 0; c2 < size.get(2); c2++)
					res.set(c0, c1, c2, this.get(new int[]{c0, c1, c2}));
		return res;
	}

	/** Returns Constructs and returns a new sub-range view. The view contains the cells from coordinates initialSubindex to
	 * initialSubindex+width-1, all inclusive. It is equivalent to viewSelection where in the i-th dimension, the ranges are (in JOM notation)
	 * (initialSubindex(i)) : (initialSubindex(i)+width(i)-1). In JOM notation Note that the view is really just a range restriction: The returned
	 * matrix is backed by this matrix, so changes in the returned matrix are reflected in this matrix, and vice-versa
	 * @param initialSubindex The coordinates of the cell that will is the upper-left corner of the subrange. This cell will have index 0 in the new
	 *                           view.
	 * @param width For each dimension, its size in the new view.
	 * @return the new view */
	public DoubleMatrixND viewPart(IntMatrix1D initialSubindex, IntMatrix1D width)
	{
		if ((int) initialSubindex.size() != this.numDim) throw new RuntimeException("Wrong number of dimensions in n-dim array");
		if ((int) width.size() != this.numDim) throw new RuntimeException("Wrong number of dimensions in n-dim array");
		if (initialSubindex.getMinLocation()[0] < 0) throw new RuntimeException("Wrong size of n-dim array");

		IntMatrix1D newSize = width.copy();
		IntMatrix1D binVectorIndexes = IntFactory1D.sparse.make(1, 1);

		for (int dim = this.numDim - 1; dim >= 0; dim--)
		{
			int originalSizeThisDim = this.size.get(dim);
			int initialSubindex_thisDim = initialSubindex.get(dim);
			int width_thisDim = width.get(dim);
			if (initialSubindex_thisDim + width_thisDim > originalSizeThisDim) throw new RuntimeException("Wrong size of n-dim array");

			IntMatrix1D binVectorCoordThisDim = IntFactory1D.sparse.make(originalSizeThisDim);
			binVectorCoordThisDim.viewPart(initialSubindex_thisDim, width_thisDim).assign(1);
			binVectorIndexes = IntMatrixND.kron(binVectorIndexes, binVectorCoordThisDim);
		}

		if (this.numElem != (int) binVectorIndexes.size()) throw new RuntimeException("Unexpected error");
		IntArrayList indexesList = new IntArrayList();
		IntArrayList auxList = new IntArrayList();
		binVectorIndexes.getNonZeros(indexesList, auxList);
		indexesList.trimToSize();
		if (IntMatrixND.prod(newSize) != indexesList.size()) throw new RuntimeException("Unexpected error");
		/* in the returned object, x will be a view */
		return new DoubleMatrixND(newSize.toArray(), this.x.viewSelection(indexesList.elements()));

	}

	/** Returns Constructs and returns a new sub-range view. The operation is equal to how subarrays are extracted in JOM. See JOM documentation for
	 *  more info. Note that the view is really just a range restriction: The returned matrix is backed by this matrix, so changes in the returned
	 *  matrix are reflected in this matrix, and vice-versa
	 * @param subindexesPerDim One IntMatrix1D object for each dimension, containing the indexes in the subrange.
	 * @return the new view */
	public DoubleMatrixND viewSelection(IntMatrix1D[] subindexesPerDim)
	{
		if (subindexesPerDim.length != this.numDim) throw new RuntimeException("Wrong number of dimensions in n-dim array");
		if (subindexesPerDim[0] != null)
		{
			if (subindexesPerDim[0].size() > 0)
				if (subindexesPerDim[0].getMinLocation()[0] < 0) throw new RuntimeException("Wrong size of n-dim array");
			if (subindexesPerDim[0].size() > 0)
				if (subindexesPerDim[0].getMaxLocation()[0] >= size.get(0)) throw new RuntimeException("Wrong size of n-dim array");
		}

		IntMatrix1D blockSizesOriginalArray = DoubleMatrixND.computeBlockSizes(size.toArray());
		IntMatrix1D newSize = this.size.copy();
		if (subindexesPerDim[0] != null) newSize.set(0, (int) subindexesPerDim[0].size());
		IntMatrix1D indexesPreviousDimensions = (subindexesPerDim[0] == null) ?
				new IntMatrixND(new int[]{1, size.get(0)}).ascending(0, 1).elements() :
				subindexesPerDim[0];
		for (int dim = 1; dim < subindexesPerDim.length; dim++)
		{
			int originalSizeThisDim = this.size.get(dim);
			int newSizeThisDim = (subindexesPerDim[dim] == null) ? originalSizeThisDim : (int) subindexesPerDim[dim].size();
			if (subindexesPerDim[dim] != null)
			{
				if (subindexesPerDim[dim].size() > 0)
					if (subindexesPerDim[dim].getMinLocation()[0] < 0) throw new RuntimeException("Wrong size of n-dim array");
				if (subindexesPerDim[dim].size() > 0)
					if (subindexesPerDim[dim].getMaxLocation()[0] >= originalSizeThisDim) throw new RuntimeException("Wrong size of n-dim array");
			}
			newSize.set(dim, newSizeThisDim);

			IntMatrix1D indexesInThisDimension = IntFactory1D.dense.make(((int) indexesPreviousDimensions.size()) * newSizeThisDim);
			for (int contNewCoord = 0; contNewCoord < newSizeThisDim; contNewCoord++)
			{
				int indexThisCoord = (subindexesPerDim[dim] == null) ? contNewCoord : subindexesPerDim[dim].get(contNewCoord);
				IntMatrix1D indexesPrevDimensionsRecomputed = indexesPreviousDimensions.copy().assign(IntFunctions.plus(indexThisCoord *
						blockSizesOriginalArray.get(dim)));
				indexesInThisDimension.viewPart(contNewCoord * (int) indexesPreviousDimensions.size(), (int) indexesPreviousDimensions.size())
						.assign(indexesPrevDimensionsRecomputed);
			}
			indexesPreviousDimensions = indexesInThisDimension;
		}
		if (IntMatrixND.prod(newSize) != indexesPreviousDimensions.size())
			throw new RuntimeException("Unexpected error: IntMatrixND.prod(newSize)" + IntMatrixND.prod(newSize) + " , indexesPreviousDimensions"
					+ ".size():" + indexesPreviousDimensions.size());
		return new DoubleMatrixND(newSize.toArray(), this.x.viewSelection(indexesPreviousDimensions.toArray()));
	}

	/** Returns Constructs and returns a new sub-range view. The elements in the new view are the ones given in the array of indexes. Note that the
	 * view is really just a range restriction: The returned matrix is backed by this matrix, so changes in the returned matrix are reflected in
	 * this matrix, and vice-versa
	 * @param newSize The new size of the array. The associated number of cells in the new array must be equal to the number of indexes in indexes
	 * @param indexes An array with the inedexes to include in the view (in the same order as they appear here)
	 * @return the new view */
	public DoubleMatrixND viewSelectionByIndexes(int[] newSize, IntMatrix1D indexes)
	{
		if (IntMatrixND.prod(newSize) != (int) indexes.size()) throw new RuntimeException("Size mismatch");
		return new DoubleMatrixND(newSize, this.x.viewSelection(indexes.toArray()));
	}

	/* Constructs and returns a new stride view which is a sub matrix consisting of every i-th cell. More specifically, the view has this.slices()
	/sliceStride slices and this.rows()/rowStride rows and this.columns()/columnStride columns holding cells this.get(k*sliceStride,i*rowStride,
	j*columnStride) for all k = 0..slices()/sliceStride - 1, i = 0..rows()/rowStride - 1, j = 0..columns()/columnStride - 1 . The returned view is
	backed by this matrix, so changes in the returned view are reflected in this matrix, and vice-versa. */
	DoubleMatrixND viewStrides(IntMatrix1D strides)
	{
		if ((int) strides.size() != this.numDim) throw new RuntimeException("Wrong number of dimensions in n-dim array");
		if (strides.getMinLocation()[0] < 0) throw new RuntimeException("Wrong size of n-dim array");

		IntMatrix1D newSize = this.size.like();
		IntMatrix1D binVectorIndexes = IntFactory1D.sparse.make(1, 1);

		for (int dim = this.numDim - 1; dim >= 0; dim--)
		{
			int originalSizeThisDim = this.size.get(dim);
			int strideThisDim = strides.get(dim);
			newSize.set(dim, originalSizeThisDim / strideThisDim);

			IntMatrix1D binVectorCoordThisDim = IntFactory1D.sparse.make(originalSizeThisDim);
			binVectorCoordThisDim.viewStrides(strideThisDim).assign(1);
			binVectorIndexes = IntMatrixND.kron(binVectorIndexes, binVectorCoordThisDim);
		}

		if (this.numElem != (int) binVectorIndexes.size()) throw new RuntimeException("Unexpected error");
		IntArrayList indexesList = new IntArrayList();
		IntArrayList auxList = new IntArrayList();
		binVectorIndexes.getNonZeros(indexesList, auxList);
		indexesList.trimToSize();
		if (IntMatrixND.prod(newSize) != indexesList.size()) throw new RuntimeException("Unexpected error");
		/* in the returned object, x will be a view */
		return new DoubleMatrixND(newSize.toArray(), this.x.viewSelection(indexesList.elements()));
	}

	/** Returns the sum of all cells in the array
	 * @return the sum */
	public double zSum()
	{
		return this.x.zSum();
	}
}
