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

import cern.colt.function.tint.IntFunction;
import cern.colt.function.tint.IntIntFunction;
import cern.colt.function.tint.IntProcedure;
import cern.colt.list.tint.IntArrayList;
import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntFactory2D;
import cern.colt.matrix.tint.IntMatrix1D;
import cern.colt.matrix.tint.IntMatrix2D;
import cern.jet.math.tint.IntPlusMultSecond;

class IntMatrixND
{
	private int numDim;

	private int numElem;

	private IntMatrix1D size;

	private IntMatrix1D x;

	public IntMatrixND(int value)
	{
		this(IntFactory1D.dense.make(new int[]{1, 1}), IntFactory1D.dense.make(1, value));
	}

	public IntMatrixND(int[] size)
	{
		this(size, "sparse");
	}

	public IntMatrixND(int[] size, double[] values)
	{
		this(size, values, "sparse");
	}

	public IntMatrixND(int[] size, double[] values, String type)
	{
		this(IntFactory1D.dense.make(size), (type.equals("dense")) ? DoubleFactory1D.dense.make(values) : DoubleFactory1D.sparse.make(values));
	}

	public IntMatrixND(int[] size, DoubleMatrix1D values)
	{
		this(IntFactory1D.dense.make(size), values);
	}

	public IntMatrixND(int[] size, int value)
	{
		this(size, value, "sparse");
	}

	public IntMatrixND(int[] size, int value, String type)
	{
		this(IntFactory1D.dense.make(size), value, type);
	}

	public IntMatrixND(int[] size, int[] values)
	{
		this(size, values, "sparse");
	}

	public IntMatrixND(int[] size, int[] values, String type)
	{
		this(IntFactory1D.dense.make(size), values, type);
	}

	public IntMatrixND(int[] size, IntMatrix1D values)
	{
		this(IntFactory1D.dense.make(size), values);
	}

	public IntMatrixND(int[] size, String type)
	{
		this(IntFactory1D.dense.make(size), type);
	}

	public IntMatrixND(int[][] values)
	{
		this(values, "sparse");
	}

	public IntMatrixND(int[][] values, String type)
	{
		this((type.equals("dense")) ? IntFactory2D.dense.make(values) : IntFactory2D.sparse.make(values));
	}

	public IntMatrixND(IntMatrix1D size)
	{
		this(size, "sparse");
	}

	public IntMatrixND(IntMatrix1D size, DoubleMatrix1D values)
	{
		this(size, new DoubleMatrixND(size.toArray(), values).cast2Int().elements());
	}

	public IntMatrixND(IntMatrix1D size, int value)
	{
		this(size, value, "sparse");
	}

	public IntMatrixND(IntMatrix1D size, int value, String type)
	{
		this(size, (type.equals("dense")) ?
				IntFactory1D.dense.make(IntMatrixND.prod(size), value) :
				IntFactory1D.sparse.make(IntMatrixND.prod(size), value));
	}

	public IntMatrixND(IntMatrix1D size, int[] values)
	{
		this(size, values, "sparse");
	}

	public IntMatrixND(IntMatrix1D size, int[] values, String type)
	{
		this(size, (type.equals("dense")) ? IntFactory1D.dense.make(values) : IntFactory1D.sparse.make(values));
	}

	/****************** CONSTRUCTORS **********************************/
	public IntMatrixND(IntMatrix1D size, IntMatrix1D values)
	{
		this.size = size.copy();
		this.size.trimToSize();
		this.numDim = (int) size.size();
		if (size.getMinLocation()[0] <= 0) throw new RuntimeException("Wrong size of n-array");
		this.numElem = IntMatrixND.prod(size);
		if (this.numElem != (int) values.size()) throw new RuntimeException("Wrong size of n-array");
		this.x = values;
		this.x.trimToSize();
	}

	public IntMatrixND(IntMatrix1D size, String type)
	{
		this(size, 0, type);
	}

	public IntMatrixND(IntMatrix2D values)
	{
		this(IntFactory1D.dense.make(new int[]{values.rows(), values.columns()}), values.copy().vectorize());
	}

	/* Converts the index of the cell into its coordinates (subindex), in an array of the given size */
	public static IntMatrix1D ind2sub(int index, IntMatrix1D size)
	{
		IntMatrix2D res = IntMatrixND.ind2sub(IntFactory1D.dense.make(new int[]{index}), size);
		return res.viewRow(0).copy();
	}

	/* Converts the indexes of the cells into their coordinates (subindex), in an array of the given size. i-th row of the returned array contains
	the subindexes of i-th cell */
	public static IntMatrix2D ind2sub(IntMatrix1D indexes, IntMatrix1D size)
	{
		int numCells = IntMatrixND.prod(size);
		int numDim = (int) size.size();
		int numIndexes = (int) indexes.size();
		if (indexes.size() > 0) if (indexes.getMaxLocation()[0] > numCells - 1){ throw new RuntimeException("Index out of bounds in n-dim array"); }
		if (indexes.size() > 0) if (indexes.getMinLocation()[0] < 0){ throw new RuntimeException("Index out of bounds in n-dim array"); }

		IntMatrix1D blockSizes = IntMatrixND.computeBlockSizes(size);
		IntMatrix2D res = IntFactory2D.dense.make((int) indexes.size(), numDim);
		IntMatrix1D indexesToReduce = indexes.copy();
		for (int contDim = numDim - 1; contDim >= 0; contDim--)
		{
			IntMatrix1D thisBlockSizeVector = IntFactory1D.dense.make(numIndexes, blockSizes.get(contDim));
			IntMatrix1D coord_n = indexesToReduce.copy().assign(thisBlockSizeVector, cern.jet.math.tint.IntFunctions.div);
			res.viewColumn(contDim).assign(coord_n);

			indexesToReduce.assign(thisBlockSizeVector, cern.jet.math.tint.IntFunctions.mod);
		}
		return res;
	}

	public static IntMatrix1D kron(IntMatrix1D a, IntMatrix1D b)
	{
		IntMatrix1D res = IntFactory1D.sparse.make((int) (a.size() * b.size()));
		IntArrayList indexList = new IntArrayList();
		IntArrayList valueList = new IntArrayList();
		a.getNonZeros(indexList, valueList);
		for (int cont = 0; cont < indexList.size(); cont++)
		{
			/* for each non-zero value in a, set a block in res = (value * b) */
			int index = indexList.get(cont);
			int value = valueList.get(cont);
			res.viewPart(index * (int) b.size(), (int) b.size()).assign(b.copy(), IntPlusMultSecond.plusMult(value));
		}
		return res;
	}

	/************************ STATIC METHODS ***************************/
	public static IntMatrix2D kron(IntMatrix2D a, IntMatrix2D b)
	{
		IntMatrix2D res = IntFactory2D.sparse.make(a.rows() * b.rows(), a.columns() * b.columns());
		IntArrayList rowList = new IntArrayList();
		IntArrayList columnList = new IntArrayList();
		IntArrayList valueList = new IntArrayList();
		a.getNonZeros(rowList, columnList, valueList);
		for (int cont = 0; cont < rowList.size(); cont++)
		{
			/* for each non-zero value in a, set a block in res = (value * b) */
			int row = rowList.get(cont);
			int column = columnList.get(cont);
			int value = valueList.get(cont);
			res.viewPart(row * b.rows(), column * b.columns(), b.rows(), b.columns()).assign(b.copy(), IntPlusMultSecond.plusMult(value));
		}
		return res;
	}

	public static int prod(int[] v)
	{
		int accum = 1;
		for (int cont = 0; cont < v.length; cont++)
			accum *= v[cont];
		return accum;
	}

	public static int prod(int[] v, int firstIndex, int numIndexes)
	{
		int accum = 1;
		for (int cont = 0; cont < numIndexes; cont++)
			accum *= v[firstIndex + cont];
		return accum;
	}

	public static int prod(IntMatrix1D v)
	{
		int accum = 1;
		for (int cont = 0; cont < v.size(); cont++)
			accum *= v.get(cont);
		return accum;
	}

	/* Converts the subindexes of a cell given by coord, into its index, in an array of the given size */
	public static int sub2ind(IntMatrix1D coord, IntMatrix1D size)
	{
		IntMatrix1D blockSizes = IntMatrixND.computeBlockSizes(size);
		return coord.zDotProduct(blockSizes);
	}

	/* Converts the subindexes of the cells (one cell per row of coord) into their indexes, in an array of the given size */
	public static IntMatrix1D sub2ind(IntMatrix2D coord, IntMatrix1D size)
	{
		IntMatrix1D blockSizes = IntMatrixND.computeBlockSizes(size);
		IntMatrix1D indexes = coord.zMult(blockSizes, null);
		return indexes;
	}

	/******************** PRIVATE METHODS *********************/

	private static IntMatrix1D computeBlockSizes(IntMatrix1D size)
	{
		IntMatrix1D blockSizes = IntFactory1D.dense.make((int) size.size());
		blockSizes.set(0, 1);
		for (int dim = 1; dim < (int) size.size(); dim++)
			blockSizes.set(dim, blockSizes.get(dim - 1) * size.get(dim - 1));
		return blockSizes;
	}

	/****************** AGGREGATE METHODS ********************************/
	/* Applies a function to each cell and aggregates the results. */
	public int aggregate(IntIntFunction aggr, IntFunction f)
	{
		return x.aggregate(aggr, f);
	}

	/* Applies a function to all cells with a given indexes and aggregates the results. */
	public int aggregate(IntIntFunction aggr, IntFunction f, IntArrayList indexList)
	{
		return x.aggregate(aggr, f, indexList);
	}

	/* Applies a function to each corresponding cell of two matrices and aggregates the results. */
	public int aggregate(IntMatrixND other, IntIntFunction aggr, IntIntFunction f)
	{
		return x.aggregate(other.elements(), aggr, f);
	}

	public IntMatrixND ascending(int first, int step)
	{
		for (int cont = 0; cont < this.numElem; cont++)
			this.x.set(cont, first + cont * step);
		return this;
	}

	/****************** ASSIGN METHODS ********************************/
	/* Sets all cells to the state specified by value. */
	public IntMatrixND assign(int value)
	{
		x.assign(value);
		return this;
	}

	/* Sets all cells to the state specified by values. values is required to have as many cells as this. The values are copied. So subsequent
	changes in values are not reflected in the matrix, and vice-versa. */
	public IntMatrixND assign(int[] values)
	{
		this.x.assign(values);
		return this;
	}

	/* Assigns the result of a function to each cell */
	public IntMatrixND assign(IntFunction f)
	{
		this.x.assign(f);
		return this;
	}

	/* Replaces all cell values of the receiver with the values of another matrix. Both matrices must have the same size. If both matrices share the
	 same cells (as is the case if they are views derived from the same matrix) and intersect in an ambiguous way, then replaces as if using an
	 intermediate auxiliary deep copy of other. */
	public IntMatrixND assign(IntMatrixND other)
	{
		if (!this.size.equals(other.getSize())) throw new RuntimeException("Wrong size of n-array");
		this.x.assign(other.elements());
		return this;
	}

	/* Assigns the result of a function to each cell; x[cell] = function(x[cell],y[cell]). */
	public IntMatrixND assign(IntMatrixND y, IntIntFunction function)
	{
		this.x.assign(y.elements(), function);
		return this;
	}

	/* Assigns the result of a function to all cells with a given indexes */
	/* function - a function object taking as first argument the current cell's value of this, and as second argument the current cell's value of
	y, */
	public IntMatrixND assign(IntMatrixND y, IntIntFunction function, IntArrayList indexList)
	{
		this.x.assign(y.elements(), function, indexList);
		return this;
	}

	/* Assigns a value to all cells that satisfy a condition. */
	public IntMatrixND assign(IntProcedure cond, int value)
	{
		this.x.assign(cond, value);
		return this;
	}

	/* Assigns the result of a function to all cells that satisfy a condition. */
	public IntMatrixND assign(IntProcedure cond, IntFunction f)
	{
		this.x.assign(cond, f);
		return this;
	}

	/*************** MISCELANEOUS METHODS ****************************/
	/* Returns the number of cells having non-zero values; ignores tolerance. */
	public int cardinality()
	{
		return x.cardinality();
	}

	/* Creates a new copy of the array, with all the elements transformed to integer */
	public DoubleMatrixND cast2Double()
	{
		DoubleMatrix1D res = DoubleFactory1D.dense.make(this.numElem);
		IntArrayList indexList = new IntArrayList();
		IntArrayList valueList = new IntArrayList();
		getNonZeros(indexList, valueList);
		indexList.trimToSize();
		int[] indexesNonZero = indexList.elements();
		for (int cont = 0; cont < indexesNonZero.length; cont++)
			res.set(indexesNonZero[cont], this.x.get(indexesNonZero[cont]));
		return new DoubleMatrixND(this.size.toArray(), res);
	}

	/* Constructs and returns a deep copy of the receiver */
	public IntMatrixND copy()
	{
		return new IntMatrixND(this.size, this.x.copy());
	}

	public IntMatrixND descending(int first, int step)
	{
		for (int cont = 0; cont < this.numElem; cont++)
			this.x.set(cont, first - cont * step);
		return this;
	}

	/* Returns the elements of this matrix. */
	public IntMatrix1D elements()
	{
		return x;
	}

	/* Returns whether all cells are equal to the given value */
	public boolean equals(double value)
	{
		return x.equals(value);
	}

	/* Compares this object against the specified object. */
	public boolean equals(IntMatrixND other)
	{
		return (this.x.equals(other.elements()) && (this.size.equals(other.getSize())));
	}

	/* Returns the matrix cell value at the given index. */
	public double get(int index)
	{
		return x.get(index);
	}

	/* Returns the matrix cell value at the given subindexes */
	public double get(int[] subindexes)
	{
		return x.get(DoubleMatrixND.sub2ind(IntFactory1D.dense.make(subindexes), this.size));
	}

	/* Returns the matrix cell value at the given subindexes */
	public double get(IntMatrix1D subindexes)
	{
		return x.get(IntMatrixND.sub2ind(subindexes, this.size));
	}

	/* Return maximum value of this matrix together with its location. If the array is empty, returns and empty array */
	public int[] getMaxLocation()
	{
		return x.size() > 0 ? this.x.getMaxLocation() : new int[]{};
	}

	/* Return minimum value of this matrix together with its location. If the array is empty, returns and empty array */
	public int[] getMinLocation()
	{
		return x.size() > 0 ? this.x.getMinLocation() : new int[]{};
	}

	/* Fills the indexes and values of cells having negative values into the specified lists. */
	public void getNegativeValues(IntArrayList indexList, IntArrayList valueList)
	{
		this.x.getNegativeValues(indexList, valueList);
	}

	/* Fills the indexes and values of cells having non zero values into the specified lists. */
	public void getNonZeros(IntArrayList indexList, IntArrayList valueList)
	{
		this.x.getNonZeros(indexList, valueList);
	}

	public int getNumDim()
	{
		return this.numDim;
	}

	public int getNumElements()
	{
		return this.numElem;
	}

	/* Fills the indexes and values of cells having positive values into the specified lists. */
	public void getPositiveValues(IntArrayList indexList, IntArrayList valueList)
	{
		this.x.getPositiveValues(indexList, valueList);
	}

	/* Returns the matrix cell value at index given. */
	public double getQuick(int index)
	{
		return this.x.getQuick(index);
	}

	public IntMatrix1D getSize()
	{
		return this.size.copy();
	}

	public int getSize(int dim)
	{
		return this.size.get(dim);
	}

	/* Returns whether the receiver is a view or not. */
	public boolean isView()
	{
		return this.x.isView();
	}

	/* Construct and returns a new empty matrix of the same dynamic type as the receiver, having the same number of slices, rows and columns. */
	public IntMatrixND like()
	{
		return new IntMatrixND(this.size, this.x.like());
	}

	/* Construct and returns a new empty matrix of the same dynamic type as the receiver, having the specified size. */
	public IntMatrixND like(IntMatrix1D newSize)
	{
		return new IntMatrixND(newSize, this.x.like(IntMatrixND.prod(newSize)));
	}

	/* Changes the shape of this array. The number of elements is the same, the elements are not copied. */
	public IntMatrixND reshape(IntMatrix1D newSize)
	{
		if (IntMatrixND.prod(newSize) != this.numElem) throw new RuntimeException("Wrong size of n-array");
		this.size = newSize.copy();
		this.numDim = (int) size.size();
		return this;
	}

	/* Sets the matrix cell at coordinate given by the subindexes to the specified value. */
	public void set(int index, int value)
	{
		x.set(index, value);
	}

	/* Sets the matrix cell at coordinate given by the subindexes to the specified value. */
	public void set(IntMatrix1D subindexes, int value)
	{
		x.set(IntMatrixND.sub2ind(subindexes, this.size), value);
	}

	/* Sets the matrix cell at coordinate given by the subindexes to the specified value. */
	public void setQuick(int index, int value)
	{
		x.setQuick(index, value);
	}

	/* Sets the matrix cell at coordinate given by the subindexes to the specified value. */
	public void setQuick(IntMatrix1D subindexes, int value)
	{
		x.setQuick(IntMatrixND.sub2ind(subindexes, this.size), value);
	}

	public Object toArray()
	{
		switch (this.numDim)
		{
			case 1:
				return this.x.toArray();
			case 2:
			{
				int[][] res = new int[size.get(0)][size.get(1)];
				int cont = 0;
				for (int i1 = 0; i1 < size.get(1); i1++)
					for (int i0 = 0; i0 < size.get(0); i0++)
						res[i0][i1] = this.x.get(cont++);
				return res;
			}
			case 3:
			{
				int[][][] res = new int[size.get(0)][size.get(1)][size.get(2)];
				int cont = 0;
				for (int i2 = 0; i2 < size.get(2); i2++)
					for (int i1 = 0; i1 < size.get(1); i1++)
						for (int i0 = 0; i0 < size.get(0); i0++)
							res[i0][i1][i2] = this.x.get(cont++);
				return res;
			}
			case 4:
			{
				int[][][][] res = new int[size.get(0)][size.get(1)][size.get(2)][size.get(3)];
				int cont = 0;
				for (int i3 = 0; i3 < size.get(3); i3++)
					for (int i2 = 0; i2 < size.get(2); i2++)
						for (int i1 = 0; i1 < size.get(1); i1++)
							for (int i0 = 0; i0 < size.get(0); i0++)
								res[i0][i1][i2][i3] = this.x.get(cont++);
			}
			case 5:
			{
				int[][][][][] res = new int[size.get(0)][size.get(1)][size.get(2)][size.get(3)][size.get(4)];
				int cont = 0;
				for (int i4 = 0; i4 < size.get(4); i4++)
					for (int i3 = 0; i3 < size.get(3); i3++)
						for (int i2 = 0; i2 < size.get(2); i2++)
							for (int i1 = 0; i1 < size.get(1); i1++)
								for (int i0 = 0; i0 < size.get(0); i0++)
									res[i0][i1][i2][i3][i4] = this.x.get(cont++);
			}
			case 6:
			{
				int[][][][][][] res = new int[size.get(0)][size.get(1)][size.get(2)][size.get(3)][size.get(4)][size.get(5)];
				int cont = 0;
				for (int i5 = 0; i5 < size.get(5); i5++)
					for (int i4 = 0; i4 < size.get(4); i4++)
						for (int i3 = 0; i3 < size.get(3); i3++)
							for (int i2 = 0; i2 < size.get(2); i2++)
								for (int i1 = 0; i1 < size.get(1); i1++)
									for (int i0 = 0; i0 < size.get(0); i0++)
										res[i0][i1][i2][i3][i4][i5] = this.x.get(cont++);
			}

		}
		throw new RuntimeException("Too many dimensions to use toArray ()");
	}

	@Override
	public String toString()
	{
		return this.cast2Double().toString();
	}

	public int toValue()
	{
		if (this.numElem != 1) throw new RuntimeException("This is not a scalar");
		return this.x.get(0);
	}

	public IntMatrix1D view1D()
	{
		if (this.numDim != 1) throw new RuntimeException("Tne N-DIM array must have one dimension");
		return IntFactory1D.dense.make(this.x.toArray());
	}

	public IntMatrix2D view2D()
	{
		if (this.numDim != 2) throw new RuntimeException("Tne N-DIM array must have two dimensions");
		return IntFactory2D.dense.make(this.x.toArray(), this.size.get(0));
	}

	public IntMatrixND viewPart(IntMatrix1D initialSubindex, IntMatrix1D width)
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
		return new IntMatrixND(newSize, this.x.viewSelection(indexesList.elements()));
	}

	public IntMatrixND viewSelection(IntMatrix1D[] subindexesPerDim)
	{
		if (subindexesPerDim.length != this.numDim) throw new RuntimeException("Wrong number of dimensions in n-dim array");

		IntMatrix1D newSize = this.size.like();
		IntMatrix1D binVectorIndexes = IntFactory1D.sparse.make(1, 1);

		for (int dim = subindexesPerDim.length - 1; dim >= 0; dim--)
		{
			int originalSizeThisDim = this.size.get(dim);
			IntMatrix1D binVectorCoordThisDim;
			if (subindexesPerDim[dim] == null)
			{
				newSize.set(dim, originalSizeThisDim);
				binVectorCoordThisDim = IntFactory1D.sparse.make(originalSizeThisDim, 1);
			} else
			{
				if (subindexesPerDim[dim].size() > 0)
					if (subindexesPerDim[dim].getMinLocation()[0] < 0) throw new RuntimeException("Wrong size of n-dim array");
				if (subindexesPerDim[dim].size() > 0)
					if (subindexesPerDim[dim].getMaxLocation()[0] >= originalSizeThisDim) throw new RuntimeException("Wrong size of n-dim array");
				newSize.set(dim, (int) subindexesPerDim[dim].size());
				binVectorCoordThisDim = IntFactory1D.sparse.make(originalSizeThisDim);
				binVectorCoordThisDim.viewSelection(subindexesPerDim[dim].toArray()).assign(1);
			}

			binVectorIndexes = IntMatrixND.kron(binVectorIndexes, binVectorCoordThisDim);
		}

		if (this.numElem != (int) binVectorIndexes.size()) throw new RuntimeException("Unexpected error");
		IntArrayList indexesList = new IntArrayList();
		IntArrayList auxList = new IntArrayList();
		binVectorIndexes.getNonZeros(indexesList, auxList);
		indexesList.trimToSize();
		if (IntMatrixND.prod(newSize) != indexesList.size()) throw new RuntimeException("Unexpected error");
		/* in the returned object, x will be a view */
		return new IntMatrixND(newSize, this.x.viewSelection(indexesList.elements()));
	}

	public IntMatrixND viewSelectionByIndexes(IntMatrix1D newSize, IntMatrix1D indexes)
	{
		if (IntMatrixND.prod(newSize) != (int) indexes.size()) throw new RuntimeException("Size mismatch");
		return new IntMatrixND(newSize, this.x.viewSelection(indexes.toArray()));
	}

	/* Constructs and returns a new stride view which is a sub matrix consisting of every i-th cell. More specifically, the view has this.slices()
	/sliceStride slices and this.rows()/rowStride rows and this.columns()/columnStride columns holding cells this.get(k*sliceStride,i*rowStride,
	j*columnStride) for all k = 0..slices()/sliceStride - 1, i = 0..rows()/rowStride - 1, j = 0..columns()/columnStride - 1 . The returned view is
	backed by this matrix, so changes in the returned view are reflected in
	 * this matrix, and vice-versa. */
	public IntMatrixND viewStrides(IntMatrix1D strides)
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
		return new IntMatrixND(newSize, this.x.viewSelection(indexesList.elements()));
	}

	/* Returns the sum of all cells */
	public int zSum()
	{
		return this.x.zSum();
	}

}
