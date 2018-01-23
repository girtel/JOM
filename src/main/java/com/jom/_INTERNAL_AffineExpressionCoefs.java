/**
 * Copyright (c) 2015 Pablo Pavon Mari�o.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * <p>
 * Contributors:
 * Pablo Pavon Mari�o - initial API and implementation
 *
 * @author Pablo Pav�n Mari�o, Technical University of Cartagena (Spain)
 * @version 0.1
 * 12/02/2013
 */

/**


 *
 * @author Pablo Pav�n Mari�o, Technical University of Cartagena (Spain)
 * @version 0.1
 * 12/02/2013
 */
package com.jom;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.Map.Entry;

import cern.colt.matrix.tdouble.DoubleFactory2D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;

/**
 * @author Pablo
 *
 */
class _INTERNAL_AffineExpressionCoefs
{
	private final OptimizationProblem model;
	private       int[] size;
	private       int numCells;
	private       int numDim;
	private       LinkedHashMap<Integer, Cell> linearCoefs; // first key is the cell index, second key is the variable index

	_INTERNAL_AffineExpressionCoefs(OptimizationProblem model, int[] size, LinkedHashMap<Integer, Cell> linearCoefs)
	{
		this.model = model;
		this.resize(size);
		this.linearCoefs = (linearCoefs == null) ? new LinkedHashMap<Integer, Cell>() : linearCoefs;
	}

	_INTERNAL_AffineExpressionCoefs(OptimizationProblem model, int[] size)
	{
		this.model = model;
		this.resize(size);
		this.linearCoefs = new LinkedHashMap<Integer, Cell>();
	}

	_INTERNAL_AffineExpressionCoefs(OptimizationProblem model, double constant)
	{
		this.model = model;
		this.resize(new int[]{1, 1});
		this.linearCoefs = new LinkedHashMap<Integer, Cell>();
		if (constant == 0) return;
		this.linearCoefs.put(0, new Cell(constant));
	}

	_INTERNAL_AffineExpressionCoefs(OptimizationProblem model, int[] size, double[] constants)
	{
		this.model = model;
		this.resize(size);
		if (constants.length != this.numCells) throw new JOMException("Affine expression operation error: Unexpected error");
		this.linearCoefs = new LinkedHashMap<Integer, Cell>();
		for (int cont = 0; cont < constants.length; cont++)
			if (constants[cont] != 0)
				this.linearCoefs.put(cont, new Cell(constants[cont]));
	}

	_INTERNAL_AffineExpressionCoefs(OptimizationProblem model, IntMatrixND varIds)
	{
		this.model = model;
		this.resize(Arrays.copyOf(varIds.getSize().toArray(), varIds.getNumDim()));
		this.linearCoefs = new LinkedHashMap<Integer, Cell>();
		for (int cellCont = 0; cellCont < numCells; cellCont++)
		{
			final int dv = (int) varIds.get(cellCont);
			LinkedHashMap<Integer, Double> v = new LinkedHashMap<Integer, Double>();
			v.put(dv, 1.0);
			this.linearCoefs.put(cellCont, new Cell(v));
		}
	}

	static _INTERNAL_AffineExpressionCoefs operator_multiply(_INTERNAL_AffineExpressionCoefs a, _INTERNAL_AffineExpressionCoefs b)
	{
		if (!a.isConstant() && !b.isConstant()) throw new JOMException("Affine expression operation error: Unexpected error");
		if ((a.numDim != 2) || (b.getNumDim() < 2)) throw new JOMException("Affine expression operation error: Unexpected error");
		if (a.size[1] != b.size[0]) throw new JOMException("Affine expression operation error: Unexpected error");
		final int aRows = a.size[0];
		final int aCols = a.size[1];
		final int bRows = b.size[0];
		final int bCols = b.size[1];
		final int resRows = aRows;
		final int resCols = bCols;
		final int[] resSize = Arrays.copyOf(b.size, b.size.length);
		resSize[0] = aRows;
		final int resNumCells = IntMatrixND.prod(resSize);
		final LinkedHashMap<Integer, Cell> aLC = a.linearCoefs;
		final LinkedHashMap<Integer, Cell> bLC = b.linearCoefs;
		//		final int [] resSize = new int [] { aRows , bCols };

		LinkedHashMap<Integer, Cell> res_linearCoefs = new LinkedHashMap<Integer, Cell>();

		/* Arrange in an array, for each column of a, the non-zero cells */
		ArrayList<LinkedList<Entry<Integer, Cell>>> a_perColumn = new ArrayList<LinkedList<Entry<Integer, Cell>>>(aCols);
		for (int cont = 0; cont < aCols; cont++) a_perColumn.add(new LinkedList<Entry<Integer, Cell>>());
		for (Entry<Integer, Cell> aCellEntry : aLC.entrySet())
		{
			final int aCol = aCellEntry.getKey() / aRows;
			a_perColumn.get(aCol).add(aCellEntry);
		}

		for (Entry<Integer, Cell> bCellEntry : bLC.entrySet())
		{
			final int bCellId = bCellEntry.getKey();
			final Cell bCell = bCellEntry.getValue();
			final IntMatrix1D subIndexesInB = DoubleMatrixND.ind2sub(bCellId, b.size);
			final int bCol = subIndexesInB.get(1);//bCellId / bRows;
			final int bRow = subIndexesInB.get(0);//bCellId % bRows;
			IntMatrix1D subIndexesInRes = subIndexesInB.copy();
			for (Entry<Integer, Cell> aCellEntry : a_perColumn.get(bRow))
			{
				final int aCellId = aCellEntry.getKey();
				final Cell aCell = aCellEntry.getValue();
				final int aRow = aCellId % aRows;
				final int aCol = aCellId / aRows;
				final int resRow = aRow;
				final int resCol = bCol;
				subIndexesInRes.set(0, resRow);
				final int resCellId = DoubleMatrixND.sub2ind(subIndexesInRes, resSize); //resRow + resRows * resCol;
				final Cell summand = aCell.deepCopy().multiply(bCell);
				if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
				Cell accumCell = res_linearCoefs.get(resCellId);
				if (accumCell == null)
					res_linearCoefs.put(resCellId, summand);
				else
					accumCell.sum(summand);
			}
		}

		return new _INTERNAL_AffineExpressionCoefs(a.model, resSize, res_linearCoefs);
	}

	private static int prod(int[] array)
	{
		int res = 1;
		for (int num : array) res *= num;
		return res;
	}

	static int[] toArrayInt(ArrayList<Integer> a)
	{
		int[] res = new int[a.size()];
		int counter = 0;
		for (int val : a) res[counter++] = val;
		return res;
	}

	static double[] toArrayDouble(ArrayList<Double> a)
	{
		double[] res = new double[a.size()];
		int counter = 0;
		for (double val : a) res[counter++] = val;
		return res;
	}

	static ArrayList<Integer> newAggregatedListNonOrderNonRepeatedNonZeros(final ArrayList<Integer> a, final ArrayList<Integer> b, final int maxNum)
	{
		ArrayList<Integer> res = new ArrayList<Integer>();
		boolean[] alreadyInserted = new boolean[maxNum];
		for (int val : a)
			if (!alreadyInserted[val])
			{
				res.add(val);
				alreadyInserted[val] = true;
			}
		for (int val : b)
			if (!alreadyInserted[val])
			{
				res.add(val);
				alreadyInserted[val] = true;
			}

		return res;
	}

	private static int[] computeBlockSizes(int[] size)
	{
		int[] blockSizes = new int[size.length];
		blockSizes[0] = 1;
		for (int dim = 1; dim < size.length; dim++)
			blockSizes[dim] = blockSizes[dim - 1] * size[dim - 1];
		return blockSizes;
	}

	private static int[] ind2sub(int index, int[] size, int[] blockSizes)
	{
		final int numDim = size.length;
		int[] res = new int[size.length];
		int indexToReduce = index;
		for (int contDim = numDim - 1; contDim >= 0; contDim--)
		{
			res[contDim] = indexToReduce / blockSizes[contDim];
			indexToReduce = indexToReduce % blockSizes[contDim];
		}
		return res;
	}

	DoubleMatrix2D getJacobian()
	{
		DoubleMatrix2D res = DoubleFactory2D.sparse.make(this.numCells, this.model.getNumScalarDecisionVariables());
		for (Entry<Integer, Cell> c : this.linearCoefs.entrySet())
			if (c.getValue().lCoefs != null)
				for (Entry<Integer, Double> v : c.getValue().lCoefs.entrySet())
					res.set(c.getKey(), v.getKey(), v.getValue());
		return res;
	}

	LinkedHashMap<Integer, Double> getScalarExpressionLinearCoefs(int cellId)
	{
		Cell c = this.linearCoefs.get(cellId);
		if (c == null) return new LinkedHashMap<Integer, Double>();
		return c.lCoefs;
	}

	LinkedHashMap<Integer, HashSet<Integer>> getNonZeroLinearCoefsCols()
	{
		LinkedHashMap<Integer, HashSet<Integer>> hm = new LinkedHashMap<Integer, HashSet<Integer>>();
		//if (this.linearCoefs == null) return hm;
		for (Entry<Integer, Cell> e : this.linearCoefs.entrySet())
		{
			if (e.getValue().lCoefs == null) continue;
			HashSet<Integer> dvs = new HashSet<Integer>(e.getValue().lCoefs.size());
			for (int dv : e.getValue().lCoefs.keySet())
				dvs.add(dv);
			hm.put(e.getKey(), dvs);
		}
		return hm;
	}

	int getNumVarIds(){ return this.model.getNumScalarDecisionVariables(); }

	int[] getSize(){ return this.size; }

	int getNumScalarExpressions(){ return this.numCells; }

	int getNumDim(){ return this.numDim; }

	OptimizationProblem getModel(){ return this.model; }

	LinkedHashMap<Integer, Double> getCellNonZeroLinearCoefs(int cellIndex)
	{
		Cell c = this.linearCoefs.get(cellIndex);
		return (c == null) ? null : c.lCoefs;
	}

	double getCellConstantCoef(int cellIndex)
	{
		if (this.linearCoefs == null) throw new JOMException("Affine expression operation error: Unexpected error");
		Cell c = this.linearCoefs.get(cellIndex);
		return (c == null) ? 0 : c.constantValue;
	}

	double[] getConstantCoefArray()
	{
		double[] res = new double[this.numCells];
		for (Entry<Integer, Cell> e : this.linearCoefs.entrySet())
			res[e.getKey()] = e.getValue().constantValue;
		return res;
	}

	int numNonZerosLinearCoefs(int cellIndex)
	{
		Cell c = this.linearCoefs.get(cellIndex);
		return (c == null) ? 0 : (c.lCoefs == null) ? 0 : c.lCoefs.size();
	}

	int numNonZerosLinearCoefs()
	{
		int res = 0;
		for (int c = 0; c < this.numCells; c++) res += numNonZerosLinearCoefs(c);
		return res;
	}

	void getNonZerosRowColVal(int[][] rows, int[][] cols, double[][] vals)
	{
		int numNonZerosLinearCoefs = 0;
		/* count number of non-zero coefficients */
		for (Cell c : this.linearCoefs.values()) if (c.lCoefs != null) numNonZerosLinearCoefs += c.lCoefs.size();

		int[] rowsArray = new int[numNonZerosLinearCoefs];
		int[] colsArray = new int[numNonZerosLinearCoefs];
		double[] valsArray = new double[numNonZerosLinearCoefs];
		int counter = 0;
		for (Entry<Integer, Cell> e : this.linearCoefs.entrySet())
		{
			final int cellId = e.getKey();
			final Cell c = e.getValue();
			if (c.lCoefs != null)
				for (Entry<Integer, Double> lc : c.lCoefs.entrySet())
				{
					rowsArray[counter] = cellId;
					colsArray[counter] = lc.getKey();
					valsArray[counter++] = lc.getValue();
				}
		}
		rows[0] = rowsArray;
		cols[0] = colsArray;
		vals[0] = valsArray;
	}

	void getNonZerosRowColValForXpressSolver(int[][] return_mnel, int[][] return_mrwind, int[][] return_mstart, double[][] return_dmatval , int numCols)
	{
		if (numDim != 2) throw new JOMException ("This call is allowed just for matrices");
		/* count number of non-zero coefficients */
		int numNonZerosLinearCoefs = 0;
		for (Cell c : this.linearCoefs.values()) if (c.lCoefs != null) numNonZerosLinearCoefs += c.lCoefs.size();
		
		int [] _mnel = new int [numCols]; // one element per column in the constraints matrix, with the number of nonzeros
		int [] _mrwind = new int [numNonZerosLinearCoefs]; // for each non-zero in the constraints matrix, the row
		int [] _mstart = new int [numCols]; // for each column, the index of the first element in mrwind and dmatval for the nonzeros of that column
        double [] _dmatval = new double [numNonZerosLinearCoefs];
        
        return_mnel [0] = _mnel;
        return_mrwind [0] = _mrwind;
        return_mstart [0] = _mstart;
        return_dmatval [0] = _dmatval;
        
        if (numNonZerosLinearCoefs == 0) return;
        
        /* First, create the _mnel array: the number of nonzero elements in each column */
		for (Cell c : this.linearCoefs.values())
			if (c != null)
				for (Integer col : c.lCoefs.keySet())
					_mnel [col] ++;
        /* Create the mstart array: the start of the first nonzero element */
		for (int col = 1 ; col < numCols ; col ++)
			_mstart[col] = _mstart[col-1] + _mnel[col-1];
		int [] numAlreadySetValuesThisColumn = new int [numCols];
		
        /* Create the mrwind and dmatval arrays */
		for (Entry<Integer, Cell> e : this.linearCoefs.entrySet())
		{
			final int row = e.getKey();
			final Cell c = e.getValue();
			if (c.lCoefs != null)
				for (Entry<Integer, Double> lc : c.lCoefs.entrySet())
				{
					final int col = lc.getKey();
					final int positionInDataArray = _mstart [col] + numAlreadySetValuesThisColumn [col];
					_mrwind [positionInDataArray] = row;
					_dmatval [positionInDataArray] = lc.getValue(); 
					numAlreadySetValuesThisColumn [col] ++;
				}
		}
	}

	
	double[] getCellLinearCoefsFull(int cellIndex)
	{
		double[] res = new double[model.getNumScalarDecisionVariables()];
		Cell c = this.linearCoefs.get(cellIndex);
		if (c == null) return res;
		if (c.lCoefs != null)
			for (Entry<Integer, Double> lc : c.lCoefs.entrySet())
				res[lc.getKey()] = lc.getValue();
		return res;
	}

	double getMinLinearCoef(int cellIndex)
	{
		Cell c = this.linearCoefs.get(cellIndex);
		if (c == null) return 0;
		if (c.lCoefs == null) return 0;
		double minimum = Double.MAX_VALUE;
		for (Entry<Integer, Double> lc : c.lCoefs.entrySet())
			minimum = (lc.getValue() < minimum) ? lc.getValue() : minimum;
		if (c.lCoefs.size() == this.model.getNumScalarDecisionVariables())
			return minimum;
		else
			return Math.min(0, minimum);
	}

	double getMaxLinearCoef(int cellIndex)
	{
		Cell c = this.linearCoefs.get(cellIndex);
		if (c == null) return 0;
		if (c.lCoefs == null) return 0;
		double maximum = -Double.MAX_VALUE;
		for (Entry<Integer, Double> lc : c.lCoefs.entrySet())
			maximum = (lc.getValue() > maximum) ? lc.getValue() : maximum;
		if (c.lCoefs.size() == this.model.getNumScalarDecisionVariables())
			return maximum;
		else
			return Math.max(0, maximum);
	}

	/* evluates the affine expression for the given value of the decision variables */
	DoubleMatrixND evaluate(final double[] valuesDVs)
	{
		double[] res = new double[this.numCells];
		for (Entry<Integer, Cell> e : this.linearCoefs.entrySet())
		{
			int cellIndex = e.getKey();
			Cell c = e.getValue();
			res[cellIndex] = c.evaluate(valuesDVs);
		}
		return new DoubleMatrixND(this.size, res);
	}

	//	static _INTERNAL_AffineExpressionCoefs operator_multiply_old (_INTERNAL_AffineExpressionCoefs a , _INTERNAL_AffineExpressionCoefs b)
	//	{
	//		if (!a.isConstant() && !b.isConstant()) throw new JOMException ("Affine expression operation error: Unexpected error");
	//		if ((a.numDim != 2) || (b.getNumDim() != 2)) throw new JOMException ("Affine expression operation error: Unexpected error");
	//		if (a.size [1] != b.size [0]) throw new JOMException ("Affine expression operation error: Unexpected error");
	//		final int aRows = a.size [0];
	//		final int aCols = a.size [1];
	//		final int bRows = b.size [0];
	//		final int bCols = b.size [1];
	//		final int resRows = aRows;
	//		final int resCols = bCols;
	//		final int resNumCells = resRows * resCols;
	//		final LinkedHashMap<Integer,Cell> aLC = a.linearCoefs;
	//		final LinkedHashMap<Integer,Cell> bLC = b.linearCoefs;
	//		final int [] resSize = new int [] { aRows , bCols };
	//
	//		LinkedHashMap<Integer,Cell> res_linearCoefs = new LinkedHashMap<Integer,Cell> ();
	//
	//		/* Arrange in an array, for each column of a, the non-zero cells */
	//		ArrayList<LinkedList<Entry<Integer,Cell>>> a_perColumn = new ArrayList<LinkedList<Entry<Integer,Cell>>> (aCols);
	//		for (int cont = 0 ; cont < aCols ; cont ++) a_perColumn.add (new LinkedList<Entry<Integer,Cell>> ());
	//		for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//			{ final int aCol = aCellEntry.getKey() / aRows; a_perColumn.get (aCol).add(aCellEntry);		}
	//
	//		for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//		{
	//			final int bCellId = bCellEntry.getKey();
	//			final Cell bCell = bCellEntry.getValue();
	//			final int bCol = bCellId / bRows;
	//			final int bRow = bCellId % bRows;
	//			for (Entry<Integer,Cell> aCellEntry : a_perColumn.get(bRow))
	//			{
	//				final int aCellId = aCellEntry.getKey();
	//				final Cell aCell = aCellEntry.getValue();
	//				final int aRow = aCellId % aRows;
	//				final int aCol = aCellId / aRows;
	//				final int resRow = aRow;
	//				final int resCol = bCol;
	//				final int resCellId = resRow + resRows * resCol;
	//				final Cell summand = aCell.deepCopy().multiply(bCell);
	//				if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//				Cell accumCell = res_linearCoefs.get(resCellId);
	//				if (accumCell == null)
	//					res_linearCoefs.put(resCellId, summand);
	//				else
	//					accumCell.sum(summand);
	//			}
	//		}
	//
	//		return new _INTERNAL_AffineExpressionCoefs (a.model , resSize , res_linearCoefs);
	//	}
	//
	//	static _INTERNAL_AffineExpressionCoefs function_matProdOld (_INTERNAL_AffineExpressionCoefs A , _INTERNAL_AffineExpressionCoefs B , boolean
	// dimensionsApplyToA , int dimensionRows0Index , int dimensionColumns0Index)
	//	{
	//		if (!A.isConstant() && !B.isConstant()) throw new JOMException ("Affine expression operation error: Unexpected error");
	//		//if ((A.numDim <= 2) && (B.getNumDim() <= 2)) throw new JOMException ("Affine expression operation error: matProd operation is
	// applicable if at least a matrix is larger than 2D");
	//		LinkedHashMap<Integer,Cell> globalRes_linearCoefs = new LinkedHashMap<Integer,Cell> ();
	//
	//		int [] globalResSize = null;
	//		final int [] sizeA = A.getSize ();
	//		final int [] sizeB = B.getSize ();
	//		if (dimensionsApplyToA)
	//		{
	//			globalResSize = Arrays.copyOf(sizeA , sizeA.length); globalResSize [dimensionColumns0Index] = sizeB [1];
	//		}
	//		else
	//		{
	//			globalResSize = Arrays.copyOf(sizeB , sizeB.length); globalResSize [dimensionRows0Index] = sizeA [0];
	//		}
	//		//IntMatrix1D auxCoordComputation = IntFactory1D.dense.make (globalResSize);
	//
	//		System.out.println ("globalResSize: " + Arrays.toString (globalResSize));
	//
	//		final int aRows = dimensionsApplyToA? dimensionRows0Index : sizeA [0];
	//		final int aCols = dimensionsApplyToA? dimensionColumns0Index : sizeA [1];
	//		final int bRows = dimensionsApplyToA? sizeB [0] : dimensionRows0Index;
	//		final int bCols = dimensionsApplyToA? sizeB [1] : dimensionColumns0Index;
	//
	//		final LinkedHashMap<Integer,Cell> aLC = A.linearCoefs;
	//		final LinkedHashMap<Integer,Cell> bLC = B.linearCoefs;
	//
	//		/* Arrange in an array, for each column of a, the non-zero cells */
	//		ArrayList<LinkedList<Entry<Integer,Cell>>> a_perColumn = new ArrayList<LinkedList<Entry<Integer,Cell>>> (aCols);
	//		ArrayList<LinkedList<Entry<Integer,Cell>>> b_perRow = new ArrayList<LinkedList<Entry<Integer,Cell>>> (bRows);
	//		if (dimensionsApplyToA)
	//		{
	//			// B is 2D => cache matrix B per rows
	//			for (int cont = 0 ; cont < bRows ; cont ++) b_perRow.add (new LinkedList<Entry<Integer,Cell>> ());
	//			for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//				{ final int bRow = bCellEntry.getKey() % bRows; b_perRow.get (bRow).add(bCellEntry);		}
	//			for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//			{
	//				final int aCellId = aCellEntry.getKey();
	//				final Cell aCell = aCellEntry.getValue();
	//				final IntMatrix1D subindexesInA = DoubleMatrixND.ind2sub(aCellId , sizeA);
	//				final int aCol = subindexesInA.get (dimensionColumns0Index);
	//				final int aRow = subindexesInA.get (dimensionRows0Index);
	//				for (Entry<Integer,Cell> bCellEntry : b_perRow.get(aCol))
	//				{
	//					final int bCellId = bCellEntry.getKey();
	//					final Cell bCell = bCellEntry.getValue();
	//					final int bRow = bCellId % bRows;
	//					final int bCol = bCellId / bRows;
	//					final int resRow = aRow;
	//					final int resCol = bCol;
	//					subindexesInA.set (dimensionRows0Index , resRow);
	//					subindexesInA.set (dimensionColumns0Index , resCol);
	//					final int globalResCellId = DoubleMatrixND.sub2ind(subindexesInA , globalResSize);
	//					final Cell summand = bCell.deepCopy().multiply(aCell);
	//					if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//					Cell accumCell = globalRes_linearCoefs.get(globalResCellId);
	//					if (accumCell == null)
	//						globalRes_linearCoefs.put(globalResCellId, summand);
	//					else
	//						accumCell.sum(summand);
	//				}
	//			}
	//		}
	//		else
	//		{
	//			// A is 2D => cache matrix A per columns
	//			for (int cont = 0 ; cont < aCols ; cont ++) a_perColumn.add (new LinkedList<Entry<Integer,Cell>> ());
	//			for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//				{ final int aCol = aCellEntry.getKey() / aRows; a_perColumn.get (aCol).add(aCellEntry);		}
	//			for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//			{
	//				final int bCellId = bCellEntry.getKey();
	//				final Cell bCell = bCellEntry.getValue();
	//				final IntMatrix1D subindexesInB = DoubleMatrixND.ind2sub(bCellId , sizeB);
	//				final int bCol = subindexesInB.get (dimensionColumns0Index);
	//				final int bRow = subindexesInB.get (dimensionRows0Index);
	//				for (Entry<Integer,Cell> aCellEntry : a_perColumn.get(bRow))
	//				{
	//					final int aCellId = aCellEntry.getKey();
	//					final Cell aCell = aCellEntry.getValue();
	//					final int aRow = aCellId % aRows;
	//					final int aCol = aCellId / aRows;
	//					final int resRow = aRow;
	//					final int resCol = bCol;
	//					subindexesInB.set (dimensionRows0Index , resRow);
	//					subindexesInB.set (dimensionColumns0Index , resCol);
	//					final int globalResCellId = DoubleMatrixND.sub2ind(subindexesInB , globalResSize);
	//					final Cell summand = aCell.deepCopy().multiply(bCell);
	//					if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//					Cell accumCell = globalRes_linearCoefs.get(globalResCellId);
	//					if (accumCell == null)
	//						globalRes_linearCoefs.put(globalResCellId, summand);
	//					else
	//						accumCell.sum(summand);
	//				}
	//			}
	//		}
	//
	//		return new _INTERNAL_AffineExpressionCoefs (A.model , globalResSize , globalRes_linearCoefs);
	//	}
	//
	//	static _INTERNAL_AffineExpressionCoefs function_matProd (_INTERNAL_AffineExpressionCoefs A , _INTERNAL_AffineExpressionCoefs B , boolean
	// dimensionsApplyToA , int dimensionRows0Index , int dimensionColumns0Index)
	//	{
	//		if (!A.isConstant() && !B.isConstant()) throw new JOMException ("Affine expression operation error: Unexpected error");
	//		//if ((A.numDim <= 2) && (B.getNumDim() <= 2)) throw new JOMException ("Affine expression operation error: matProd operation is
	// applicable if at least a matrix is larger than 2D");
	//		LinkedHashMap<Integer,Cell> globalRes_linearCoefs = new LinkedHashMap<Integer,Cell> ();
	//
	//		int [] globalResSize = null;
	//		final int [] sizeA = A.getSize ();
	//		final int [] sizeB = B.getSize ();
	//		if (dimensionsApplyToA)
	//		{
	//			globalResSize = Arrays.copyOf(sizeA , sizeA.length); globalResSize [dimensionColumns0Index] = sizeB [1];
	//		}
	//		else
	//		{
	//			globalResSize = Arrays.copyOf(sizeB , sizeB.length); globalResSize [dimensionRows0Index] = sizeA [0];
	//		}
	//		//IntMatrix1D auxCoordComputation = IntFactory1D.dense.make (globalResSize);
	//
	//		System.out.println ("globalResSize: " + Arrays.toString (globalResSize));
	//
	//		final int aRows = dimensionsApplyToA? dimensionRows0Index : sizeA [0];
	//		final int aCols = dimensionsApplyToA? dimensionColumns0Index : sizeA [1];
	//		final int bRows = dimensionsApplyToA? sizeB [0] : dimensionRows0Index;
	//		final int bCols = dimensionsApplyToA? sizeB [1] : dimensionColumns0Index;
	//
	//		final LinkedHashMap<Integer,Cell> aLC = A.linearCoefs;
	//		final LinkedHashMap<Integer,Cell> bLC = B.linearCoefs;
	//
	//		/* Arrange in an array, for each column of a, the non-zero cells */
	//		ArrayList<LinkedList<Entry<Integer,Cell>>> a_perColumn = new ArrayList<LinkedList<Entry<Integer,Cell>>> (aCols);
	//		ArrayList<LinkedList<Entry<Integer,Cell>>> b_perRow = new ArrayList<LinkedList<Entry<Integer,Cell>>> (bRows);
	//		if (dimensionsApplyToA)
	//		{
	//			// B is 2D => cache matrix B per rows
	//			for (int cont = 0 ; cont < bRows ; cont ++) b_perRow.add (new LinkedList<Entry<Integer,Cell>> ());
	//			for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//				{ final int bRow = bCellEntry.getKey() % bRows; b_perRow.get (bRow).add(bCellEntry);		}
	//			for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//			{
	//				final int aCellId = aCellEntry.getKey();
	//				final Cell aCell = aCellEntry.getValue();
	//				final IntMatrix1D subindexesInA = DoubleMatrixND.ind2sub(aCellId , sizeA);
	//				final int aCol = subindexesInA.get (dimensionColumns0Index);
	//				final int aRow = subindexesInA.get (dimensionRows0Index);
	//				for (Entry<Integer,Cell> bCellEntry : b_perRow.get(aCol))
	//				{
	//					final int bCellId = bCellEntry.getKey();
	//					final Cell bCell = bCellEntry.getValue();
	//					final int bRow = bCellId % bRows;
	//					final int bCol = bCellId / bRows;
	//					final int resRow = aRow;
	//					final int resCol = bCol;
	//					subindexesInA.set (dimensionRows0Index , resRow);
	//					subindexesInA.set (dimensionColumns0Index , resCol);
	//					final int globalResCellId = DoubleMatrixND.sub2ind(subindexesInA , globalResSize);
	//					final Cell summand = bCell.deepCopy().multiply(aCell);
	//					if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//					Cell accumCell = globalRes_linearCoefs.get(globalResCellId);
	//					if (accumCell == null)
	//						globalRes_linearCoefs.put(globalResCellId, summand);
	//					else
	//						accumCell.sum(summand);
	//				}
	//			}
	//		}
	//		else
	//		{
	//			// A is 2D => cache matrix A per columns
	//			for (int cont = 0 ; cont < aCols ; cont ++) a_perColumn.add (new LinkedList<Entry<Integer,Cell>> ());
	//			for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//				{ final int aCol = aCellEntry.getKey() / aRows; a_perColumn.get (aCol).add(aCellEntry);		}
	//			for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//			{
	//				final int bCellId = bCellEntry.getKey();
	//				final Cell bCell = bCellEntry.getValue();
	//				final IntMatrix1D subindexesInB = DoubleMatrixND.ind2sub(bCellId , sizeB);
	//				final int bCol = subindexesInB.get (dimensionColumns0Index);
	//				final int bRow = subindexesInB.get (dimensionRows0Index);
	//				for (Entry<Integer,Cell> aCellEntry : a_perColumn.get(bRow))
	//				{
	//					final int aCellId = aCellEntry.getKey();
	//					final Cell aCell = aCellEntry.getValue();
	//					final int aRow = aCellId % aRows;
	//					final int aCol = aCellId / aRows;
	//					final int resRow = aRow;
	//					final int resCol = bCol;
	//					subindexesInB.set (dimensionRows0Index , resRow);
	//					subindexesInB.set (dimensionColumns0Index , resCol);
	//					final int globalResCellId = DoubleMatrixND.sub2ind(subindexesInB , globalResSize);
	//					final Cell summand = aCell.deepCopy().multiply(bCell);
	//					if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//					Cell accumCell = globalRes_linearCoefs.get(globalResCellId);
	//					if (accumCell == null)
	//						globalRes_linearCoefs.put(globalResCellId, summand);
	//					else
	//						accumCell.sum(summand);
	//				}
	//			}
	//		}
	//
	//		return new _INTERNAL_AffineExpressionCoefs (A.model , globalResSize , globalRes_linearCoefs);
	//	}
	//
	//	static _INTERNAL_AffineExpressionCoefs function_matProd (_INTERNAL_AffineExpressionCoefs A , _INTERNAL_AffineExpressionCoefs B)
	//	{
	//		if (!A.isConstant() && !B.isConstant()) throw new JOMException ("Affine expression operation error: Unexpected error");
	//		//if ((A.numDim <= 2) && (B.getNumDim() <= 2)) throw new JOMException ("Affine expression operation error: matProd operation is
	// applicable if at least a matrix is larger than 2D");
	////		LinkedHashMap<Integer,Cell> globalRes_linearCoefs = new LinkedHashMap<Integer,Cell> ();
	//		final int aRows = A.size [0];
	//		final int aCols = A.size [1];
	//		final int bRows = B.size [0];
	//		final int bCols = B.size [1];
	//
	//		int [] globalResSize = Arrays.copyOf(B.size, B.size.length); globalResSize [0] = aRows;
	//		//IntMatrix1D auxCoordComputation = IntFactory1D.dense.make (globalResSize);
	//
	//		System.out.println ("globalResSize: " + Arrays.toString (globalResSize));
	//
	//		final LinkedHashMap<Integer,Cell> aLC = A.linearCoefs;
	//		final LinkedHashMap<Integer,Cell> bLC = B.linearCoefs;
	//
	//		/* Arrange in an array, for each column of a, the non-zero cells */
	//		LinkedHashMap<Integer,Cell> res_linearCoefs = new LinkedHashMap<Integer,Cell> ();
	//
	//		/* Arrange in an array, for each column of a, the non-zero cells */
	//		ArrayList<LinkedList<Entry<Integer,Cell>>> a_perColumn = new ArrayList<LinkedList<Entry<Integer,Cell>>> (aCols);
	//		for (int cont = 0 ; cont < aCols ; cont ++) a_perColumn.add (new LinkedList<Entry<Integer,Cell>> ());
	//		for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//			{ final int aCol = aCellEntry.getKey() / aRows; a_perColumn.get (aCol).add(aCellEntry);		}
	//
	//		for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//		{
	//			final int bCellId = bCellEntry.getKey();
	//			final Cell bCell = bCellEntry.getValue();
	//			final int bCol = bCellId / bRows;
	//			final int bRow = bCellId % bRows;
	//			for (Entry<Integer,Cell> aCellEntry : a_perColumn.get(bRow))
	//			{
	//				final int aCellId = aCellEntry.getKey();
	//				final Cell aCell = aCellEntry.getValue();
	//				final int aRow = aCellId % aRows;
	//				final int aCol = aCellId / aRows;
	//				final int resRow = aRow;
	//				final int resCol = bCol;
	//				final int resCellId = resRow + resRows * resCol;
	//				final Cell summand = aCell.deepCopy().multiply(bCell);
	//				if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//				Cell accumCell = res_linearCoefs.get(resCellId);
	//				if (accumCell == null)
	//					res_linearCoefs.put(resCellId, summand);
	//				else
	//					accumCell.sum(summand);
	//			}
	//		}
	//
	//
	//
	//
	//
	//		if (dimensionsApplyToA)
	//		{
	//			// B is 2D => cache matrix B per rows
	//			for (int cont = 0 ; cont < bRows ; cont ++) b_perRow.add (new LinkedList<Entry<Integer,Cell>> ());
	//			for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//				{ final int bRow = bCellEntry.getKey() % bRows; b_perRow.get (bRow).add(bCellEntry);		}
	//			for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//			{
	//				final int aCellId = aCellEntry.getKey();
	//				final Cell aCell = aCellEntry.getValue();
	//				final IntMatrix1D subindexesInA = DoubleMatrixND.ind2sub(aCellId , sizeA);
	//				final int aCol = subindexesInA.get (dimensionColumns0Index);
	//				final int aRow = subindexesInA.get (dimensionRows0Index);
	//				for (Entry<Integer,Cell> bCellEntry : b_perRow.get(aCol))
	//				{
	//					final int bCellId = bCellEntry.getKey();
	//					final Cell bCell = bCellEntry.getValue();
	//					final int bRow = bCellId % bRows;
	//					final int bCol = bCellId / bRows;
	//					final int resRow = aRow;
	//					final int resCol = bCol;
	//					subindexesInA.set (dimensionRows0Index , resRow);
	//					subindexesInA.set (dimensionColumns0Index , resCol);
	//					final int globalResCellId = DoubleMatrixND.sub2ind(subindexesInA , globalResSize);
	//					final Cell summand = bCell.deepCopy().multiply(aCell);
	//					if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//					Cell accumCell = globalRes_linearCoefs.get(globalResCellId);
	//					if (accumCell == null)
	//						globalRes_linearCoefs.put(globalResCellId, summand);
	//					else
	//						accumCell.sum(summand);
	//				}
	//			}
	//		}
	//		else
	//		{
	//			// A is 2D => cache matrix A per columns
	//			for (int cont = 0 ; cont < aCols ; cont ++) a_perColumn.add (new LinkedList<Entry<Integer,Cell>> ());
	//			for (Entry<Integer,Cell> aCellEntry : aLC.entrySet())
	//				{ final int aCol = aCellEntry.getKey() / aRows; a_perColumn.get (aCol).add(aCellEntry);		}
	//			for (Entry<Integer,Cell> bCellEntry : bLC.entrySet())
	//			{
	//				final int bCellId = bCellEntry.getKey();
	//				final Cell bCell = bCellEntry.getValue();
	//				final IntMatrix1D subindexesInB = DoubleMatrixND.ind2sub(bCellId , sizeB);
	//				final int bCol = subindexesInB.get (dimensionColumns0Index);
	//				final int bRow = subindexesInB.get (dimensionRows0Index);
	//				for (Entry<Integer,Cell> aCellEntry : a_perColumn.get(bRow))
	//				{
	//					final int aCellId = aCellEntry.getKey();
	//					final Cell aCell = aCellEntry.getValue();
	//					final int aRow = aCellId % aRows;
	//					final int aCol = aCellId / aRows;
	//					final int resRow = aRow;
	//					final int resCol = bCol;
	//					subindexesInB.set (dimensionRows0Index , resRow);
	//					subindexesInB.set (dimensionColumns0Index , resCol);
	//					final int globalResCellId = DoubleMatrixND.sub2ind(subindexesInB , globalResSize);
	//					final Cell summand = aCell.deepCopy().multiply(bCell);
	//					if ((summand.constantValue == 0) && (summand.lCoefs == null)) continue;
	//					Cell accumCell = globalRes_linearCoefs.get(globalResCellId);
	//					if (accumCell == null)
	//						globalRes_linearCoefs.put(globalResCellId, summand);
	//					else
	//						accumCell.sum(summand);
	//				}
	//			}
	//		}
	//
	//		return new _INTERNAL_AffineExpressionCoefs (A.model , globalResSize , globalRes_linearCoefs);
	//	}
	//

	String toStringLinearScalarExpression(int cellIndex)
	{
		Cell c = this.linearCoefs.get(cellIndex);
		if (c == null) throw new JOMException("Affine expression operation error: Unexpected error");

		double constant = c.constantValue;
		String s = "" + constant;
		if (c.lCoefs == null) return s;
		for (Entry<Integer, Double> lc : c.lCoefs.entrySet())
		{
			int varId = lc.getKey();
			double coef = lc.getValue();
			_INTERNAL_DecisionVariableArray varInfo = model.getVarInfo(varId);
			String varName = varInfo.getName();
			int indexInVariable = varId - varInfo.getFirstVarId();
			IntMatrix1D subindexesInVariable = IntMatrixND.ind2sub(indexInVariable, varInfo.getSize());
			s = s + " + " + coef + "*" + varName + Arrays.toString(subindexesInVariable.toArray());
		}
		return s;
	}

	_INTERNAL_AffineExpressionCoefs function_appendRows(_INTERNAL_AffineExpressionCoefs a)
	{
		if ((this.numDim != 2) || (a.getNumDim() != 2))
			throw new JOMException("Affine expression operation error: Append rows function only applies to 2-D expressions, in the form of column "
					+ "vector");
		if ((this.size[1] != 1) || (a.getSize()[1] != 1))
			throw new JOMException("Affine expression operation error: Append rows function only applies to 2-D expressions, in the form of column "
					+ "vector");

		final int oldNumCellsThis = this.numCells;
		this.resize(new int[]{this.numCells + a.getSize()[0], 1});
		for (Entry<Integer, Cell> e : a.linearCoefs.entrySet())
			this.linearCoefs.put(e.getKey() + oldNumCellsThis, e.getValue());
		return this;
	}

	_INTERNAL_AffineExpressionCoefs function_diag()
	{
		if (this.getNumDim() != 2)
			throw new JOMException("Affine expression operation error: Diag function is for 2-D expressions in the form of row or column vector");
		if ((this.getSize()[0] != 1) && (this.getSize()[1] != 1))
			throw new JOMException("Affine expression operation error: Diag function is for 2-D expressions in the form of row or column vector");
		int numElemDiag = this.getNumScalarExpressions();
		int resNumCells = numElemDiag * numElemDiag;
		LinkedHashMap<Integer, Cell> newLinearCoefs = new LinkedHashMap<Integer, Cell>();
		for (Entry<Integer, Cell> e : this.linearCoefs.entrySet())
		{
			final int oldCellIndex = e.getKey();
			final Cell c = e.getValue();
			newLinearCoefs.put(oldCellIndex * (numElemDiag + 1), c);
		}
		this.resize(new int[]{numElemDiag, numElemDiag});
		this.linearCoefs = newLinearCoefs;
		return this;
	}

	_INTERNAL_AffineExpressionCoefs function_onesOrZeros(double numberOneOrZero)
	{
		if (this.getNumDim() != 2)
			throw new JOMException("Affine expression operation error: Ones and zeros functions should receive a row or column vector of constants "
					+ "indicanting the array size");
		if ((this.getSize()[0] != 1) && (this.getSize()[1] != 1))
			throw new JOMException("Affine expression operation error: Ones and zeros functionsshould receive a row or column vector of constants "
					+ "indicanting the array size");
		if (!this.isConstant())
			throw new JOMException("Affine expression operation error: Ones and zeros functions should receive a row or column vector of constants "
					+ "indicanting the array size");

		int[] size_i = new int[this.getNumScalarExpressions()];
		for (int cont = 0; cont < size_i.length; cont++) size_i[cont] = (int) this.linearCoefs.get(cont).constantValue;
		this.resize(size_i);
		this.linearCoefs = new LinkedHashMap<Integer, Cell>();

		/* If zeros, just return */
		if (numberOneOrZero == 0) return this;

		/* If ones, create cells */
		for (int cellId = 0; cellId < this.numCells; cellId++)
			this.linearCoefs.put(cellId, new Cell(1));
		return this;
	}

	_INTERNAL_AffineExpressionCoefs function_eye(int a, int b)
	{
		final int numElementsDiagonal = (a < b) ? a : b;
		this.resize(new int[]{a, b});
		this.linearCoefs = new LinkedHashMap<Integer, Cell>();
		for (int cont = 0; cont < numElementsDiagonal; cont++)
			this.linearCoefs.put(cont * (numElementsDiagonal + 1), new Cell(1));
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_transpose()
	{
		if (this.numDim != 2) throw new JOMException("Affine expression operation error: Unexpected error");
		final int oldRows = this.size[0];
		final int oldCols = this.size[1];
		this.resize(new int[]{this.size[1], this.size[0]});
		LinkedHashMap<Integer, Cell> newLinearCoefs = new LinkedHashMap<Integer, Cell>(); // first key is the cell index, second key is the variable
		// index

		for (Entry<Integer, Cell> cellEntry : this.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			final Cell cell = cellEntry.getValue();
			final int oldCol = cellId / oldRows;
			final int oldRow = cellId % oldRows;
			newLinearCoefs.put(oldCol + oldCols * oldRow, cell);
		}
		this.linearCoefs = newLinearCoefs;
		return this;
	}

	_INTERNAL_AffineExpressionCoefs function_permute(int[] permutedDimensionsIndex0)
	{
		if (this.numDim != permutedDimensionsIndex0.length) throw new JOMException("Affine expression operation error: Unexpected error");
		final int[] oldSize = Arrays.copyOf(this.size, size.length);
		final int[] newSize = new int[permutedDimensionsIndex0.length];
		for (int dim = 0; dim < newSize.length; dim++) newSize[dim] = this.size[permutedDimensionsIndex0[dim]];
		this.resize(newSize);
		LinkedHashMap<Integer, Cell> newLinearCoefs = new LinkedHashMap<Integer, Cell>(); // first key is the cell index, second key is the variable
		// index

		for (Entry<Integer, Cell> cellEntry : this.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			final Cell cell = cellEntry.getValue();
			final IntMatrix1D oldSubindexes = DoubleMatrixND.ind2sub(cellId, oldSize);
			final IntMatrix1D newSubindexes = IntFactory1D.dense.make(this.numDim);
			for (int dim = 0; dim < newSize.length; dim++) newSubindexes.set(dim, oldSubindexes.get(permutedDimensionsIndex0[dim]));
			final int newIndex = DoubleMatrixND.sub2ind(newSubindexes, newSize);
			newLinearCoefs.put(newIndex, cell);
		}
		this.linearCoefs = newLinearCoefs;
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_appendColumns(_INTERNAL_AffineExpressionCoefs a)
	{
		if ((this.numDim != 2) || (a.getNumDim() != 2))
			throw new JOMException("Affine expression operation error: Append columns function only applies to 2-D expressions");
		if (this.size[0] != a.getSize()[0])
			throw new JOMException("Affine expression operation error: Append columns function only applies to 2-D expressions, with the same number"
					+ " of rows");
		final int nElem_this = this.numCells;
		for (Entry<Integer, Cell> cellEntry : a.linearCoefs.entrySet())
			this.linearCoefs.put(nElem_this + cellEntry.getKey(), cellEntry.getValue());
		this.resize(new int[]{this.size[0], this.size[1] + a.getSize()[1]});
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_appendRows(_INTERNAL_AffineExpressionCoefs a)
	{
		if ((this.numDim != 2) || (a.getNumDim() != 2))
			throw new JOMException("Affine expression operation error: Append rows function only applies to 2-D expressions");
		if (this.size[1] != a.getSize()[1])
			throw new JOMException("Affine expression operation error: Append rows function only applies to 2-D expressions, with the same number of"
					+ " columns");

		final int resNumCols = this.size[1];
		final int numRows_this = this.size[0];
		final int numRows_a = a.getSize()[0];
		final int resNumRows = numRows_this + numRows_a;

		LinkedHashMap<Integer, Cell> newLinearCoefs = new LinkedHashMap<Integer, Cell>(); // first key is the cell index, second key is the variable
		// index
		for (Entry<Integer, Cell> cellEntry : this.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			final Cell cell = cellEntry.getValue();
			final int thisCol = cellId / numRows_this;
			final int thisRow = cellId % numRows_this;
			newLinearCoefs.put(thisRow + thisCol * resNumRows, cell);
		}
		for (Entry<Integer, Cell> cellEntry : a.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			final Cell cell = cellEntry.getValue();
			final int aCol = cellId / numRows_a;
			final int aRow = cellId % numRows_a;
			newLinearCoefs.put((aRow + numRows_this) + aCol * resNumRows, cell);
		}
		this.linearCoefs = newLinearCoefs;
		this.resize(new int[]{resNumRows, resNumCols});
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_neg()
	{
		for (Cell c : this.linearCoefs.values())
			c.neg();
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_e2eMultiply(_INTERNAL_AffineExpressionCoefs a, final boolean iAmConstant, final boolean aIsConstant)
	{
		if (!Arrays.equals(this.size, a.getSize()))
			throw new JOMException("Affine expression operation error: 2e2 multiply function. Wrong size of array");
		if (!iAmConstant && !aIsConstant) throw new JOMException("Affine expression operation error: 2e2 multiply function. Result is not linear");

		LinkedList<Integer> this_cellIdsToRemove = new LinkedList<Integer>();
		for (Entry<Integer, Cell> cellEntry : this.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			Cell cell_this = cellEntry.getValue();
			Cell cell_a = a.linearCoefs.get(cellId);
			if (cell_a == null) this_cellIdsToRemove.add(cellId);
			else cell_this.multiply(cell_a);
		}
		for (Integer idToRemove : this_cellIdsToRemove) this.linearCoefs.remove(idToRemove);
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_e2eDivide(_INTERNAL_AffineExpressionCoefs a)
	{
		if (!Arrays.equals(this.size, a.getSize()))
			throw new JOMException("Affine expression operation error: 2e2 divide function. Wrong size of array");
		if (!a.isConstant())
			throw new JOMException("Affine expression operation error: 2e2 divide function. Denominator must be a constant to produce a linear "
					+ "expression");
		if (a.linearCoefs.size() != this.numCells) throw new JOMException("Affine expression operation error: 2e2 divide function. Divide by zero");

		for (Entry<Integer, Cell> cellEntry : this.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			Cell cell_this = cellEntry.getValue();
			Cell cell_a = a.linearCoefs.get(cellId);
			if (cell_a.constantValue == 0) throw new JOMException("Affine expression operation error: 2e2 divide function. Divide by zero");
			cell_this.multiplyConstant(1 / cell_a.constantValue);
		}
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_plus(_INTERNAL_AffineExpressionCoefs a) // antiguo plusMinus
	{
		for (Entry<Integer, Cell> cellEntry : a.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			Cell cell_a = cellEntry.getValue();
			Cell cell_this = this.linearCoefs.get(cellId);
			if (cell_this == null) this.linearCoefs.put(cellId, cell_a);
			else cell_this.sum(cell_a);
		}
		return this;
	}

	_INTERNAL_AffineExpressionCoefs operator_substract(_INTERNAL_AffineExpressionCoefs a) // antiguo plusMinus
	{
		for (Entry<Integer, Cell> cellEntry : a.linearCoefs.entrySet())
		{
			final int cellId = cellEntry.getKey();
			Cell cell_a = cellEntry.getValue();
			Cell cell_this = this.linearCoefs.get(cellId);
			if (cell_this == null) this.linearCoefs.put(cellId, cell_a.neg());
			else cell_this.substract(cell_a);
		}
		return this;
	}

	_INTERNAL_AffineExpressionCoefs function_repeatExpression(int[] newSize)
	{
		if (this.numCells != 1)
			throw new JOMException("Affine expression operation error: Repeat expression function. Repeat function is for scalar expressions");
		Cell c = this.linearCoefs.get(0);
		if (c == null) c = new Cell(); // this means a constant zero

		this.resize(newSize);
		for (int cont = 1; cont < numCells; cont++)
			this.linearCoefs.put(cont, c.deepCopy());
		return this;
	}

	_INTERNAL_AffineExpressionCoefs function_repeatDiagonalExpression(int rowsAndColumns)
	{
		if (this.numCells != 1)
			throw new JOMException("Affine expression operation error: Repeat diagonal expression function. Repeat function is for scalar "
					+ "expressions");
		Cell c = this.linearCoefs.get(0);
		if (c == null) throw new JOMException("Affine expression operation error: Repeat expression function. Unexpected error");
		this.resize(new int[]{rowsAndColumns, rowsAndColumns});
		for (int cont = 1; cont < rowsAndColumns; cont++)
			this.linearCoefs.put(cont * (rowsAndColumns + 1), c.deepCopy());
		return this;
	}

	_INTERNAL_AffineExpressionCoefs function_sum(int dimensionInt, int[] newSize)
	{
		LinkedHashMap<Integer, Cell> newLinearCoefs = new LinkedHashMap<Integer, Cell>(); // first key is the cell index, second key is the variable
		// index

		/* sum all dimensions producing a scalar  */
		if (dimensionInt == -1)
		{
			this.resize(newSize);
			Cell accumCell = new Cell();
			for (Cell c : this.linearCoefs.values())
				accumCell.sum(c);
			newLinearCoefs.put(0, accumCell);
			this.linearCoefs = newLinearCoefs;
			return this;
		}

		/* Sum in a dimension */
		final int[] oldSize = this.size;
		final int[] blocksOldSize = computeBlockSizes(oldSize);
		final int[] blocksNewSize = computeBlockSizes(newSize);

		for (Entry<Integer, Cell> cellEntry : this.linearCoefs.entrySet())
		{
			final int oldCellId = cellEntry.getKey();
			final Cell cellToSum = cellEntry.getValue();
			final int[] oldCoordinates = ind2sub(oldCellId, oldSize, blocksOldSize);
			/* compute the index of the new cell */
			int newCellIndex = 0;
			for (int cont = 0; cont <= dimensionInt - 1; cont++) newCellIndex += blocksNewSize[cont] * oldCoordinates[cont];
			for (int cont = dimensionInt + 1; cont < oldSize.length; cont++) newCellIndex += blocksNewSize[cont - 1] * oldCoordinates[cont];
			/* Accumulate in that cell */
			Cell accumCell = newLinearCoefs.get(newCellIndex);
			if (accumCell == null)
				newLinearCoefs.put(newCellIndex, cellToSum);
			else
				accumCell.sum(cellToSum);
		}

		this.resize(newSize);
		this.linearCoefs = newLinearCoefs;
		return this;
	}

	final _INTERNAL_AffineExpressionCoefs reshape(int[] newSize)
	{
		if (prod(newSize) != prod(this.size))
			throw new JOMException("Affine expression operation error: Reshape expression function. Wrong array size. Original size: " + Arrays.toString(this.size) + ", new size: " + Arrays.toString(newSize));
		this.resize(Arrays.copyOf(newSize, newSize.length));
		return this;
	}

	boolean isConstant()
	{
		for (Cell c : this.linearCoefs.values())
			if (c.lCoefs != null) return false;
		return true;
	}

	private void resize(int[] size)
	{
		this.size = size;
		this.numDim = size.length;
		this.numCells = 1;
		for (int c : size) this.numCells *= c;
	}

	private class Cell
	{
		double constantValue;
		LinkedHashMap<Integer, Double> lCoefs;

		Cell()
		{
			this.constantValue = 0;
			this.lCoefs = null;
		}

		Cell(double val)
		{
			this.constantValue = val;
			this.lCoefs = null;
		}

		Cell(LinkedHashMap<Integer, Double> lCoefs)
		{
			this.constantValue = 0;
			this.lCoefs = lCoefs;
		}

		Cell(double constantValue, LinkedHashMap<Integer, Double> lCoefs)
		{
			this.constantValue = constantValue;
			this.lCoefs = lCoefs;
		}

		double evaluate(double[] dvValues)
		{
			double res = this.constantValue;
			if (this.lCoefs == null) return res;
			for (Entry<Integer, Double> lc : this.lCoefs.entrySet()) res += lc.getValue() * dvValues[lc.getKey()];
			return res;
		}

		Cell sum(Cell c)
		{
			this.constantValue += c.constantValue;
			if (c.lCoefs == null) return this;
			if (this.lCoefs == null)
			{
				this.lCoefs = c.lCoefs;
				return this;
			}
			for (Entry<Integer, Double> lcC : c.lCoefs.entrySet())
			{
				Double thisLcVal = this.lCoefs.get(lcC.getKey());
				if (thisLcVal == null) this.lCoefs.put(lcC.getKey(), lcC.getValue());
				else this.lCoefs.put(lcC.getKey(), thisLcVal + lcC.getValue());
			}
			return this;
		}

		Cell substract(Cell c)
		{
			this.constantValue -= c.constantValue;
			if (c.lCoefs == null) return this;
			if (this.lCoefs == null){ this.lCoefs = new LinkedHashMap<Integer, Double>(); }
			for (Entry<Integer, Double> lcC : c.lCoefs.entrySet())
			{
				Double thisLcVal = this.lCoefs.get(lcC.getKey());
				if (thisLcVal == null) this.lCoefs.put(lcC.getKey(), -lcC.getValue());
				else this.lCoefs.put(lcC.getKey(), thisLcVal - lcC.getValue());
			}
			return this;
		}

		Cell multiplyConstant(double c)
		{
			if (c == 0)
			{
				this.constantValue = 0;
				this.lCoefs = null;
				return this;
			}
			this.constantValue *= c;
			if (this.lCoefs == null) return this;
			for (Entry<Integer, Double> e : this.lCoefs.entrySet())
				e.setValue(c * e.getValue());
			return this;
		}

		Cell multiply(Cell c) // does not modify c
		{
			if (c == null)
			{
				this.constantValue = 0;
				this.lCoefs = null;
				return this;
			}
			if ((c.lCoefs != null) && (this.lCoefs != null))
				throw new JOMException("Affine expression operation error: multiplication of two non constant linear expressions is not a linear "
						+ "expression");
			if (c.lCoefs == null)
			{
				multiplyConstant(c.constantValue);
				return this;
			} // c is a constant
			/* c is not constant, i am constant */
			this.lCoefs = new LinkedHashMap<Integer, Double>();
			if (this.constantValue == 0) return this;
			for (Entry<Integer, Double> e : c.lCoefs.entrySet())
				this.lCoefs.put(e.getKey(), e.getValue() * this.constantValue);
			this.constantValue *= c.constantValue;
			return this;
		}

		Cell neg()
		{
			this.constantValue = -this.constantValue;
			if (this.lCoefs == null) return this;
			LinkedHashMap<Integer, Double> res = new LinkedHashMap<Integer, Double>();
			for (Entry<Integer, Double> e : this.lCoefs.entrySet())
				res.put(e.getKey(), -e.getValue());
			this.lCoefs = res;
			return this;
		}

		Cell deepCopy()
		{
			if (this.lCoefs == null) return new Cell(this.constantValue);
			else
			{
				LinkedHashMap<Integer, Double> res = new LinkedHashMap<Integer, Double>();
				for (Entry<Integer, Double> e : this.lCoefs.entrySet()) res.put(e.getKey(), e.getValue());
				return new Cell(this.constantValue, res);
			}
		}

		@Override
        public String toString()
		{
			String s = "" + this.constantValue;
			if (this.lCoefs != null) for (Entry<Integer, Double> e : this.lCoefs.entrySet())
				s = s + " + " + e.getValue() + "*" + model.getVarInfo(e.getKey()).getName() + "(" + (e.getKey() - model.getVarInfo(e.getKey())
						.getFirstVarId()) + ")";
			return s;
		}

		;
	}

}
