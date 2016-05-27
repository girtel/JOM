/*******************************************************************************
 * Copyright (c) 2015 Pablo Pavon Mariño.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * 
 * Contributors:
 *     Pablo Pavon Mariño - initial API and implementation
 ******************************************************************************/



 




package com.jom;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.regex.Pattern;

import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;

import com.jom.javaluator.AbstractEvaluator;
import com.jom.javaluator.BracketPair;
import com.jom.javaluator.Constant;
import com.jom.javaluator.Function;
import com.jom.javaluator.Operator;
import com.jom.javaluator.Parameters;

class _INTERNAL_ExpressionParser extends AbstractEvaluator<Expression> 
{
	private final HashMap<String, _INTERNAL_DecisionVariableArray> decisionVariables;
	private final HashMap<String, DoubleMatrixND> inputParameters;
	private final OptimizationProblem model;
	private final static BracketPair bracketPairs = BracketPair.PARENTHESES;
	private final static HashMap<String, Object []> listFunctions; // min num operans, max num operamds , file name, expand scalar->array
	private final static HashMap<String, Object []> listOperators; // symbol, num operans, associativity, priority, class to load, expand scalar->array
	
	static /* Initialize the operators information */
	{
		/* Define operators. Key in HashMap = symbol + " " + num operands */
		listOperators = new HashMap<String, Object []> ();
		listOperators.put(": 2", new Object [] {":", 2, Operator.Associativity.LEFT, 1, "com.jom._OPERATOR_CONTIGUOUSRANGE" , false}); // assignment
		listOperators.put("[ 1", new Object [] {"[", 1, Operator.Associativity.RIGHT, 4, "com.jom._OPERATOR_OPENARRAY" , false}); // e2eMultiply
		listOperators.put("; 2", new Object [] {";", 2, Operator.Associativity.LEFT, 4, "com.jom._OPERATOR_APPENDCOLUMNS" , false}); // e2eMultiply
		listOperators.put(";; 2", new Object [] {";;", 2, Operator.Associativity.LEFT, 4, "com.jom._OPERATOR_APPENDROWS" , false}); // e2eMultiply

		listOperators.put("+ 2", new Object [] {"+", 2, Operator.Associativity.LEFT, 1, "com.jom._OPERATOR_PLUS" , true}); // plus
		listOperators.put("- 2", new Object [] {"-", 2, Operator.Associativity.LEFT, 1, "com.jom._OPERATOR_SUBSTRACT" , true}); // minus
		listOperators.put("- 1", new Object [] {"-", 1, Operator.Associativity.RIGHT, 3, "com.jom._OPERATOR_NEG" , false}); // negate
		listOperators.put("* 2", new Object [] {"*", 2, Operator.Associativity.LEFT, 2, "com.jom._OPERATOR_MULTIPLY" , true}); // e2eMultiply
		listOperators.put(".* 2", new Object [] {".*", 2, Operator.Associativity.LEFT, 2, "com.jom._OPERATOR_E2EMULTIPLY" , true}); // e2eMultiply
		listOperators.put("/ 2", new Object [] {"/", 2, Operator.Associativity.LEFT, 2, "com.jom._OPERATOR_E2EDIVIDE" , true}); // e2eMultiply
		listOperators.put("./ 2", new Object [] {"./", 2, Operator.Associativity.LEFT, 2, "com.jom._OPERATOR_E2EDIVIDE" , true}); // e2eMultiply
		listOperators.put("^ 2", new Object [] {"^", 2, Operator.Associativity.LEFT, 3, "com.jom._OPERATOR_POW" , true}); // pow
		listOperators.put(".^ 2", new Object [] {"^", 2, Operator.Associativity.LEFT, 3, "com.jom._OPERATOR_POW" , true}); // pow
		listOperators.put("' 1", new Object [] {"'", 1, Operator.Associativity.LEFT, 4, "com.jom._OPERATOR_TRANSPOSE" , false}); // e2eMultiply

		listFunctions = new HashMap<String, Object []> ();
		listFunctions.put ("sum" , new Object [] {1 , 2 , "com.jom._FUNCTION_SUM" , false});
		listFunctions.put ("ln" , new Object [] { 1 , 1 , "com.jom._FUNCTION_LN" , false});
		listFunctions.put ("exp" , new Object [] { 1 , 1 , "com.jom._FUNCTION_EXP" , false});
		listFunctions.put ("sqrt" , new Object [] {1 , 1 , "com.jom._FUNCTION_SQRT" , false});
		listFunctions.put ("sin" , new Object [] { 1 , 1 , "com.jom._FUNCTION_SIN" , false});
		listFunctions.put ("cos" , new Object [] { 1 , 1 , "com.jom._FUNCTION_COS" , false});
		listFunctions.put ("tan" , new Object [] { 1 , 1 , "com.jom._FUNCTION_TAN" , false});
		listFunctions.put ("asin" , new Object [] { 1 , 1 , "com.jom._FUNCTION_ASIN" , false});
		listFunctions.put ("acos" , new Object [] { 1 , 1 , "com.jom._FUNCTION_ACOS" , false});
		listFunctions.put ("atan" , new Object [] { 1 , 1 , "com.jom._FUNCTION_ATAN" , false});
		listFunctions.put ("erlangB" , new Object [] {2 , 2 , "com.jom._FUNCTION_ERLANGB_TRAF" , true});
//		listFunctions.put ("matProd" , new Object [] {4 , 4 , "com.jom._FUNCTION_MATPROD" , false});
		listFunctions.put ("permute" , new Object [] {2 , 2 , "com.jom._FUNCTION_PERMUTE" , false});
		listFunctions.put ("diag" , new Object [] {1 , 1 , "com.jom._FUNCTION_DIAG" , false});
		listFunctions.put ("ones" , new Object [] {1 , 1 , "com.jom._FUNCTION_ONES" , false});
		listFunctions.put ("zeros" , new Object [] {1 , 1 , "com.jom._FUNCTION_ZEROS" , false});
		listFunctions.put ("eye" , new Object [] {1 , 2 , "com.jom._FUNCTION_EYE" , false});
}
	
	_INTERNAL_ExpressionParser(OptimizationProblem model, HashMap<String, DoubleMatrixND> inputParameters , HashMap<String, _INTERNAL_DecisionVariableArray> decisionVariables)
	{
		super(computeParametersObject (inputParameters , decisionVariables));
		this.model = model;
		this.inputParameters = inputParameters;
		this.decisionVariables = decisionVariables;
	}

	private static Parameters computeParametersObject (HashMap<String, DoubleMatrixND> inputParameters , HashMap<String, _INTERNAL_DecisionVariableArray> decisionVariables)
	{
		Parameters param = new Parameters();

		/* Define constants of the parser */
		ArrayList<Constant> listConstants = new ArrayList<Constant>();
		listConstants.add(new Constant("all"));
		for (int cont = 0; cont < listConstants.size(); cont++)
			param.add(listConstants.get(cont));

		/* Input parameters are functions (to define ranges) */
		Iterator it = inputParameters.entrySet().iterator();
		while (it.hasNext())
		{
			Map.Entry pairs = (Map.Entry) it.next();
			String name = (String) pairs.getKey();
			param.add(new Function(name, 0, OptimizationProblem.MAX_NUMBER_DIMENSIONS_INPUTPARAMETER));
			param.add(new Function(name, 0, OptimizationProblem.MAX_NUMBER_DIMENSIONS_INPUTPARAMETER));
		}

		/* Decision variables are both constants (full array) and functions (to define ranges) */
		it = decisionVariables.entrySet().iterator();
		while (it.hasNext())
		{
			Map.Entry pairs = (Map.Entry) it.next();
			String name = (String) pairs.getKey();
			_INTERNAL_DecisionVariableArray decisionVariable = (_INTERNAL_DecisionVariableArray) pairs.getValue();
			param.add (new Function(name, 0, (int) decisionVariable.getVarIds().getSize().size()));
		}

		/* Add the operators */
		for (Entry<String, Object []> entry : listOperators.entrySet())
		{
			String key = entry.getKey();
			Object [] var = entry.getValue();
		// symbol, num operans, associativity, priority, class to load, expand scalar->array
			String symbol = (String) var [0]; 
			int operandCount = (Integer) var [1];
			Operator.Associativity associativity = (Operator.Associativity) var [2];
			int precedence = (Integer) var [3];
			String className = (String) var [4];
			boolean expandScalarToArray = (Boolean) var [5];

			param.add(new Operator(symbol, operandCount, associativity, precedence));
		}

		/* *Add the functions of the parser */
		for (Entry<String,Object []> entry : listFunctions.entrySet())
		{
			String name = entry.getKey();
			Object [] var = entry.getValue();
			int minNumParameters = (Integer) var [0];
			int maxNumParameters = (Integer) var [1];
			String functionClassName = (String) var [2];
			boolean expandScalarToArray = (Boolean) var [3];
			param.add(new Function(name, minNumParameters , maxNumParameters));
		}

		param.addExpressionBracket(bracketPairs);
		param.addFunctionBracket(bracketPairs);
		param.setFunctionArgumentSeparator(',');

		return param;
	}
		
	@Override
	public Expression evaluate(String expression)
	{
		String newExp = addEmptyParenthesisWhereRequired(expression);
		newExp = addParenthesisForExplicitMatrices(newExp);
		return super.evaluate(newExp);
	}

	@Override
	public Expression evaluate(String expression, Object evaluationContext)
	{
		String newExp = addEmptyParenthesisWhereRequired(expression);
		newExp = addParenthesisForExplicitMatrices(newExp);
		return super.evaluate(newExp, evaluationContext);
	}

	//

	@Override
	protected Expression evaluate(Constant constant, Object evaluationContext)
	{
		if (constant.getName().equals("all")) return new _FUNCTION_LINEAREXPRESSION(model, -1.0);

		throw new JOMException ("Syntax error in the JOM expression. Unknown identifier '" + constant.getName() + "'");
	}

	@Override
	protected Expression evaluate(Function function, Iterator arguments, Object evaluationContext)
	{
		/* Compute the parameters list  */
		final String name = function.getName();
		ArrayList<Expression> argList = new ArrayList<Expression> ();
		while (arguments.hasNext()) argList.add((Expression) arguments.next());
		argList.trimToSize();
		
		final DoubleMatrixND inputParameterArray = model.getInputParameter(name);
		if (inputParameterArray != null)
		{
			/* No arguments: the complete constant */
			if (argList.size() == 0)
			{
				Expression aux = new _FUNCTION_LINEAREXPRESSION(model, inputParameterArray);
				return aux;
			}

			/* One argument: a range of indexes. The output is an array of the same size as the indexes array, with the given coordinates */
			if (argList.size() == 1)
			{
				Expression rangeIndexes = argList.get(0);
				if (!rangeIndexes.isConstant()) throw new JOMException ("Syntax error in the JOM expression. The indexes ranges must be constants");
				IntMatrix1D indexes = rangeIndexes.evaluateConstant().cast2Int().elements();
				/* The "all" clausule */
				int [] newSize = rangeIndexes.getSize(); 
				if ((indexes.size() == 1) && (indexes.get(0) == -1))
				{
					newSize = new int [] {1 , inputParameterArray.getNumElements()};
					indexes = new IntMatrixND (newSize).ascending(0, 1).elements();
				}
				
				Expression aux = new _FUNCTION_LINEAREXPRESSION(model, inputParameterArray.viewSelectionByIndexes(newSize , indexes));
				return aux;
			}

			if (argList.size() != inputParameterArray.getNumDim()) throw new JOMException ("Syntax error in the JOM expression. Wrong number of dimensions");

			/* Arguments: as many as dimensions, they are the constant ranges */
			String stringExpression = "(" + name + "(";
			IntMatrix1D[] coordinatesPerDim = new IntMatrix1D[(int) inputParameterArray.getSize().size()];
			for (int contDim = 0; contDim < inputParameterArray.getNumDim () ; contDim++)
			{
				Expression rangeThisDim = argList.get(contDim);
				if (!rangeThisDim.isConstant()) throw new JOMException ("Syntax error in the JOM expression. The subarray ranges must be constants");
				if (!rangeThisDim.isRowOrColumn2DExpression ()) throw new JOMException ("Syntax error in the JOM expression. The subarray ranges must be 1D expressions or scalars");
				if (contDim != 0) stringExpression = stringExpression + ",";
				if (rangeThisDim.isScalar() && (rangeThisDim.evaluateConstant().get(0) == -1))
				{
					coordinatesPerDim[contDim] = null; // "all" keyword
					stringExpression = stringExpression + "all";
				} else
				{
					coordinatesPerDim[contDim] = rangeThisDim.evaluateConstant().cast2Int().elements();
					stringExpression = stringExpression; // + rangeThisDim.getExpressionString();
				}
			}
			stringExpression = stringExpression + ")";

			DoubleMatrixND res = inputParameterArray.viewSelection(coordinatesPerDim);
			Expression aux = new _FUNCTION_LINEAREXPRESSION(this.model, res);
			return aux;

		}

		_INTERNAL_DecisionVariableArray dv = model.getDecisionVariable(name);
		if (dv != null)
		{
			/* No arguments: the complete decision variable */
			if (argList.size() == 0)
			{
				Expression aux = new _FUNCTION_LINEAREXPRESSION(model, dv);
				return aux;
			}

			/* One argument: a range of indexes. The output is an array of the same size as the indexes array, with the given coordinates of the DV */
			if (argList.size() == 1)
			{
				_FUNCTION_LINEAREXPRESSION rangeIndexes = (_FUNCTION_LINEAREXPRESSION) argList.get(0);
				if (!rangeIndexes.isConstant()) throw new JOMException ("Syntax error in the JOM expression. The indexes ranges must be constants");
				IntMatrix1D indexes = rangeIndexes.evaluateConstant().cast2Int().elements();
				IntMatrix1D newSize = IntFactory1D.dense.make(rangeIndexes.getSize());
				/* "all" clausule => all the decision variables */
				if ((indexes.size() == 1) && (indexes.get(0) == -1))
				{
					newSize = IntFactory1D.dense.make(new int [] {1 , dv.getNumDecVar()});
					indexes = new IntMatrixND (newSize.toArray()).ascending(0, 1).elements();
				}
				String stringExpression = "(" + name + "(" + (((indexes.size() == 1) && (indexes.get(0) == -1))? "all" : "expression string") + ") )";
				IntMatrixND originalVarIds = dv.getVarIds();
				IntMatrixND res = originalVarIds.viewSelectionByIndexes (newSize , indexes);
				Expression aux = new _FUNCTION_LINEAREXPRESSION(this.model, res);
				return aux;
			}

			if (argList.size() != dv.getNumDim()) throw new JOMException ("Syntax error in the JOM expression. Wrong number of dimensions");

			/* Arguments: as many as dimensions, they are the constant ranges */
			String stringExpression = "(" + name + "(";
			IntMatrix1D[] coordinatesPerDim = new IntMatrix1D[(int) dv.getSize().size()];
			for (int contDim = 0; contDim < dv.getNumDim(); contDim++)
			{
				_FUNCTION_LINEAREXPRESSION rangeThisDim = (_FUNCTION_LINEAREXPRESSION) argList.get(contDim);
				if (!rangeThisDim.isConstant()) throw new JOMException ("Syntax error in the JOM expression. The subarray ranges must be constants");
				if (!rangeThisDim.isRowOrColumn2DExpression ()) throw new JOMException ("Syntax error in the JOM expression. The subarray ranges must be 1D expressions or scalars");
				if (contDim != 0) stringExpression = stringExpression + ",";
				if (rangeThisDim.isScalar() && (rangeThisDim.getAffineExpression().getCellConstantCoef(0) == -1))
				{
					coordinatesPerDim[contDim] = null;
					stringExpression = stringExpression + "all";
				} else
				{
					coordinatesPerDim[contDim] = rangeThisDim.evaluateConstant().cast2Int().elements();
					stringExpression = stringExpression; // + rangeThisDim.getExpressionString();
				}
			}
			stringExpression = stringExpression + ")";

			IntMatrixND originalVarIds = dv.getVarIds();
			IntMatrixND res = originalVarIds.viewSelection(coordinatesPerDim);
			Expression aux = new _FUNCTION_LINEAREXPRESSION(this.model, res);
			return aux;
		}


		
		/* Here if it is a JOM function (sum, ln, ...) */
		Object [] functionInfo = listFunctions.get(name);
		if (functionInfo == null) throw new JOMException ("Syntax error in the JOM expression. Unkonwn function: '" + name + "'");
		int minNumParameters = (Integer) functionInfo [0];
		int maxNumParameters = (Integer) functionInfo [1];
		int numParameters = argList.size(); if ((numParameters > maxNumParameters) || (numParameters < minNumParameters)) throw new JOMException ("Syntax error in the JOM expression. Wrong number of parameters");  
		String functionFileName = (String) functionInfo [2];
		boolean expandScalarToArray = (Boolean) functionInfo [3];

		/* Expand scalars to arrays if needed */
		if (expandScalarToArray) expandScalarArgumentsToArray (argList);
		
		/* Instantiate the class that executes the operand */
		Class functionClass;
		try { functionClass = _INTERNAL_ExpressionParser.class.getClassLoader().loadClass(functionFileName); } catch (Exception e)
		{
			e.printStackTrace();
			throw new JOMException ("Unexpected exception trying to load class of '" + functionFileName + "'");
		}
		Expression res;
		try
		{
			if (numParameters == 1)
				res = (Expression) functionClass.getDeclaredConstructor(new Class[] { OptimizationProblem.class, Expression.class }).newInstance(this.model, argList.get(0));
			else if (numParameters == 2)
				res = (Expression) functionClass.getDeclaredConstructor(new Class[] { OptimizationProblem.class, Expression.class, Expression.class }).newInstance(this.model, argList.get(0), argList.get(1));
			else if (numParameters == 4)
				res = (Expression) functionClass.getDeclaredConstructor(new Class[] { OptimizationProblem.class, Expression.class, Expression.class , Expression.class , Expression.class }).newInstance(this.model, argList.get(0), argList.get(1) , argList.get(2) , argList.get(3)); // PABLO
			else
				throw new JOMException ("Syntax error in the JOM expression. Operators with more than two operands are not supported");
			return res;
		} catch (Exception e)
		{
			if (e instanceof RuntimeException) throw (RuntimeException) e;
			e.printStackTrace();
			throw new JOMException ("Unexpected exception trying to load class of '" + functionFileName + "'");
		}
		
	}

	@Override
	protected Expression evaluate(Operator operator, Iterator<Expression> operands, Object evaluationContext)
	{
		ArrayList<Expression> arguments = new ArrayList<Expression> ();
		while (operands.hasNext()) arguments.add((Expression) operands.next());
		if (arguments.size () != operator.getOperandCount()) throw new JOMException ("Unexpected error");
		
		/* Retrieve the operator info */
		Object [] opInfo = this.listOperators.get(operator.getSymbol() + " " + operator.getOperandCount());
		if (opInfo == null) throw new JOMException ("Syntax error in the JOM expression. Unkonwn operator '" + operator.getSymbol() + "'");
		String operatorSymbol = (String) opInfo [0]; if (!operatorSymbol.equals(operator.getSymbol())) throw new JOMException ("Unexpected error");
		int operandCount = (Integer) opInfo [1]; if (operandCount != operator.getOperandCount()) throw new JOMException ("Unexpected error");
		String functionClassName = (String) opInfo [4];
		boolean expandScalarToArray = (Boolean) opInfo [5];

		/* For operators + , - , * , / , .* , ./ , .^ => scalar operators should be expanded to array operators */
		/* In the case of * : In this case the operator is changed to .*  */
		boolean scalarToExpressionChange = false;
		if (expandScalarToArray) 
			scalarToExpressionChange = expandScalarArgumentsToArray (arguments);

		if (scalarToExpressionChange && operator.getSymbol ().equals ("*") ) { operatorSymbol = ".*"; functionClassName = "com.jom._OPERATOR_E2EMULTIPLY"; }

		/* Instantiate the class that executes the operand */
		Class operatorGenerationClasss;
		try { operatorGenerationClasss = _INTERNAL_ExpressionParser.class.getClassLoader().loadClass(functionClassName); } catch (Exception e)
		{
			e.printStackTrace();
			throw new JOMException ("Exception trying to load class of '" + functionClassName + "'");
		}
		Expression res;
		try
		{
			if (operandCount == 1)
				res = (Expression) operatorGenerationClasss.getDeclaredConstructor(new Class[] { OptimizationProblem.class, Expression.class }).newInstance(this.model, arguments.get(0));
			else if (operandCount == 2)
				res = (Expression) operatorGenerationClasss.getDeclaredConstructor(new Class[] { OptimizationProblem.class, Expression.class, Expression.class }).newInstance(this.model, arguments.get(0), arguments.get(1));
			else
				throw new JOMException ("Syntax error in the JOM expression. Operators with more than two operands are not supported");

			return res;
		} catch (Exception e)
		{
			if (e instanceof RuntimeException) throw (RuntimeException) e;
			e.printStackTrace();
			throw new JOMException ("Exception executing the class '" + functionClassName + "'");
		}

	}

	@Override
	protected Expression toValue(String literal, Object evaluationContext)
	{
		String trimmedLiteral = literal.trim();

		double number;
		// if it is a number, return a number
		number = Double.parseDouble(trimmedLiteral);
		_FUNCTION_LINEAREXPRESSION res = new _FUNCTION_LINEAREXPRESSION(model, number);
		return res;
	}

	private String addEmptyParenthesisWhereRequired(String exp)
	{
		final Iterator<String> tokens = this.tokenizer.tokenize(exp);
		String leftBracketsAllowed = bracketPairs.getOpen();

		String modifiedExpression = "";

		while (tokens.hasNext())
		{
			String token = tokens.next();
			modifiedExpression = modifiedExpression + token;

			if ((model.getDecisionVariable(token.trim()) != null) || (model.getInputParameter(token.trim()) != null))
			{
				if (!tokens.hasNext())
					modifiedExpression = modifiedExpression + "()";
				else
				{
					String nextToken = tokens.next();

					if (leftBracketsAllowed.contains(nextToken.trim().substring(0, 1)))
						modifiedExpression = modifiedExpression + nextToken;
					else
						modifiedExpression = modifiedExpression + " () " + nextToken;
				}
			}

		}

		return modifiedExpression;
	}
	
	private String addParenthesisForExplicitMatrices (String exp)
	{
		String s = exp.replaceAll(Pattern.quote("["), "([( ( ("); // open parenthesis for: all before ,[, all after [, the row, the cell
		s = s.replaceAll(Pattern.quote("]"), ") ) ) )"); // close parenthesis for: all before ,[, all after [, the row, the cell
		s = s.replaceAll(Pattern.quote(";"), ");("); // close & open parenthesis for the cell
		s = s.replaceAll(Pattern.quote(");();("), "));;(("); // close & open parenthesis for the cell and the row
		return s;
	}


	private boolean expandScalarArgumentsToArray (ArrayList<Expression> arguments)
	{
		if (arguments.size() != 2) throw new JOMException ("Unexpected error");
		if ((arguments.get(0) == null) || (arguments.get(1) == null)) throw new JOMException ("Unexpected error");
		
		if (arguments.get(0).isScalar() && !arguments.get(1).isScalar())
		{
			arguments.set(0, new _FUNCTION_REPEATEXPRESSION (arguments.get(0).getModel() , arguments.get(0), arguments.get(1).getSize()) );
			return true;
		}
		if (arguments.get(1).isScalar() && !arguments.get(0).isScalar())
		{
			arguments.set(1, new _FUNCTION_REPEATEXPRESSION (arguments.get(1).getModel() , arguments.get(1), arguments.get(0).getSize()) );
			return true;
		}
		if (arguments.get(1).isScalar() && arguments.get(0).isScalar())
			return true;
		return false;
	}	
}
