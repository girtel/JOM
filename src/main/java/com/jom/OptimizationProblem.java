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

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;

import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tdouble.DoubleMatrix2D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.jet.math.tdouble.DoubleFunctions;

/** This class contains the methods for handling optimization problems, defining their input parameters (if any), decision variables,
 * objetive function and constraints, choosing and calling a solver to obtain a numerical solution, and retrieving that solution.
 * @author Pablo Pavon Mariño
 * @see <a href="http://www.net2plan.com/jom">http://www.net2plan.com/jom</a>
 */

public class OptimizationProblem
{
	static { System.loadLibrary("JVmipcl"); }
	
	final static int MAX_NUMBER_DIMENSIONS_INPUTPARAMETER = 10;
	private Map<String, _INTERNAL_ConstraintArray> constraints;
	private Map<String, _INTERNAL_DecisionVariableArray> decisionVariables;
	private _INTERNAL_ExpressionParser evaluator;
	private Map<String, DoubleMatrixND> initialSolution;
	private Map<String, DoubleMatrixND> inputParameters;
	private int lpNumScalarConstraints;
	private int nlpNumScalarConstraints;
	private int numScalarDecisionVariables;
	private Expression objectiveFunction;
	private _INTERNAL_SolverIO solverIO;
	private TimeInfo timeInfo;
	private boolean toMinimize;

	public enum JOMSolver
    {
        glpk, ipopt, cplex, xpress
    }

	/** Returns the list of the names of the installed solvers, in the form that is accepted as the solverName inpt parameter of solve method
	 * @return see above
	 */
	public static List<String> getAcceptedSolverNames ()
    {
        final List<String> acceptedNames = new ArrayList<>();
        for (JOMSolver solver : JOMSolver.values())
            acceptedNames.add(solver.name());
        return acceptedNames;
    }

	/** Creates an optimization problem object
	 *
	 */
	public OptimizationProblem()
	{
		this.inputParameters = new LinkedHashMap<String, DoubleMatrixND>();
		this.decisionVariables = new LinkedHashMap<String, _INTERNAL_DecisionVariableArray>();
		this.initialSolution = new LinkedHashMap<String, DoubleMatrixND>();
		this.constraints = new LinkedHashMap<String, _INTERNAL_ConstraintArray>();
		this.evaluator = null;
		this.solverIO = null;
		this.timeInfo = new TimeInfo();

		this.toMinimize = false;
		this.objectiveFunction = null;
		this.numScalarDecisionVariables = 0;

		this.lpNumScalarConstraints = 0;
		this.nlpNumScalarConstraints = 0;

		/* So that the JVM does not execute its JNA library (if any), but the one of jna.jar */
		System.setProperty("jna.nosys", "true");
	}

	/** Returns an informative string about the JOM library
	 * @return the string
	 */
	public String about(){ return "Java Optimization Modeler v. 0.1.13\nAuthor: Pablo Pavon Mari�o. \nDate: September 2015. \nJOM is open-source, "
			+ "licensed by GNU Lesser General Public License v3.0\n"; }

	
	/** Removes from the model the constraints indicated by the given identifier
	 * @param identifier the identifier, as indicated in the addConstraint method
	 * @return true if the constraints existed, and thus were removed, false otherwise
	 */
	public boolean removeConstraint (String identifier)
	{
		final _INTERNAL_ConstraintArray removedElement = this.constraints.remove(identifier);
		return removedElement != null;
	}

	/** Indicates if a constraint set with the given identifier exists in the model
	 * @param identifier the identifier of the constraint to check
	 * @return see above
	 */
	public boolean hasConstraint (String identifier)
	{
		return this.constraints.containsKey(identifier);
	}

	/** Indicates if a decision variable with the given identifier exists in the model
	 * @param identifier the identifier of the constraint to check
	 * @return see above
	 */
	public boolean hasDecisionVariable (String identifier)
	{
		return this.decisionVariables.containsKey(identifier);
	}

	/** Adds (an array of) constraints to the optimization problem. No identifier is assigned, so after the problem is solved, it is not possible to
	 *  access the multipliers and other constraint-related info.
	 * Constraints have the form: expression-left-hand-side connector expression-right-hand-side.
	 * Expressions are strings that can be evaluated to arrayed expressions of the defined input parameters and decision variables. Connectors must
	 * be of the type less-equal {@literal (<= or =<)}, greater-equal {@literal (>= or =>)} or equal {@literal (==)}.
	 * @param expression The string expression of the constraint.
	 * @return The Expression object holding the expression (lhs - rhs), where lhs and rhs are the expression in the left-hand-side and
	 * right-hand-side of the constraint
	 */
	public Expression addConstraint(String expression)
	{
		String identifier = new Double(Math.random()).toString();
		while (this.constraints.containsKey(identifier))
			identifier = new Double(Math.random()).toString();
		return this.addConstraint(expression, identifier);
	}

	/** Adds (an array of) constraints to the optimization problem.
	 * Constraints have the form: expression-left-hand-side connector expression-right-hand-side.
	 * Expressions are strings that can be evaluated to arrayed expressions of the defined input parameters and decision variables. Connectors must
	 * be of the type less-equal {@literal (<= or =<)}, greater-equal {@literal (>= or =>)} or equal {@literal (==)}.
	 * @param expression The string expression of the constraint.
	 * @param identifier A unique string given to identify this constraint
	 * @return The Expression object holding the expression (lhs - rhs), where lhs and rhs are the expression in the left-hand-side and
	 * right-hand-side of the constraint
	 */
	public Expression addConstraint(String expression, String identifier)
	{
		long initTime = System.nanoTime();

		if (evaluator == null) this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);
		if (this.constraints.containsKey(identifier)) throw new JOMException("This constraint identifier already exists");

		final _INTERNAL_ConstraintArray cs = new _INTERNAL_ConstraintArray(evaluator, identifier, expression, this.lpNumScalarConstraints + this
				.nlpNumScalarConstraints); // PABLO: Changed this

		this.constraints.put(identifier, cs);
		if (cs.isLinear())
			this.lpNumScalarConstraints += cs.getNumScalarConstraints();
		else
			this.nlpNumScalarConstraints += cs.getNumScalarConstraints();

		checkProblemSizeLimitation(this.numScalarDecisionVariables, this.lpNumScalarConstraints, this.nlpNumScalarConstraints);
		timeInfo.timeSettingConstraints.add(System.nanoTime() - initTime);
		return cs.getLhsMinusRhs();
	}

	/** Same as addDecisionVariable(name, isInteger, size, null, null);
	 * @param name Name of the array of decision variables, as it will used in the expressions.
	 * @param isInteger True if the decision variables in the array are all constrained to be integer
	 * @param size One coordinate per dimension of the array, each coordinate is the size of the array in the respective dimension.
	 */
	public void addDecisionVariable(String name, boolean isInteger, int[] size)
	{
		addDecisionVariable(name, isInteger, size, (DoubleMatrixND) null, (DoubleMatrixND) null);
	}

	/** addDecisionVariable(String name, boolean isInteger, int[] size, DoubleMatrixND x_l, DoubleMatrixND x_u), but now the lower (upper) bounds
	 * for all the decision variables in the array are equal to x_l (x_u).
	 * @param name Name of the array of decision variables, as it will used in the expressions.
	 * @param isInteger True if the decision variables in the array are all constrained to be integer
	 * @param size One coordinate per dimension of the array, each coordinate is the size of the array in the respective dimension.
	 * @param x_l The lower bounds of the variables are created as new DoubleMatrixND (size,x_l)
	 * @param x_u The upper bounds of the variables are created as new DoubleMatrixND (size,x_u)
	 */
	public void addDecisionVariable(String name, boolean isInteger, int[] size, double x_l, double x_u)
	{
		addDecisionVariable(name, isInteger, size, new DoubleMatrixND(size, x_l, "sparse"), new DoubleMatrixND(size, x_u, "sparse"));
	}

	/** Adds an array of decision variables to the optimization problem
	 * @param name Name of the array of decision variables, as it will used in the expressions.
	 * @param isInteger True if the decision variables in the array are all constrained to be integer
	 * @param size One coordinate per dimension of the array, each coordinate is the size of the array in the respective dimension.
	 * @param x_l The lower bounds of the variables are created as new DoubleMatrixND (size,x_l)
	 * @param x_u The upper bounds of the variables are created as new DoubleMatrixND (size,x_u)
	 */
	public void addDecisionVariable(String name, boolean isInteger, int[] size, double[] x_l, double[] x_u)
	{
		if (objectiveFunctionOrConstraintsAdded())
			throw new JOMException("Decision variables cannot be added after objective function or constraints are set");

		final int numberNewElements = IntMatrixND.prod(IntFactory1D.dense.make(size));
		final _INTERNAL_DecisionVariableArray v = new _INTERNAL_DecisionVariableArray(name, isInteger, IntFactory1D.dense.make(size), new DoubleMatrixND
				(size, x_l), new DoubleMatrixND(size, x_u), this.numScalarDecisionVariables);
		decisionVariables.put(name, v);
		this.numScalarDecisionVariables += numberNewElements;
		checkProblemSizeLimitation(this.numScalarDecisionVariables, this.lpNumScalarConstraints, this.nlpNumScalarConstraints);
	}

	/** Adds an array of decision variables to the optimization problem
	 * @param name Name of the array of decision variables, as it will used in the expressions.
	 * @param isInteger True if the decision variables in the array are all constrained to be integer
	 * @param size One coordinate per dimension of the array, each coordinate is the size of the array in the respective dimension.
	 * @param x_l The lower bounds of the variables are created as new DoubleMatrixND (size,x_l)
	 * @param x_u The upper bounds of the variables are created as new DoubleMatrixND (size,x_u)
	 */
	public void addDecisionVariable(String name, boolean isInteger, int[] size, DoubleMatrix1D x_l, DoubleMatrix1D x_u)
	{
		addDecisionVariable(name, isInteger, size, x_l.toArray(), x_u.toArray());
	}

	/** Adds an array of decision variables to the optimization problem
	 * @param name Name of the array of decision variables, as it will used in the expressions.
	 * @param isInteger True if the decision variables in the array are all constrained to be integer
	 * @param size One coordinate per dimension of the array, each coordinate is the size of the array in the respective dimension.
	 * @param x_l An array of doubles of the same size of the decision variables, with the lower bounds of the variables.
	 * @param x_u An array of doubles of the same size of the decision variables, with the upper bounds of the variables.
	 */
	public void addDecisionVariable(String name, boolean isInteger, int[] size, DoubleMatrixND x_l, DoubleMatrixND x_u)
	{
		if (objectiveFunctionOrConstraintsAdded())
			throw new JOMException("Decision variables cannot be added after objective function or constraints are set");

		final int numberNewElements = IntMatrixND.prod(IntFactory1D.dense.make(size));
		final _INTERNAL_DecisionVariableArray v = new _INTERNAL_DecisionVariableArray(name, isInteger, IntFactory1D.dense.make(size), x_l, x_u, this
				.numScalarDecisionVariables);
		decisionVariables.put(name, v);
		this.numScalarDecisionVariables += numberNewElements;
		checkProblemSizeLimitation(this.numScalarDecisionVariables, this.lpNumScalarConstraints, this.nlpNumScalarConstraints);
	}

	/** Sets the isInteger flag for the decision variables array identified by the given name
	 * @param decisionVariableName the name, as introduced with addDecisionVariable
	 * @param isInteger true if the variables are contrained to be integer
	 */
	public void setIsIntegerDecisionVariable (String decisionVariableName , boolean isInteger)
	{
		if (!decisionVariables.containsKey(decisionVariableName)) 
			throw new JOMException ("Unknown decision variables");
		decisionVariables.get(decisionVariableName).setIsInteger(isInteger);
	}

	/** Sets the lower bound limit of the indicated decision variable array
	 * @param decisionVariableName the name, as introduced with addDecisionVariable
	 * @param x_l an array with the same size as the decision variable, with the new bound
	 */
	public void setLowerBoundDecisionVariable (String decisionVariableName , DoubleMatrixND x_l)
	{
		if (!decisionVariables.containsKey(decisionVariableName)) 
			throw new JOMException ("Unknown decision variables");
		decisionVariables.get(decisionVariableName).set_x_l(x_l);
	}

   /** Sets the lower bound limit of the indicated decision variable array
     * @param decisionVariableName the name, as introduced with addDecisionVariable
     * @param x_l a constant, the new bound
     */
    public void setLowerBoundDecisionVariable (String decisionVariableName , double x_l)
    {
        if (!decisionVariables.containsKey(decisionVariableName)) 
            throw new JOMException ("Unknown decision variables");
        decisionVariables.get(decisionVariableName).set_x_l(new DoubleMatrixND(decisionVariables.get(decisionVariableName).getSize().toArray() , x_l));
    }

	/** Sets the upper bound limit of the indicated decision variable array
	 * @param decisionVariableName the name, as introduced with addDecisionVariable
	 * @param x_u an array with the same size as the decision variable, with the new bound
	 */
	public void setUpperBoundDecisionVariable (String decisionVariableName , DoubleMatrixND x_u)
	{
		if (!decisionVariables.containsKey(decisionVariableName)) 
			throw new JOMException ("Unknown decision variables");
		decisionVariables.get(decisionVariableName).set_x_u(x_u);
	}
	
	/** Returns true if the problem was attempted to be solved, and the solver declares that the problem has no feasible solutions. If the method
	 * "solve" was not invoked yet, throws an exception
	 * @return See above
	 */
	public boolean feasibleSolutionDoesNotExist()
	{
		if (!this.solverIO.problemAlreadyAttemptedTobeSolved) throw new JOMException("The problem has not been solved yet");
		else return this.solverIO.out.feasibleSolutionDoesNotExist;
	}

	/** Returns true if the problem was attempted to be solved, and the solver declares that the problem is unbounded. If the method "solve" was not
	 *  invoked yet, throws an exception
	 * @return See above
	 */
	public boolean foundUnboundedSolution()
	{
		if (!this.solverIO.problemAlreadyAttemptedTobeSolved) throw new JOMException("The problem has not been solved yet");
		else return this.solverIO.out.foundUnboundedSolution;
	}

	/** Obtains the value (an array of doubles) of the input parameter previously set
	 * @param name Name of the input parameter as it appears in the expressions
	 * @return The value of the input parameter
	 */
	public DoubleMatrixND getInputParameter(String name)
	{
		return inputParameters.get(name);
	}

	/** Returns the multipliers of the automatic constraints: \(x_l \leq varName\) added to the problem, where \(x_l\) is the array of lower bounds provided,
	 * associated to the decision variable.
	 * @param decisionVariableName The name of the decision variable
	 * @return The multipliers as an array of the same size of the decision variables
	 */
	public DoubleMatrixND getMultiplierOfLowerBoundConstraintToPrimalVariables(String decisionVariableName)
	{
		_INTERNAL_DecisionVariableArray dv = this.decisionVariables.get(decisionVariableName);
		if (dv == null) throw new JOMException("Unknown decision variable name: " + decisionVariableName);
		DoubleMatrix1D values = this.solverIO.out.multiplierOfLowerBoundConstraintToPrimalVariables.viewSelection(dv.getVarIds().elements().toArray
				());
		return new DoubleMatrixND(dv.getVarIds().getSize(), values.copy());
	}

	/** Returns the multipliers of the automatic constraints: \(varName \leq x_u\) added to the problem, where \(x_u\) is the array of upper bounds provided,
	 * associated to the decision variable.
	 * @param decisionVariableName The name of the decision variable
	 * @return The multipliers as an array of the same size of the decision variables
	 */
	public DoubleMatrixND getMultiplierOfUpperBoundConstraintToPrimalVariables(String decisionVariableName)
	{
		_INTERNAL_DecisionVariableArray dv = this.decisionVariables.get(decisionVariableName);
		if (dv == null) throw new JOMException("Unknown decision variable name: " + decisionVariableName);
		DoubleMatrix1D values = this.solverIO.out.multiplierOfUpperBoundConstraintToPrimalVariables.viewSelection(dv.getVarIds().elements().toArray
				());
		return new DoubleMatrixND(dv.getVarIds().getSize(), values.copy());
	}

	/** Returns the multipliers of the given constraint
	 * @param constraintIdentifier The unique identifier of the constraint
	 * @return The multipliers as an array of the same size of the constraints
	 */
	public DoubleMatrixND getMultipliersOfConstraint(String constraintIdentifier)
	{
		_INTERNAL_ConstraintArray cs = this.constraints.get(constraintIdentifier);
		if (cs == null) throw new JOMException("Unknown constraint set name: " + constraintIdentifier);
		return new DoubleMatrixND(cs.getSize(), this.solverIO.out.multiplierOfConstraint.viewPart(cs.getIndex_0_FirstConstraint(), cs
				.getNumScalarConstraints()));
	}

	/** Gets the number of scalar constraints in the problem that are linear (with integer variables or not). Arrayed constraints are counted as
	 * many times as its number of constituting scalar elements.
	 * @return The number of linear scalar constraints
	 */
	public int getNumLinearScalarConstraints()
	{
		return this.lpNumScalarConstraints;
	}

	/** Gets the number of scalar constraints in the problem that are non linear (with integer variables or not). Arrayed constraints are counted as
	 *  many times as its number of constituting scalar elements.
	 * @return The number of non-linear scalar constraints
	 */
	public int getNumNonLinearScalarConstraints()
	{
		return this.nlpNumScalarConstraints;
	}

	/** Gets the number of scalar decision variables. Arrayed decision variables are counted as many times as its number of constituting scalar
	 * elements.
	 * @return The number of scalar decision variables
	 */
	public int getNumScalarDecisionVariables()
	{
		return this.numScalarDecisionVariables;
	}

	/** Gets the objective function expression set in the problem
	 * @return The objective expression of the objective function
	 */
	public Expression getObjectiveFunction()
	{
		return this.objectiveFunction;
	}

	/** Returns the cost of the solution obtained by the solver
	 * @return The primal solution cost
	 */
	public double getOptimalCost()
	{
		return this.solverIO.out.primalCost;
	}

	/** Returns the best optimality bound found by the solver. This is a lower bound of the value of an optimum solution in minimization
	 * problems, and an upper bound in maximization problems. Some solvers are able to provide this information, even
	 * if a feasible solution is not found. If this information is not available, {@code -Double.MAX_VALUE} is returned in
	 * minimization problems, and {@code Double.MAX_VALUE} is returned in maximization problems.
	 * @return The optimality bound
	 */
	public double getBestOptimalityBound()
	{
		return this.solverIO.out.bestOptimalityBound;
	}

	/** Returns the primal solution obtained after solving the problem for the given (array of) decision variables. If the
	 * decision variable is constrained to be integer, it is rounded. In this case, we are assuming that non-integer values returned by the solver
	 * were caused by negligible numerical errors. Recall that, after the solver ends, JOM checks if its exit code signals an error-free execution.
	 * Thus, if the solver returns an "everything ok" code, numerical errors in integer variables should be negligible.
	 * @param decisionVariable The name of the decision variable
	 * @return An array of the same size of the decision variables, with the solution
	 */
	public DoubleMatrixND getPrimalSolution(String decisionVariable)
	{
		_INTERNAL_DecisionVariableArray dv = this.decisionVariables.get(decisionVariable);
		if (dv == null) throw new JOMException("Unknown decision variable name: " + decisionVariable);
		DoubleMatrix1D values = this.solverIO.out.primalSolution.viewSelection(dv.getVarIds().elements().toArray());
		if (dv.getIsInteger()) values.assign(DoubleFunctions.round(1));
		return new DoubleMatrixND(dv.getVarIds().getSize(), values.copy());
	}

	/** Returns the slack of the given constraints. It should be zero for equality constraints.
	 * @param constraintIdentifier Id of the constrain
	 * @return An array of the same size of the constraints, with the slack values
	 */
	public DoubleMatrixND getSlackOfConstraint(String constraintIdentifier)
	{
		_INTERNAL_ConstraintArray cs = this.constraints.get(constraintIdentifier);

		if (cs == null) throw new JOMException("Unknown constraint set name: " + constraintIdentifier);
		return new DoubleMatrixND(cs.getSize(), this.solverIO.out.primalValuePerConstraint.viewPart(cs.getIndex_0_FirstConstraint(), cs
				.getNumScalarConstraints()));
	}

	/** Returns true if the given name corresponds to an input parameter defined for the problem
	 * @param name Name of the input parameter
	 * @return True if the previous statement is true, false otherwise
	 */
	public boolean isInputParameter(String name)
	{
		return inputParameters.containsKey(name);
	}

	/** Returns true if the problem has one or more decision variables constrained to be integer, false otherwise
	 * @return See above.
	 */
	public boolean isIntegerProblem()
	{
		for (Entry<String, _INTERNAL_DecisionVariableArray> entry : decisionVariables.entrySet())
			if (entry.getValue().getIsInteger()) return true;
		return false;
	}

	/** Returns true if the problem is linear (with some/all integer variables or not). False if the objective function or one or more constraints
	 * are non-linear.
	 * @return See above
	 */
	public boolean isLinearProblem()
	{
		if (!this.objectiveFunction.isLinear()) return false;
		return (this.getNumNonLinearScalarConstraints() == 0);
	}

	/** Returns true if the current objective function has already been defined, and the problem is set as a minimization problem
	 * @return True if the previous statement is true, false otherwise
	 */
	public boolean isToMinimize()
	{
		return this.toMinimize;
	}

	/** Help function (not needed usually for solving optimization problems) to generate an Expression object from its String representation.
	 * @param expression The string expression.
	 * @return The Expression object
	 */
	public Expression parseExpression(String expression)
	{
		if (evaluator == null) this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);
		return this.evaluator.evaluate(expression, new Object());
	}

	/** Resets the timers that track the amount of time involved in each part of the computation, to write the time report with the
	 */
	public void resetTimer(){ timeInfo.reset(); }

	/** Sets the initial solution of the optimization algorithm solving the problem (can be used for setting the starting solution in IPOPT solver).
	 * @param decisionVariableName Name of the decision variable
	 * @param initialValue Initial solution. All the scalar decision variables in the array of decision variables have this same initial value.
	 */
	public void setInitialSolution(String decisionVariableName, double initialValue)
	{
		_INTERNAL_DecisionVariableArray dv = this.decisionVariables.get(decisionVariableName);
		if (dv == null) throw new JOMException("Unknown decision variable name: " + decisionVariableName);
		this.initialSolution.put(decisionVariableName, new DoubleMatrixND(dv.getSize().toArray(), initialValue, "dense"));
	}

	/** Sets the initial solution of the optimization algorithm solving the problem (to be used by some solvers like IPOPT).
	 * @param decisionVariableName Name of the decision variable
	 * @param initialValue Initial solution. Array of the same size as the decision variable
	 */
	public void setInitialSolution(String decisionVariableName, DoubleMatrixND initialValue)
	{
		_INTERNAL_DecisionVariableArray dv = this.decisionVariables.get(decisionVariableName);
		if (dv == null) throw new JOMException("Unknown decision variable name: " + decisionVariableName);
		if (!initialValue.getSize().equals(dv.getSize()))
			throw new JOMException("Size mismatch: the decision variable size is different to the provided constant size");
		this.initialSolution.put(decisionVariableName, initialValue);
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param value The scalar value of the input parameter, converted to an N-DIM array object of size (1,1)
	 */
	public void setInputParameter(String name, double value)
	{
		boolean newInputParameterName = !(inputParameters.containsKey(name));
		inputParameters.put(name, new DoubleMatrixND(value));
		if ((this.evaluator != null) && newInputParameterName)
			this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables); // create again the evaluator, since it has a
		// new parameter name to consider in the parsing
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param values The 1D array of values of the input parameter, converted to an N-DIM array object of size (1,values.length)
	 * @param rowColumnVector "row" for a row vector, "column" for a column vector
	 */
	public void setInputParameter(String name, double[] values, String rowColumnVector)
	{
		boolean newInputParameterName = !(inputParameters.containsKey(name));
		if (rowColumnVector.toLowerCase().equals("row"))
			inputParameters.put(name, new DoubleMatrixND(new int[]{1, values.length}, values));
		else if (rowColumnVector.toLowerCase().equals("column"))
			inputParameters.put(name, new DoubleMatrixND(new int[]{values.length, 1}, values));
		else
			throw new JOMException("Indicate \"row\" or \"column\" in the rowColumnVector parameter");
		if ((this.evaluator != null) && newInputParameterName)
			this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);  // create again the evaluator, since it has
		// a new parameter name to consider in the parsing
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param values The 1D array of values of the input parameter, converted to an N-DIM array object of size (1,values.length)
	 * @param rowColumnVector "row" for a row vector, "column" for a column vector
	 */
	public void setInputParameter(String name, DoubleMatrix1D values, String rowColumnVector)
	{
		setInputParameter(name, values.toArray(), rowColumnVector);
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param values A collection with values of type Number (e.g. Double, Integer, Long...) of the input parameter
	 * @param rowColumnVector "row" for a row vector, "column" for a column vector
	 */
	public void setInputParameter(String name, Collection<? extends Number> values, String rowColumnVector)
	{
		double[] doubleVals = new double[values.size()];
		int counter = 0;
		for (Number val : values) doubleVals[counter++] = val.doubleValue();
		setInputParameter(name, doubleVals, rowColumnVector);
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param values The 2D array of values of the input parameter, converted to an N-DIM array object of the same size
	 */
	public void setInputParameter(String name, double[][] values)
	{
		boolean newInputParameterName = !(inputParameters.containsKey(name));
		inputParameters.put(name, new DoubleMatrixND(values));
		if ((this.evaluator != null) && newInputParameterName)
			this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);  // create again the evaluator, since it has
		// a new parameter name to consider in the parsing
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param array The array of values of the input parameter
	 */
	public void setInputParameter(String name, DoubleMatrixND array)
	{
		if (array.getNumDim() > MAX_NUMBER_DIMENSIONS_INPUTPARAMETER)
			throw new JOMException("The maximum number of dimensions for an input parameter array is " + MAX_NUMBER_DIMENSIONS_INPUTPARAMETER + ". "
					+ "You can increase this number setting the MAX_NUMBER_DIMENSIONS_INPUTPARAMETER option");
		boolean newInputParameterName = !(inputParameters.containsKey(name));
		inputParameters.put(name, array);
		if ((this.evaluator != null) && newInputParameterName)
			this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);  // create again the evaluator, since it has
		// a new parameter name to consider in the parsing
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param array The array of values of the input parameter
	 */
	public void setInputParameter(String name, DoubleMatrix2D array)
	{
		setInputParameter(name, new DoubleMatrixND(array));
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param value The scalar value of the input parameter, casted to double, and converted to an N-DIM array object of size (1,1)
	 */
	public void setInputParameter(String name, int value)
	{
		boolean newInputParameterName = !(inputParameters.containsKey(name));
		inputParameters.put(name, new DoubleMatrixND(value));
		if ((this.evaluator != null) && newInputParameterName)
			this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables); // create again the evaluator, since it has a
		// new parameter name to consider in the parsing
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. If an input parameter with the same name was previously set, the old
	 *  values are lost
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param values The 1D array of values of the input parameter, casted to double, and converted to an N-DIM array object of size (1,values
	 *                  .length)
	 * @param rowColumnVector "row" for a row vector, "column" for a column vector
	 */
	public void setInputParameter(String name, int[] values, String rowColumnVector)
	{
		boolean newInputParameterName = !(inputParameters.containsKey(name));
		if (rowColumnVector.toLowerCase().equals("row"))
			inputParameters.put(name, new DoubleMatrixND(new int[]{1, values.length}, values));
		else if (rowColumnVector.toLowerCase().equals("column"))
			inputParameters.put(name, new DoubleMatrixND(new int[]{values.length, 1}, values));
		else
			throw new JOMException("Indicate \"row\" or \"column\" in the rowColumnVector parameter");
		if ((this.evaluator != null) && newInputParameterName)
			this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);  // create again the evaluator, since it has
		// a new parameter name to consider in the parsing
	}

	/** Sets the value of the (arrayed) input parameter identified by its name. The String should be evaluated to a constant Expression.
	 * @param name Name of the input parameter as it appears in the expressions.
	 * @param expression array The expression is evaluated using a common parseExpression. If it is not constant, an Exception is raised.
	 */
	public void setInputParameter(String name, String expression)
	{
		boolean newInputParameterName = !(inputParameters.containsKey(name));
		inputParameters.put(name, this.parseExpression(expression).evaluateConstant());
		if ((this.evaluator != null) && newInputParameterName)
			this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);  // create again the evaluator, since it has
		// a new parameter name to consider in the parsing
	}

	/** Sets the objective function of the optimization problem, and its direction (maximice or minimice). If an objective function and direction
	 * was previously set, they are lost.
	 * @param minMax "Minimize" if the target is to minimize this expression, "Maximize" if the target is to maximize this expression
	 * @param expression The string expression of the objective function, that should successfully parsed as an scalar expression of the decision
	 *                      variables and input parameters of the problem
	 * @return The Expression object associated to the objective function
	 */
	public Expression setObjectiveFunction(String minMax, String expression)
	{
		long initTime = System.nanoTime();

		if (minMax.toLowerCase() == "minimize")
			toMinimize = true;
		else if (minMax.toLowerCase() == "maximize")
			toMinimize = false;
		else
			throw new JOMException("Invalid optimization target. Please choose \"minimize\" or \"maximize\"");

		if (evaluator == null) this.evaluator = new _INTERNAL_ExpressionParser(this, inputParameters, decisionVariables);

		Expression exp = this.evaluator.evaluate(expression, new Object());
		if (exp.getNumScalarExpressions() != 1) throw new JOMException("The objective function must be an scalar function");
		this.objectiveFunction = exp;

		timeInfo.accumTimeSettingObj += System.nanoTime() - initTime;
		checkProblemSizeLimitation(this.numScalarDecisionVariables, this.lpNumScalarConstraints, this.nlpNumScalarConstraints);
		return exp;
	}

	/** Returns true if the problem was solved, and the solution obtained is feasible. If the method "solve" was not invoked yet, throws an exception
	 * @return See above
	 */
	public boolean solutionIsFeasible()
	{
		if (!this.solverIO.problemAlreadyAttemptedTobeSolved) throw new JOMException("The problem has not been solved yet");
		else return this.solverIO.out.solutionIsFeasible;
	}

	/** Returns true if the problem was solved, and the solution obtained is optimal according to the solver. If the method "solve" was not invoked
	 * yet, throws an exception
	 * @return See above
	 */
	public boolean solutionIsOptimal()
	{
		if (!this.solverIO.problemAlreadyAttemptedTobeSolved) throw new JOMException("The problem has not been solved yet");
		else return this.solverIO.out.solutionIsOptimal;
	}

	/** Calls the indicated solver to solve the optimization problem in its current state (e.g. current objective function and input parameters
	 * set).
	 * @param solverName The name of the solver to be used. Current supported solvers are GLPK (solverName "glpk") and CPLEX (solverName "cplex")
	 *                      for mixed integer linear problems
	 * and IPOPT (solverName "ipopt") for non-linear non-integer problems.
	 * @param paramValuePairs Parameters to be passed to the solver to tune its operation. An even number of Object objects is
	 * to be passed. For each object pair, first Object must be a String with the name of the parameter, second a general Object with its value.
	 * If no name-value pairs are set, solver default values are used. For the "solverLibraryName" parameter, passing a value equal to an empty
	 *                           string ""
	 * is equivalent to not passing this parameter.
	 * See the documentation for further information.
	 */
	public void solve(String solverName, Object... paramValuePairs)
	{
		long initTime = System.nanoTime();

		this.solverIO = new _INTERNAL_SolverIO(this, this.objectiveFunction, this.toMinimize, this.decisionVariables, this.constraints, this
				.initialSolution);

		final Map<String, Object> params = createParametersMapInitilizingDefaults(solverName, paramValuePairs, this.solverIO);

		if (solverName.equalsIgnoreCase("glpk"))
		{
			if (this.solverIO.in.lhsMinusRhsAccumulatedConstraint != null)
			{
				if (!this.solverIO.in.lhsMinusRhsAccumulatedConstraint.isLinear()) throw new JOMException("The problem is not linear");
				if ((this.nlpNumScalarConstraints != 0) && this.solverIO.in.lhsMinusRhsAccumulatedConstraint.isLinear())
					throw new JOMException("Unexpected error");
			}
			if (!this.solverIO.in.objectiveFunction.isLinear()) throw new JOMException("The problem is not linear");
			_SOLVERWRAPPER_GLPK solver = new _SOLVERWRAPPER_GLPK(solverIO, params);
			int errorCode = solver.solve();
			solver.glpk_errorMessage_solve(solverIO, errorCode);
		} else if (solverName.equalsIgnoreCase("cplex"))
		{
			if (this.solverIO.in.lhsMinusRhsAccumulatedConstraint != null)
			{
				if (!this.solverIO.in.lhsMinusRhsAccumulatedConstraint.isLinear()) throw new JOMException("The problem is not linear");
				if ((this.nlpNumScalarConstraints != 0) && this.solverIO.in.lhsMinusRhsAccumulatedConstraint.isLinear())
					throw new JOMException("Unexpected error");
			}
			if (!this.solverIO.in.objectiveFunction.isLinear()) throw new JOMException("The problem is not linear");
			_SOLVERWRAPPER_CPLEX solver = new _SOLVERWRAPPER_CPLEX(solverIO, params);
			int errorCode = solver.solve();
		} else if (solverName.equalsIgnoreCase("ipopt"))
		{
			// int errorCode = ExampleJNA_IPOPT.solveIPOPTExample ();
			_SOLVERWRAPPER_IPOPT solver = new _SOLVERWRAPPER_IPOPT(solverIO, params);
			// Debug.pln(solver.printOptions ());
			solver.solve();
		} else if (solverName.equalsIgnoreCase("xpress"))
		{
			_SOLVERWRAPPER_XPRESS solver = new _SOLVERWRAPPER_XPRESS(solverIO, params);
			int errorCode = solver.solve();
		} else if (solverName.equalsIgnoreCase("mipcl"))
		{
			_SOLVERWRAPPER_MIPCL solver = new _SOLVERWRAPPER_MIPCL(solverIO, params);
			int errorCode = solver.solve();
		} else
			throw new JOMException("Unknown solver type");

		if (this.solutionIsFeasible()) solverIO.checkFeasibility();
		if (this.solutionIsOptimal()) solverIO.checkOptimality();

		timeInfo.accumTimeSolvingTheProblem += System.nanoTime() - initTime;
	}

	/** Return a printed report with the time spended in the different parts of the problem: defining the variables, objective, constraints and
	 * solving.
	 * @return a formatted string with the report
	 */
	public String timeReport(){ return this.timeInfo.toString(); }

	void checkProblemSizeLimitation(int numVarIds, int lpConstraints, int nlpConstraints)
	{
		if (numVarIds < 0) throw new JOMException("The number of decision variables must be a number between 0 and " + Integer.MAX_VALUE);
		if (lpConstraints < 0) throw new JOMException("The number of linear constraints must be a number between 0 and " + Integer.MAX_VALUE);
		if ((nlpConstraints > 0) && (numVarIds * (lpConstraints + nlpConstraints) >= Integer.MAX_VALUE))
			throw new JOMException("For problems with non linear constraints, the number of decision variables multiplied by the number of "
					+ "constraints is limited to " + Integer.MAX_VALUE);
	}

	_INTERNAL_ConstraintArray getConstraintSetInfo(int constraintIndex)
	{
		/* Decision variables are both constants (full array) and functions (to define ranges) */
		Iterator it = constraints.entrySet().iterator();
		while (it.hasNext())
		{
			Map.Entry pairs = (Map.Entry) it.next();
			String name = (String) pairs.getKey();
			_INTERNAL_ConstraintArray cs = (_INTERNAL_ConstraintArray) pairs.getValue();
			if ((cs.getIndex_0_FirstConstraint() <= constraintIndex) && (cs.getIndex_0_LastConstraint() >= constraintIndex)) return cs;
		}
		return null;
	}	/** Returns a formatted string with some information of the optimization problem
	 * @return The string
	 */
	@Override
    public String toString()
	{
		String s = "--------------------- PRINT MODEL ---\n";
		s = s + "----------- DECISION VARIABLES: \n";
		if (decisionVariables == null) s = s + "NO DECISION VARIABLES DEFINED";
		else
		{
			for (Entry<String, _INTERNAL_DecisionVariableArray> entry : decisionVariables.entrySet())
			{
				String key = entry.getKey();
				_INTERNAL_DecisionVariableArray var = entry.getValue();
				s = s + var + "\n";
			}
		}
		s = s + "----------- OBJECTIVE FUNCTION: \n";
		if (this.objectiveFunction == null)
			s += "NULL \n";
		else
			s += this.objectiveFunction.toString() + "\n";

		if (solverIO == null) return s;

		s = s + "----------- CONSTRAINTS: number = " + this.solverIO.in.numConstraints + "(linear: " + this.solverIO.in.numLpConstraints + ", "
				+ "non-linear: " + solverIO.in.numNlpConstraints + ")\n";

		/* Decision variables are both constants (full array) and functions (to define ranges) */
		if (this.constraints == null) s = s + "-- no constraints defined --";
		else
		{
			for (Entry<String, _INTERNAL_ConstraintArray> entry : this.constraints.entrySet())
			{
				_INTERNAL_ConstraintArray cs = entry.getValue();
				if (cs == null)
					s = s + " -- null constraint set --";
				else
					s = s + cs.getName() + ": " + cs.getOriginalExpression() + " (" + cs.getNumScalarConstraints() + " scalar constraints, indexes: "
							+ "" + cs.getIndex_0_FirstConstraint() + "->" + cs.getIndex_0_LastConstraint() + ")\n";
			}
		}

		s = s + "----------- LP SCALAR CONSTRAINTS: number = " + this.lpNumScalarConstraints + "\n";
		/* Decision variables are both constants (full array) and functions (to define ranges) */
		for (int cont = 0; cont < this.lpNumScalarConstraints; cont++)
		{
			double lb = solverIO.in.constraintLowerBound.get(cont);
			double ub = solverIO.in.constraintUpperBound.get(cont);
			_INTERNAL_ConstraintArray cs = getConstraintSetInfo(cont);
			s = s + cs.getName() + ":  " + ((lb == -Double.MAX_VALUE) ?
					"-inf" :
					lb) + " <= " + solverIO.in.lhsMinusRhsAccumulatedConstraint.getAffineExpression().toStringLinearScalarExpression(cont) + " <= "
					+ (
					(ub == Double.MAX_VALUE) ?
							"inf" :
							ub) + "\n";
		}

		s = s + "----------- SOLUTION: \n";
		if (!solverIO.problemAlreadyAttemptedTobeSolved)
		{
			s = s + " -- NOT SOLVED YET --";
			return s;
		}

		s = s + "* statusCode = " + solverIO.out.statusCode + "\n";
		s = s + "* statusMessage = " + solverIO.out.statusMessage + "\n";
		s = s + "* foundFeasibleSolution = " + solverIO.out.solutionIsFeasible + "\n";
		s = s + "* solutionIsOptimal = " + solverIO.out.solutionIsOptimal + "\n";
		s = s + "* feasibleSolutionDoesNotExist = " + solverIO.out.feasibleSolutionDoesNotExist + "\n";
		s = s + "* foundUnboundedSolution = " + solverIO.out.foundUnboundedSolution + "\n";
		if (solverIO.out.solutionIsFeasible) s = s + "* primal cost = " + solverIO.out.primalCost + "\n";
		s = s + "----------- PRIMAL SOLUTION: \n";

		if (solverIO.out.solutionIsFeasible)
		{
			for (Entry<String, _INTERNAL_DecisionVariableArray> entry : decisionVariables.entrySet())
			{
				String key = entry.getKey();
				_INTERNAL_DecisionVariableArray var = entry.getValue();
				int firstVarId = var.getFirstVarId();
				int lastVarId = var.getLastVarId();
				int numDecVar = var.getNumDecVar();
				s = s + var + ", Values = " + Arrays.toString(solverIO.out.primalSolution.viewPart(firstVarId, numDecVar).toArray()) + "\n";
			}
		}
		return s;
	}

	/************ PACKAGE METHODS ********************/
	_INTERNAL_DecisionVariableArray getDecisionVariable(String name)
	{
		return decisionVariables.get(name);
	}

	_INTERNAL_DecisionVariableArray getVarInfo(int varIdThisCoef)
	{
		/* Decision variables are both constants (full array) and functions (to define ranges) */
		Iterator it = decisionVariables.entrySet().iterator();
		while (it.hasNext())
		{
			Map.Entry pairs = (Map.Entry) it.next();
			String name = (String) pairs.getKey();
			_INTERNAL_DecisionVariableArray dv = (_INTERNAL_DecisionVariableArray) pairs.getValue();
			if ((dv.getFirstVarId() <= varIdThisCoef) && (dv.getLastVarId() >= varIdThisCoef)) return dv;
		}
		return null;
	}

	private Map<String, Object> createParametersMapInitilizingDefaults(String solverName, Object[] paramValuePairs, _INTERNAL_SolverIO solverIO)
	{
		int numParameters = paramValuePairs.length / 2;
		if ((((double) paramValuePairs.length) / 2) != numParameters) throw new JOMException("A parameter has not assigned its value");
		Map<String, Object> params = new HashMap<String, Object>();
		for (int contParam = 0; contParam < numParameters; contParam++)
		{
			Object paramName = paramValuePairs[contParam * 2];
			Object value = paramValuePairs[contParam * 2 + 1];
			if (!(paramName instanceof String)) throw new JOMException("Parameter name must be a String");
			params.put((String) paramName, value);
		}

		if (!solverName.equalsIgnoreCase("glpk") && !solverName.equalsIgnoreCase("mipcl") && !solverName.equalsIgnoreCase("ipopt") && !solverName.equalsIgnoreCase("cplex") && !solverName.equalsIgnoreCase("xpress"))
			throw new JOMException("Unknown solver name");

		/* Check if the chosen solver can solve the type of problem */
		if (solverName.equalsIgnoreCase("cplex") && !solverIO.isLinearProblem()) throw new JOMException("CPLEX cannot solve non linear problems");
		if (solverName.equalsIgnoreCase("glpk") && !solverIO.isLinearProblem()) throw new JOMException("GLPK cannot solve non linear problems");
		if (solverName.equalsIgnoreCase("ipopt") && solverIO.in.hasIntegerVariables) throw new JOMException("IPOPT cannot solve problems with integer variables");
		if (solverName.equalsIgnoreCase("xpress") && !solverIO.isLinearProblem()) throw new JOMException("XPRESS cannot solve non linear problems");

		/* Default COMMON: String solverLibraryName */
		/* A "" in solverLibraryName is equivalent to not having this parameter */
		String val_solverLibraryName = (String) params.get("solverLibraryName");
		if (val_solverLibraryName != null) if (val_solverLibraryName.trim().equalsIgnoreCase("")) params.remove("solverLibraryName");
		if (!params.containsKey("solverLibraryName"))
		{
			if (solverName.equalsIgnoreCase("glpk"))
			{
				String os = System.getProperty("os.name");
				if (os.startsWith("Windows")) params.put("solverLibraryName", "glpk.dll");
				else params.put("solverLibraryName", "libglpk"); // Linux
			} else if (solverName.equalsIgnoreCase("ipopt"))
			{
				String os = System.getProperty("os.name");
				if (os.startsWith("Windows")) params.put("solverLibraryName", "ipopt.dll");
				else params.put("solverLibraryName", "libipopt"); // Linux
			} else if (solverName.equalsIgnoreCase("cplex"))
			{
				String os = System.getProperty("os.name");
				if (os.startsWith("Windows")) params.put("solverLibraryName", "cplex.dll");
				else params.put("solverLibraryName", "cplex"); // Linux
			} else if (solverName.equalsIgnoreCase("xpress"))
				params.put("solverLibraryName", "xpauth.xpr"); // default name of the license file in XPRESS
		}
//		/* Default COMMON: String maxSolverTimeInSeconds */
//		/* A "" in maxSolverTimeInSeconds is equivalent to not having this parameter */
//		Double val_maxSolverTimeInSeconds = params.get("maxSolverTimeInSeconds") == null ?
//				null :
//				((Number) params.get("maxSolverTimeInSeconds")).doubleValue();
//		if (val_maxSolverTimeInSeconds != null)
//			if (solverName.equalsIgnoreCase("glpk"))
//			{
//				final Integer numMilisecondsTimeLimit = (int) Math.ceil(1000 * val_maxSolverTimeInSeconds);
//				if (val_maxSolverTimeInSeconds > 0) params.put("tm_lim", numMilisecondsTimeLimit);
//				params.remove("maxSolverTimeInSeconds");
//			} else if (solverName.equalsIgnoreCase("ipopt"))
//			{
//				if (val_maxSolverTimeInSeconds > 0) params.put("max_cpu_time", val_maxSolverTimeInSeconds);
//				params.remove("maxSolverTimeInSeconds");
//			} else if (solverName.equalsIgnoreCase("cplex"))
//			{
//				if (val_maxSolverTimeInSeconds > 0) params.put("" + _JNA_CPLEX.CPX_PARAM_TILIM, val_maxSolverTimeInSeconds);
//				params.remove("maxSolverTimeInSeconds");
//			} else if (solverName.equalsIgnoreCase("xpress"))
//			{
//				if (val_maxSolverTimeInSeconds > 0) params.put("" + XPRSconstants.MAXTIME, (int) -Math.ceil(val_maxSolverTimeInSeconds)); // negative, if not it waits until first solution is found
//				params.remove("maxSolverTimeInSeconds");
//			}


		/* Default GLPK SPECIFIC: String glpkSolverType */
		if (solverName.equalsIgnoreCase("glpk"))
		{
			if (!params.containsKey("msg_lev")) params.put("msg_lev", 0);
			if (!params.containsKey("glpkSolverType"))
				if (!solverIO.in.hasIntegerVariables)
					//params.put("glpkSolverType", "interior-point");
					params.put("glpkSolverType", "simplex");
		}

		/* Default IPOPT SPECIFIC: print_level , derivative_test */
		if (solverName.equalsIgnoreCase("ipopt"))
		{
			if (!params.containsKey("print_level")) params.put("print_level", 0); // print_level defaults to 0
			if (!params.containsKey("derivative_test")) params.put("derivative_test", "first-order"); // derivative_test defaults to first-order
		}

		/* Default CPLEX SPECIFIC: print_level , derivative_test */
		if (solverName.equalsIgnoreCase("cplex"))
		{
			/* Print level defaults to 0 */
			String key_tuningDisplay = new Integer(_JNA_CPLEX.CPX_PARAM_TUNINGDISPLAY).toString();
			String key_screenOnOff = new Integer(_JNA_CPLEX.CPX_PARAM_SCRIND).toString();

			if (!params.containsKey(key_tuningDisplay)) params.put(key_tuningDisplay, new Integer(0)); // print_level defaults to 0
			if (!params.containsKey(key_screenOnOff))
				if ((Integer) params.get(key_tuningDisplay) != 0)
					params.put(key_screenOnOff, new Integer(_JNA_CPLEX.CPX_ON));
				else
					params.put(key_screenOnOff, new Integer(_JNA_CPLEX.CPX_OFF));
		}

		return params;
	}

	private boolean objectiveFunctionOrConstraintsAdded(){ return ((this.objectiveFunction != null) || (this.constraints.size() != 0));}

	private class TimeInfo
	{
		List<Long> accumTimeSettingDVs;
		long       accumTimeSettingObj;
		long       accumTimeSolvingTheProblem;
		long       creationTime;
		List<Long> timeSettingConstraints;

		TimeInfo(){ reset(); }

		@Override
        public String toString()
		{
			double secsSettingDVs = 0;
			for (Long num : accumTimeSettingDVs) secsSettingDVs += num * 1E-9;
			double secsSettingCons = 0;
			for (Long num : timeSettingConstraints) secsSettingCons += num * 1E-9;
			String s = "Setting DVs (secs): " + secsSettingDVs + "\n";
			s = s + "Setting Obj (secs): " + accumTimeSettingObj * 1e-9 + "\n";
			s = s + "Setting Constraints (secs): " + secsSettingCons + ": [";
			for (Long num : timeSettingConstraints) s = s + num * 1E-9 + " ";
			s = s + "] \n";
			s = s + "Solving the problem (secs): " + accumTimeSolvingTheProblem * 1e-9 + "\n";
			s = s + "Elapsed time since the object was created: " + (System.nanoTime() - creationTime) * 1E-9 + "\n";
			return s;
		}

		void reset()
		{
			creationTime = System.nanoTime();
			accumTimeSettingDVs = new ArrayList<Long>();
			accumTimeSettingObj = 0;
			timeSettingConstraints = new ArrayList<Long>();
			accumTimeSolvingTheProblem = 0;
		}
	}



}
