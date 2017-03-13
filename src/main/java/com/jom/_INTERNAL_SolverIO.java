/**
 * Copyright (c) 2015 Pablo Pavon Mari�o.
 * All rights reserved. This program and the accompanying materials
 * are made available under the terms of the GNU Lesser Public License v2.1
 * which accompanies this distribution, and is available at
 * http://www.gnu.org/licenses/lgpl.html
 * <p>
 * Contributors:
 * Pablo Pavon Mari�o - initial API and implementation
 */

/**
 *
 */
package com.jom;

import cern.colt.matrix.tdouble.DoubleFactory1D;
import cern.colt.matrix.tdouble.DoubleMatrix1D;
import cern.colt.matrix.tint.IntFactory1D;
import cern.colt.matrix.tint.IntMatrix1D;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Map.Entry;

/** @author Pablo */
class _INTERNAL_SolverIO
{
	InputInfo in;

	/* problem set up. separate lp constraints from nlp constraints */

	OutputInfo out;
	boolean problemAlreadyAttemptedTobeSolved;
	private OptimizationProblem model;

	_INTERNAL_SolverIO(OptimizationProblem model, Expression objectiveFunction, boolean toMinimize, HashMap<String, _INTERNAL_DecisionVariableArray>
			decisionVariables, HashMap<String, _INTERNAL_ConstraintArray> constraints, HashMap<String, DoubleMatrixND> initialSolution)
	{
		if (objectiveFunction == null) throw new JOMException("An objective function must be defined");

		this.model = model;
		this.problemAlreadyAttemptedTobeSolved = false;
		this.in = new InputInfo();
		this.out = new OutputInfo();

		this.in.toMinimize = toMinimize;
		this.in.objectiveFunction = objectiveFunction;
		this.in.numDecVariables = model.getNumScalarDecisionVariables();
		this.in.numConstraints = model.getNumLinearScalarConstraints() + model.getNumNonLinearScalarConstraints();
		this.in.numLpConstraints = model.getNumLinearScalarConstraints();
		this.in.numNlpConstraints = model.getNumNonLinearScalarConstraints();

		this.in.primalSolutionLowerBound = DoubleFactory1D.dense.make(model.getNumScalarDecisionVariables());
		this.in.primalSolutionUpperBound = DoubleFactory1D.dense.make(model.getNumScalarDecisionVariables());
		this.in.primalInitialSolution = DoubleFactory1D.dense.make(model.getNumScalarDecisionVariables(), 0.0);
		this.in.primalSolutionIsInteger = IntFactory1D.dense.make(model.getNumScalarDecisionVariables());
		this.in.hasIntegerVariables = false;
		for (Entry<String, _INTERNAL_DecisionVariableArray> entry : decisionVariables.entrySet())
		{
			_INTERNAL_DecisionVariableArray var = entry.getValue();
			boolean isInteger = var.getIsInteger();
			if (isInteger) this.in.hasIntegerVariables = true;
			DoubleMatrixND x_l = var.get_x_l();
			DoubleMatrixND x_u = var.get_x_u();
			int firstVarId = var.getFirstVarId();
			int numScalarVariablesThisDV = var.getNumDecVar();

			in.primalSolutionIsInteger.viewPart(firstVarId, numScalarVariablesThisDV).assign(isInteger ? 1 : 0);
			in.primalSolutionLowerBound.viewPart(firstVarId, numScalarVariablesThisDV).assign(x_l.elements().toArray()); // values are copied
			in.primalSolutionUpperBound.viewPart(firstVarId, numScalarVariablesThisDV).assign(x_u.elements().toArray()); // values are copied
			DoubleMatrixND initialSolutionThisDV = initialSolution.get(var.getName());
			if (initialSolutionThisDV != null)
			{
				in.primalInitialSolution.viewPart(firstVarId, numScalarVariablesThisDV).assign(initialSolutionThisDV.elements().toArray());
			}
		}

		this.in.lhsMinusRhsAccumulatedConstraint = null;
		this.in.constraintLowerBound = DoubleFactory1D.sparse.make(model.getNumLinearScalarConstraints() + model.getNumNonLinearScalarConstraints());
		this.in.constraintUpperBound = DoubleFactory1D.sparse.make(model.getNumLinearScalarConstraints() + model.getNumNonLinearScalarConstraints());

		/* Create the constraints matrix */
		for (Entry<String, _INTERNAL_ConstraintArray> entry : constraints.entrySet())
		{
			_INTERNAL_ConstraintArray constraintSet = entry.getValue();
			Expression exp = constraintSet.getLhsMinusRhs();
			// reshape the expression to a column vector
			exp.reshape(new int[]{exp.getNumScalarExpressions(), 1});
			int firstConstraintId = constraintSet.getIndex_0_FirstConstraint();
			int numConstraints = constraintSet.getNumScalarConstraints();
			DoubleMatrix1D lowerBound = (constraintSet.getConnector().equals(_INTERNAL_ConstraintArray.SYMBOL_LESSEQUAL)) ?
					DoubleFactory1D.dense.make(numConstraints, -Double.MAX_VALUE) :
					DoubleFactory1D.dense.make(numConstraints, 0);
			DoubleMatrix1D upperBound = (constraintSet.getConnector().equals(_INTERNAL_ConstraintArray.SYMBOL_GREATEREQUAL)) ?
					DoubleFactory1D.dense.make(numConstraints, Double.MAX_VALUE) :
					DoubleFactory1D.dense.make(numConstraints, 0);
			in.constraintLowerBound.viewPart(firstConstraintId, numConstraints).assign(lowerBound);
			in.constraintUpperBound.viewPart(firstConstraintId, numConstraints).assign(upperBound);
			if (this.in.lhsMinusRhsAccumulatedConstraint == null)
				this.in.lhsMinusRhsAccumulatedConstraint = exp;
			else
				this.in.lhsMinusRhsAccumulatedConstraint = new _FUNCTION_APPENDROWS(model, this.in.lhsMinusRhsAccumulatedConstraint, exp);
		}

		if (in.lhsMinusRhsAccumulatedConstraint == null)
		{
			if (in.numConstraints != 0) throw new JOMException("Unexpected error");
		} else if (in.lhsMinusRhsAccumulatedConstraint.getNumScalarExpressions() != in.numConstraints) throw new JOMException("Unexpected error");
	}

	@Override
	public String toString()
	{
		String s = "";
		s += "** SOLVER IO: outputSolutionSet: " + problemAlreadyAttemptedTobeSolved + "\n";
		s += "     INPUT STRUCT -----\n";
		s += "        toMinimize: " + in.toMinimize + "\n";
		s += "        numDecVariables: " + in.numDecVariables + "\n";
		s += "        numConstraints (LP / NLP): " + in.numConstraints + "(" + in.numLpConstraints + "," + in.numNlpConstraints + ")\n";
		s += "        objectiveFunction: " + in.objectiveFunction + "\n";
		s += "        primalSolutionLowerBound: " + Arrays.toString(in.primalSolutionLowerBound.toArray()) + "\n";
		s += "        primalSolutionUpperBound: " + Arrays.toString(in.primalSolutionUpperBound.toArray()) + "\n";
		s += "        primalSolutionIsInteger: " + Arrays.toString(in.primalSolutionIsInteger.toArray()) + "\n";
		if (!in.lhsMinusRhsAccumulatedConstraint.isLinear())
		{
			s += "        lhsMinusRhsAccumulatedConstraint: " + in.lhsMinusRhsAccumulatedConstraint + "\n";
			s += "        constraintLowerBound: " + Arrays.toString(in.constraintLowerBound.toArray()) + "\n";
			s += "        constraintUpperBound: " + Arrays.toString(in.constraintUpperBound.toArray()) + "\n";
		} else
		{
			s += "        lhsMinusRhsAccumulatedConstraint: " + "\n";
			_INTERNAL_AffineExpressionCoefs coefs = in.lhsMinusRhsAccumulatedConstraint.getAffineExpression();
			//			s += DoubleFactory2D.dense.make(coefs.toArray()) + "\n";
			OptimizationProblem model = in.lhsMinusRhsAccumulatedConstraint.getModel();
			for (int contCell = 0; contCell < in.lhsMinusRhsAccumulatedConstraint.getNumScalarExpressions(); contCell++)
				s = s + in.constraintLowerBound.get(contCell) + " <= " + coefs.toStringLinearScalarExpression(contCell) + " <= " + in
						.constraintUpperBound.get(contCell) + "\n";
		}
		s += "        hasIntegerVariables: " + in.hasIntegerVariables + "\n";
		if (!problemAlreadyAttemptedTobeSolved)
		{
			s += "     OUTPUT STRUCT --> NOT FILLED YET -----\n";
			return s;
		}
		s += "        statusCode: " + out.statusCode + "\n";
		s += "        statusMessage: " + out.statusMessage + "\n";
		s += "        foundFeasibleSolution: " + out.solutionIsFeasible + "\n";
		s += "        solutionIsOptimal: " + out.solutionIsOptimal + "\n";
		s += "        feasibleSolutionDoesNotExist: " + out.feasibleSolutionDoesNotExist + "\n";
		s += "        foundUnboundedSolution: " + out.foundUnboundedSolution + "\n";
		s += "        primalCost: " + out.primalCost + "\n";
		s += "        primalSolution: " + Arrays.toString(out.primalSolution.toArray()) + "\n";
		s += "        multiplierOfLowerBoundConstraintToPrimalVariables: " + ((out.multiplierOfLowerBoundConstraintToPrimalVariables != null) ?
				Arrays.toString(out.multiplierOfLowerBoundConstraintToPrimalVariables.toArray()) :
				"null") + "\n";
		s += "        multiplierOfUpperBoundConstraintToPrimalVariables: " + ((out.multiplierOfUpperBoundConstraintToPrimalVariables != null) ?
				Arrays.toString(out.multiplierOfUpperBoundConstraintToPrimalVariables.toArray()) :
				"null") + "\n";
		s += "        primalValuePerConstraint: " + Arrays.toString(out.primalValuePerConstraint.toArray()) + "\n";
		s += "        multiplierOfConstraint: " + ((out.multiplierOfConstraint != null) ?
				Arrays.toString(out.multiplierOfConstraint.toArray()) :
				"null") + "\n";

		return s;
	}

	void checkFeasibility()
	{
		final double ABS_TOLERANCE = 1e-2;
		if (!problemAlreadyAttemptedTobeSolved) throw new JOMException("The problem has not been solved");

		/* Check feasibility of the solution --> primal variable limits */
		for (int cont = 0; cont < in.numDecVariables; cont++)
		{
			if ((out.primalSolution.get(cont) < in.primalSolutionLowerBound.get(cont) - ABS_TOLERANCE))
				throw new JOMException("Check feasibility error. Decision variable ID " + cont + " has value: " + out.primalSolution.get(cont) + ", "
						+ "which is outside the feasibility limits [" + in.primalSolutionLowerBound.get(cont) + "," + in.primalSolutionUpperBound
						.get(cont) + "]");
			if ((out.primalSolution.get(cont) > in.primalSolutionUpperBound.get(cont) + ABS_TOLERANCE))
				throw new JOMException("Check feasibility error. Decision variable ID " + cont + " has value: " + out.primalSolution.get(cont) + ", "
						+ "which is outside the feasibility limits [" + in.primalSolutionLowerBound.get(cont) + "," + in.primalSolutionUpperBound
						.get(cont) + "]");
		}
		/* Check feasibility of the solution --> constraint limits */
		DoubleMatrix1D c_constraintValues = (in.numConstraints == 0) ?
				null :
				in.lhsMinusRhsAccumulatedConstraint.evaluate_internal(out.primalSolution.toArray()).elements();
		if (in.numConstraints != 0)
		{
			for (int cont = 0; cont < in.numConstraints; cont++)
			{
				if ((c_constraintValues.get(cont) < in.constraintLowerBound.get(cont) - ABS_TOLERANCE))
					throw new JOMException("Check feasibility error. Constraint number " + cont + " has a value (left hand side minus right hand "
							+ "side): " + c_constraintValues.get(cont) + ", which is outside the feasibility limits [" + in.constraintLowerBound.get
							(cont) + "," + in.constraintUpperBound.get(cont) + "]");
				if ((c_constraintValues.get(cont) > in.constraintUpperBound.get(cont) + ABS_TOLERANCE))
					throw new JOMException("Check feasibility error. Constraint number " + cont + " has a value (left hand side minus right hand "
							+ "side): " + c_constraintValues.get(cont) + ", which is outside the feasibility limits [" + in.constraintLowerBound.get
							(cont) + "," + in.constraintUpperBound.get(cont) + "]");
			}
		}

		/* Check the primal cost is ok */
		final double c_primalCostFromPrimalSolution = in.objectiveFunction.evaluate_internal(out.primalSolution.toArray()).get(0);
		if (Math.abs(out.primalCost - c_primalCostFromPrimalSolution) > ABS_TOLERANCE)
		{
			System.out.println("Warning error. The primal cost returned by the solver does not match the objetive function evaluated in the returned"
					+ " optimal point. Solver primal cost: " + out.primalCost + ". Evaluated point: " + c_primalCostFromPrimalSolution + ". This "
					+ "happens sometimes in IPOPT solver, so what I do here is set the primal cost as the one coming from evaluating the primal "
					+ "solution returned, and drop the optimum cost computed by IPOPT");
			out.primalCost = c_primalCostFromPrimalSolution;
			//throw new JOMException ("Check error. The primal cost returned by the solver does not match the objetive function evaluated in the
			// returned optimal point. Solver primal cost: " + out.primalCost + ". Evaluated point: " + c_primalCost);
		}
	}

	void checkOptimality()
	{
		final double ABS_TOLERANCE = 1e-2;

		if (!problemAlreadyAttemptedTobeSolved) throw new JOMException("The problem has not been solved");

		/* Rest of the checks do not apply to integer problems */
		if (in.hasIntegerVariables) return;

		/* Check the complementary slackness in constraints */
		DoubleMatrix1D c_constraintValues = (in.numConstraints == 0) ?
				null :
				in.lhsMinusRhsAccumulatedConstraint.evaluate_internal(out.primalSolution.toArray()).elements();
		if (in.numConstraints != 0)
		{
			for (int cont = 0; cont < in.numConstraints; cont++)
			{
				double slack_inf = Math.abs(c_constraintValues.get(cont) - in.constraintLowerBound.get(cont));
				double slack_sup = Math.abs(c_constraintValues.get(cont) - in.constraintUpperBound.get(cont));
				double mult = Math.abs(out.multiplierOfConstraint.get(cont));
				if ((mult > ABS_TOLERANCE) && (slack_inf > ABS_TOLERANCE) && (slack_sup > ABS_TOLERANCE))
					throw new JOMException("Check error. Complementary slackness conditions in constrain number " + cont + " are not satisfied. "
							+ "Multiplier: " + mult + ". Slack with lower bound: " + slack_inf + ". Slack with upper bound: " + slack_sup);
			}
		}
		/* Check the complementary slackness in variable limits */
		for (int cont = 0; cont < in.numDecVariables; cont++)
		{
			double slack_inf = Math.abs(out.primalSolution.get(cont) - in.primalSolutionLowerBound.get(cont));
			double slack_sup = Math.abs(out.primalSolution.get(cont) - in.primalSolutionUpperBound.get(cont));
			double mult_inf = Math.abs(out.multiplierOfLowerBoundConstraintToPrimalVariables.get(cont));
			double mult_sup = Math.abs(out.multiplierOfUpperBoundConstraintToPrimalVariables.get(cont));

			if ((mult_inf > ABS_TOLERANCE) && (slack_inf > ABS_TOLERANCE))
				throw new JOMException("Check error. Complementary slackness condition in lower bound constraint to decision variable number " +
						cont + " is not satisfied. Decision variable value: " + out.primalSolution.get(cont) + ". Lower bound: " + in
						.primalSolutionLowerBound.get(cont) + ". Slackness: " + slack_inf + ". Multiplier: " + mult_inf);
			if ((mult_sup > ABS_TOLERANCE) && (slack_sup > ABS_TOLERANCE))
				throw new JOMException("Check error. Complementary slackness condition in upper bound constraint to decision variable number " +
						cont + " is not satisfied. Decision variable value: " + out.primalSolution.get(cont) + ". Upper bound: " + in
						.primalSolutionUpperBound.get(cont) + ". Slackness: " + slack_sup + ". Multiplier: " + mult_sup);
		}

	}

	boolean isLinearProblem()
	{
		if (this.in.lhsMinusRhsAccumulatedConstraint == null)
			return this.in.objectiveFunction.isLinear();
		else
			return (this.in.lhsMinusRhsAccumulatedConstraint.isLinear() && this.in.objectiveFunction.isLinear());
	}

	private String printConstraintInfo(int contraintId)
	{
		if (this.in.lhsMinusRhsAccumulatedConstraint.isLinear())
			return this.in.lhsMinusRhsAccumulatedConstraint.affineExp.toStringLinearScalarExpression(contraintId);

		_INTERNAL_ConstraintArray cs = this.model.getConstraintSetInfo(contraintId);
		int cellWithinContraintArray = contraintId - cs.getIndex_0_FirstConstraint();
		return "Problem with non-linear constraints. Cell: " + cellWithinContraintArray + " within expression: " + cs.getOriginalExpression();
	}

	final class InputInfo
	{
		DoubleMatrix1D constraintLowerBound;
		DoubleMatrix1D constraintUpperBound;
		boolean        hasIntegerVariables;
		Expression     lhsMinusRhsAccumulatedConstraint; // not final since it is
		int            numConstraints;
		int            numDecVariables;
		int            numLpConstraints;
		int            numNlpConstraints;
		Expression     objectiveFunction;
		DoubleMatrix1D primalInitialSolution;
		IntMatrix1D    primalSolutionIsInteger;
		// augmented inside the
		// constructor
		DoubleMatrix1D primalSolutionLowerBound;
		DoubleMatrix1D primalSolutionUpperBound;
		boolean        toMinimize;
	}

	final class OutputInfo
	{
		boolean feasibleSolutionDoesNotExist;
		boolean foundUnboundedSolution;

		DoubleMatrix1D multiplierOfConstraint;
		DoubleMatrix1D multiplierOfLowerBoundConstraintToPrimalVariables;
		DoubleMatrix1D multiplierOfUpperBoundConstraintToPrimalVariables;
		double         primalCost;
		double bestOptimalityBound;

		DoubleMatrix1D primalSolution;
		DoubleMatrix1D primalValuePerConstraint;
		boolean        solutionIsFeasible;
		boolean        solutionIsOptimal;
		int            statusCode;
		String         statusMessage;
		public String toString ()
		{
			if (feasibleSolutionDoesNotExist) return "A feasible solution does nor exist";
			if (foundUnboundedSolution) return "An unbounded solution was found";
			if (!solutionIsFeasible) return "A feasible solution was not found, but may exist";
			StringBuffer s = new StringBuffer ();
			s.append("Status code: " + statusCode + ", message: " + statusMessage + String.format("\n"));
			s.append("Solution is optimal? " + solutionIsOptimal + String.format("\n"));
			s.append("Primal solution: " + Arrays.toString(primalSolution.toArray()) + String.format("\n"));
			s.append("Primal cost: " + primalCost + String.format("\n"));
			s.append("Best optimality bound: " + bestOptimalityBound + String.format("\n"));
			return s.toString();
		}
	}
}
