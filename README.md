# JOM
[![](https://jitpack.io/v/girtel/JOM.svg)](https://jitpack.io/#girtel/JOM)

Java Optimization Modeler


JOM (Java Optimization Modeler) is a free and open-source Java library targeted to allowing Java programs modeling and solving optimization problems, defining their input parameters, decision variables, objective function, and constraints. Mathematical expressions combining decision variables, constant parameters, operators and functions, are defined using a simple MATLAB-like syntax. JOM allows handling multidimensional arrays of variables and constants. This enables, for instance, defining arrays of constraints in a single line of code.

To solve the problem, JOM can interface with installed solvers. Current JOM version can interface with GLPK and MIPCL (free), CPLEX and XPRESS (commercial) solvers for mixed integer linear problems, and IPOPT (free) for non-linear differentiable problems.

JOM has been thought as a tool to assist the teaching and research activities which are focused on optimization, or where optimization problems are solved. Students/researchers can focus on solving optimization problems and analyzing the solutions obtained, getting rid of the burden of interfacing to the solvers, and making use of the flexibility that Java provides to handle the problem data. JOM is just accessed from standard Java programs, calling the methods in the JOM class OptimizationProblem. In addition, the JOM modeling language is much simpler than (although not as powerful as) commercial alternatives for modeling as AMPL or GAMS.

Visit http://www.net2plan.com/jom

# Building

The library is currently build under Maven. Run the following goal to build the library:

`clean package`

The result is a JOM-${VERSION} package contained under the _target/assembly_ folder of the project.

