/**
 * @file Function.java  Interface for the `Function` class.
 * <table border="2">
 * <tr>
 *   <td>  <b>Author</b> </td>  <td> Nicolai N. Pisaruk </td>
 * </tr> 
 * <tr>
 *  <td> <b> e-mail </b> </td>  <td> nicolaipisaruk@gmail.com </td>
 * </tr>
 * <tr>
 *  <td> <b>home page</b> </td> <td> wwww.mipcl-cpp.appspot.com </td>                |
 * </tr>
 * <tr>
 *    <td colspan="2"> <b>&copy; 2017</b> Nicolai N. Pisaruk </td>
 * </tr>
 * </table>
 */

package by.bsu.JVmipshell;

import java.util.List;
import by.bsu.JVmipshell.Var;
import by.bsu.JVmipshell.MIPshell;


/**
 * This class represents points in the plane.
 */
class Point
{
    /**
     * Point coordinates.
     */
	double x, y;

	public Point(double x, double y)
	{
		this.x = x;
		this.y = y;
	}

	public double getx()
	{
		return x;
	}

	public double gety()
	{
		return y;
	}

}

/**
 * The 'Function' objects represent piece-wise linear functions _y = f(x)_.
 */
public class Function
{
    /**
     * Function argument.
     */
	Var x;

    /**
     *  Function value.
     */
        Var y;

    /**
     *  List of pairs of real numbers: function break points.
     */
        List<Point> points;

    /**
     * The handle that uniquely identifies the function constraint.
     */
        int hd;

    /**
     * Set the handle of `this` variable.
     *
     * @param[in] hd new handle for `this` variable.
     */
	public void sethd(int hd)
	{
		this.hd = hd;
	}

    /**
     * The constructor that initializes a `Function` object.
     *
     * @param[in] prob problem to which `this` function be added. 
     * @param[in] hd function handle;
     * @param[in] y  function value;
     * @param[in] x function argument;
     * @param[in] points function break points.
     */
	public Function(MIPshell prob, int hd, Var y, Var x, List<Point> points)
	{
 		this.hd = hd;
		this.y = y;
		this.x = x;
		this.points = points;
		prob.addFunc(this);
 	}

    /**
     * @return argument variable.
     */
	public Var getx()
	{
		return x;
	}

    /**
     * @return function value variable.
     */
	public Var gety()
	{
		return y;
	}

    /**
     * @return reference to the list of breakpoints.
     */
	public List<Point> getpoints()
	{
		return points;
	}
}

