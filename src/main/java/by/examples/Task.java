package by.examples;

import java.util.Collections;
import java.util.Map;
import java.util.Set;

class Machine {
	/** Machine name */
	String name;

	/** `available[t]` is number of machines available in period `t`. */
	int[] available;

	/** `cost[t]` is cost of using one machine of this type in period `t`. */
	double[] cost;

	int getMachineNum(int t)
	{return available[t];}

	double getcost(int t)
	{return cost[t];}
} // end od Machine

class Product {
	/** Product name */
	String name;

	/** Init stock */
	double initStock;

	/** Final stock */
	double finalStock;

	/** Storage capacity */
	double capacity;

	/** Demands for product */
	double[] demands;

	/** Unit production costs */
	double[] costs;

	/** Holding costs */
	double[] holdingCosts;

	/** Processing times.
         *
	 *  Keys are machines capable of producing this product;
	 *  the value of `pocTimes.get(i)` is time needed to produce one unit of this product on machine `i`.
	 */
	Map<Integer,Double> procTimes;

	/** Resources.
         *
	 *  Keys are products used in production of this product;
	 *  the value of `resources.get(k)` is number of unit of resource `k` used to produce one of this product.
	 */
	Map<Integer,Double> resources;

// Attributes 
	double getinitStock()
	{return initStock;}

	double getfinalStock()
	{return finalStock;}

	double getcapacity()
	{return capacity;}

	double[] getdemands()
	{return demands;}

	double[] getcosts()
	{return costs;}

	double[] getholdingCosts()
	{return holdingCosts;}

	Map<Integer,Double> getprocTimes()
	{return procTimes;}

	Map<Integer,Double> getresources()
	{return resources;}
} // end of Product

public class Task
{
	/** Time horizon of `T=l.length` periods, `l[t]` is length (duration) of period `t`. */
	int[] l;

	/** Array of machines. */
	Machine[] machines;

	/** Products to be produced */  
	Product[] products;

// Attributes
    /** @return number of period in the planning horizon. */
	int getHorizonLength()
	{return l.length;}

    /**
     * @param[in] t period.
     * @return length (duration) of period `t`.
     */
	int getPeriodLength(int t)
	{return l[t];}

    /**
     * @return number of machines.
     */
	int getMachineNum()
	{return machines.length;}

    /**
       @param[in] i machine type;
       @param[in] t period.
       @return number of machines of type `i` available in period `t`.
     */
	int getMachineNum(int i, int t)
	{return machines[i].getMachineNum(t);}

    /**
       @param[in] i machine type;
       @param[in] t period.
       @return cost of using one machine of type `i` in period `t`.
     */
	double getMachineCost(int i, int t)
	{return machines[i].getcost(t);}

    /**
     * @return number of products.
     */
	int getProductNum()
	{return products.length;}

    /**
     * @param[in] j product index.
     * @return initial stock of product `j`.
     */
	double getinitStock(int j)
	{return products[j].getinitStock();}

    /**
     * @param[in] j product index.
     * @return final stock of product `j`.
     */
	double getfinalStock(int j)
	{return products[j].getfinalStock();}

    /**
     * @param[in] j product index.
     * @return final stock of product `j`.
     */
	double getcapacity(int j)
	{return products[j].getcapacity();}

    /**
     * @param[in] j product index;
     * @param[in] t period.
     * @return cost of producing one unit of product `j` in period `t`.
     */
	double getcost(int j, int t)
	{return products[j].getcosts()[t];}

    /**
     * @param[in] j product index;
     * @param[in] t period.
     * @return holding cost of storing one unit of product `j` during period `t`.
     */
	double getholdingCost(int j, int t)
	{return products[j].getholdingCosts()[t];}

    /**
     * @param[in] j product index;
     * @param[in] t period.
     * @return demand for product `j` in period `t`.
     */
	double getdemand(int j, int t)
	{return products[j].getdemands()[t];}

    /**
     * @param[in] j product index.
     * @return set of machines used for producing product `j`.
     */
	Set<Integer> M(int j) {
		Map<Integer,Double> times = products[j].getprocTimes(); 
		return (times != null)? times.keySet(): Collections.emptySet();
	}

    /**
     * @param[in] j product index;
     * @param[in] i machine from `M(j)`.
     * @return processing time of one unit of product `j` on machine `i`.
     */
	double getprocTime(int j, int i)
	{return products[j].getprocTimes().get(i);}

    /**
     * @param[in] j product index.
     * @return set of products used for producing product `j`.
     */
	Set<Integer> R(int j) {
		Map<Integer,Double> res=products[j].getresources();
		return (res != null)? res.keySet(): Collections.emptySet();
	}

    /**
     * @param[in] j product index;
     * @param[in] k product from `R(j)`.
     * @return amount of product `k` used for producing one unit of product `j`.
     */
	double getResUsed(int j, int k)
	{return products[j].getresources().get(k);}
} // end of Task

