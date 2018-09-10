/**
 * MOEAD_DRA.java
 * 
 * This is the main implementation of MOEA/D-DRA.
 * 
 * Reference:
 * 
 * 		Qingfu Zhang, Wudong Liu, Hui Li, The performance of a new version of 
 * 		MOEA/D on CEC09 unconstrained MOP test instances. IEEE Congress on 
 * 		Evolutionary Computation 2009: 203-208
 */

package MOEAD_ARA;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;
import java.util.StringTokenizer;
import jmetal.util.*;

import java.util.Vector;

import jmetal.core.*;
import jmetal.util.PseudoRandom;

public class MOEAD_ARA extends Algorithm {

	private int populationSize_;
	private SolutionSet population_; // Population repository
	private Solution[] savedValues_; // Individual repository

	int T_; // Neighborhood size
	int nr_; // Maximal number of solutions replaced by each child solution
	double delta_; // Probability that parent solutions are selected from
					// neighborhood
	int evaluations_; // Counter for the number of function evaluations
	int maxEvaluations_;

	double[] z_; // Z vector (ideal point)
	int[][] neighborhood_; // Neighborhood matrix
	double[][] lambda_; // Lambda vectors

	private double[] utility_;
	private int[] frequency_;

	Solution[] indArray_;
	String functionType_;

	Operator crossover_;
	Operator mutation_;

	String dataDirectory_;

	/**
	 * Constructor
	 * 
	 * @param Problem
	 *            to solve
	 */
	public MOEAD_ARA(Problem problem) {
		super(problem);

		// functionType_ = "_TCHE1";
		functionType_ = "_TCHE2";
	} // DMOEA

	public SolutionSet execute() throws JMException, ClassNotFoundException {

		evaluations_ = 0;
		maxEvaluations_ = ((Integer) this.getInputParameter("maxEvaluations")).intValue();
		populationSize_ = ((Integer) this.getInputParameter("populationSize")).intValue();
		dataDirectory_ = this.getInputParameter("dataDirectory").toString();

		population_ = new SolutionSet(populationSize_);
		savedValues_ = new Solution[populationSize_];
		utility_ = new double[populationSize_];
		frequency_ = new int[populationSize_];
		for (int i = 0; i < utility_.length; i++) {
			utility_[i] = 1.0;
			frequency_[i] = 0;
		}
		indArray_ = new Solution[problem_.getNumberOfObjectives()];

		T_ = 20;
		delta_ = 0.9;
		nr_ = (int) (0.01 * populationSize_);
		//nr_ = 2;

		neighborhood_ = new int[populationSize_][T_];

		z_ = new double[problem_.getNumberOfObjectives()];
		lambda_ = new double[populationSize_][problem_.getNumberOfObjectives()];

		crossover_ = operators_.get("crossover"); // default: DE crossover
		mutation_ = operators_.get("mutation"); // default: polynomial mutation

		// STEP 1. Initialization
		// STEP 1.1. Compute Euclidean distances between weight vectors and find
		// T
		initUniformWeight();
		initNeighborhood();

		// STEP 1.2. Initialize population
		initPopulation();

		// STEP 1.3. Initialize z_
		initIdealPoint();

		int gen = 0;
		// STEP 2. Update
	//	Util_ARA.initProcessPrint();
		do {
	//		Util_ARA.processPrint(evaluations_, population_);
			int[] permutation = new int[populationSize_];
			Utils.randomPermutation(permutation, populationSize_);
			List<Integer> order = tour_selection(10);

			for (int i = 0; i < order.size(); i++) {
				int n = order.get(i);
				frequency_[n]++;

				int type;
				double rnd = PseudoRandom.randDouble();

				// STEP 2.1. Mating selection based on probability
				if (rnd < delta_) // if (rnd < realb)
				{
					type = 1; // neighborhood
				} else {
					type = 2; // whole population
				}
				Vector<Integer> p = new Vector<Integer>();
				matingSelection(p, n, 2, type);

				// STEP 2.2. Reproduction
				Solution child;
				Solution[] parents = new Solution[3];

				parents[0] = population_.get(p.get(0));
				parents[1] = population_.get(p.get(1));
				parents[2] = population_.get(n);

				// Apply DE crossover
				child = (Solution) crossover_.execute(new Object[] { population_.get(n), parents });

				// Apply mutation
				mutation_.execute(child);

				// Evaluation
				problem_.evaluate(child);
				evaluations_++;

				// STEP 2.3. Repair. Not necessary

				// STEP 2.4. Update z_
				updateReference(child);

				// STEP 2.5. Update of solutions
				updateProblem(child, n, type);

			} // for

			gen++;
			if (gen % 50 == 0) {
				comp_utility();
			}
		} while (evaluations_ < maxEvaluations_);
		/*
		 * for (int i = 0; i < populationSize_; i++) {
		 * System.out.println(frequency_[i]); }
		 */

		return population_;
	}

	/**
	 * Initialize the weight vectors, this function only can read from the
	 * existing data file, instead of generating itself.
	 * 
	 */
	public void initUniformWeight() {
		String dataFileName;
		dataFileName = "W" + problem_.getNumberOfObjectives() + "D_" + populationSize_ + ".dat";

		try {
			// Open the file
			FileInputStream fis = new FileInputStream(dataDirectory_ + "/" + dataFileName);
			InputStreamReader isr = new InputStreamReader(fis);
			BufferedReader br = new BufferedReader(isr);

			int i = 0;
			int j = 0;
			String aux = br.readLine();
			while (aux != null) {
				StringTokenizer st = new StringTokenizer(aux);
				j = 0;
				while (st.hasMoreTokens()) {
					double value = (new Double(st.nextToken())).doubleValue();
					lambda_[i][j] = value;
					j++;
				}
				aux = br.readLine();
				i++;
			}
			br.close();
		} catch (Exception e) {
			System.out
					.println("initUniformWeight: failed when reading for file: " + dataDirectory_ + "/" + dataFileName);
			e.printStackTrace();
		}
	} // initUniformWeight

	/**
	 * Update the utilities of subproblems
	 * 
	 */
	public void comp_utility() {
		double f1, f2, uti, delta;

		for (int i = 0; i < populationSize_; i++) {
			f1 = fitnessFunction(population_.get(i), lambda_[i]);
			f2 = fitnessFunction(savedValues_[i], lambda_[i]);

			delta = (f2 - f1) / f2;
			if (delta > 0.001)
				utility_[i] = 1.0;
			else {
				uti = 0.95 * (1.0 + delta / 0.001) * utility_[i];
				utility_[i] = uti < 1.0 ? uti : 1.0;
			}
			savedValues_[i] = new Solution(population_.get(i));
		}
	}

	/**
	 * Initialize the neighborhood matrix of subproblems, based on the Euclidean
	 * distances between different weight vectors
	 * 
	 */
	public void initNeighborhood() {
		int[] idx = new int[populationSize_];
		double[] x = new double[populationSize_];

		for (int i = 0; i < populationSize_; i++) {
			/* Calculate the distances based on weight vectors */
			for (int j = 0; j < populationSize_; j++) {
				x[j] = Utils.distVector(lambda_[i], lambda_[j]);
				idx[j] = j;
			}
			/* Find 'niche' nearest neighboring subproblems */
			Utils.minFastSort(x, idx, populationSize_, T_);

			for (int k = 0; k < T_; k++) {
				neighborhood_[i][k] = idx[k];
			}
		}
	} // initNeighborhood

	/**
	 * Initialize the population, random sampling from the search space
	 * 
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	public void initPopulation() throws JMException, ClassNotFoundException {
		for (int i = 0; i < populationSize_; i++) {
			Solution newSolution = new Solution(problem_);

			problem_.evaluate(newSolution);
			evaluations_++;
			population_.add(newSolution);
			savedValues_[i] = new Solution(newSolution);
		}
	} // initPopulation

	/**
	 * Initialize the ideal point, the best objective function value for each
	 * individual objective
	 * 
	 * @throws JMException
	 * @throws ClassNotFoundException
	 */
	void initIdealPoint() throws JMException, ClassNotFoundException {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			z_[i] = 1.0e+30;
			indArray_[i] = new Solution(problem_);
			problem_.evaluate(indArray_[i]);
		//	evaluations_++;
		}

		for (int i = 0; i < populationSize_; i++)
			updateReference(population_.get(i));
	} // initIdealPoint

	/**
	 * Select the mating parents, depending on the selection 'type'
	 * 
	 * @param list
	 *            : the set of the indexes of selected mating parents
	 * @param cid
	 *            : the id of current subproblem
	 * @param size
	 *            : the number of selected mating parents
	 * @param type
	 *            : 1 - neighborhood; otherwise - whole population
	 */
	public void matingSelection(Vector<Integer> list, int cid, int size, int type) {
		int ss;
		int r;
		int p;

		ss = neighborhood_[cid].length;
		while (list.size() < size) {
			if (type == 1) {
				r = PseudoRandom.randInt(0, ss - 1);
				p = neighborhood_[cid][r];
			} else {
				p = PseudoRandom.randInt(0, populationSize_ - 1);
			}
			boolean flag = true;
			for (int i = 0; i < list.size(); i++) {
				if (list.get(i) == p) // p is in the list
				{
					flag = false;
					break;
				}
			}

			if (flag) {
				list.addElement(p);
			}
		}
	} // matingSelection

	/**
	 * Tournament selection
	 * 
	 * @param depth
	 * @return
	 */
	public List<Integer> tour_selection(int depth) {

		// selection based on utility
		List<Integer> selected = new ArrayList<Integer>();
		List<Integer> candidate = new ArrayList<Integer>();

		for (int k = 0; k < problem_.getNumberOfObjectives(); k++)
			selected.add(k); // select first m weights
		for (int n = problem_.getNumberOfObjectives(); n < populationSize_; n++)
			candidate.add(n); // set of unselected weights

		while (selected.size() < (int) (populationSize_ / 5.0)) {
			int best_idd = (int) (PseudoRandom.randDouble() * candidate.size());
			int i2;
			int best_sub = candidate.get(best_idd);
			int s2;
			for (int i = 1; i < depth; i++) {
				i2 = (int) (PseudoRandom.randDouble() * candidate.size());
				s2 = candidate.get(i2);

				if (utility_[s2] > utility_[best_sub]) {
					best_idd = i2;
					best_sub = s2;
				}
			}
			selected.add(best_sub);
			candidate.remove(best_idd);
		}
		return selected;
	}

	/**
	 * Update the current ideal point
	 * 
	 * @param individual
	 */
	void updateReference(Solution individual) {
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
			if (individual.getObjective(i) < z_[i])
				z_[i] = individual.getObjective(i);
		}
	} // updateReference

	/**
	 * Update the population by the current offspring
	 * 
	 * @param indiv:
	 *            current offspring
	 * @param id:
	 *            index of current subproblem
	 * @param type:
	 *            update solutions in - neighborhood (1) or whole population
	 *            (otherwise)
	 */
	void updateProblem(Solution indiv, int id, int type) {
		int size;  
		size = population_.size();
		int[] perm = new int[size];
		Utils.randomPermutation(perm, size);
		int[] prom_ind=new int[lambda_[0].length];
		for(int j=0;j<prom_ind.length;j++){
			prom_ind[j]=j;
		}

		for (int i = 0; i < size; i++) {
			int k;
			k = perm[i];
			for(int j=0;j<prom_ind.length;j++){
				if(p_l_angle(indiv, lambda_[k])<p_l_angle(indiv, lambda_[prom_ind[j]])){
					int temp=prom_ind[j];
					prom_ind[j]=k;
					k=temp;
				}
			}	
		}
		findRegion(indiv,prom_ind);
		
	} // updateProblem
	
	void findRegion(Solution ind_new,int[] id) {
		
		 double[] point=new double[id.length];
		 double[] f=new double[id.length];
		for(int i=0;i<id.length;i++){
			f[i]=fitnessFunction(population_.get(id[i]),lambda_[id[i]]);
		}
		for(int i=0;i<lambda_[0].length;i++){
		//	double min=f[0]*lambda_[id[0]][i];
			double min=Double.POSITIVE_INFINITY;
			for(int j=0;j<id.length;j++){
				if(min>f[j]*lambda_[id[j]][i]){
					min=f[j]*lambda_[id[j]][i];
				}
			}
			point[i]=min;
		}
		double f_new;
		for(int i=0;i<id.length;i++){
			f_new = fitnessFunction(ind_new, lambda_[id[i]]);
			if(f_new<f[i]){
				double rndd = PseudoRandom.randDouble();
				double dd=1;
				if(rndd<dd){
					if(p_l_angle(ind_new, lambda_[id[i]])<point_l_angle(point,lambda_[id[i]])){
						population_.replace(id[i], new Solution(ind_new));
					}
				}else{
					population_.replace(id[i], new Solution(ind_new));
				}
			}
			
		/*	if((f_new<f[i])&&
					(p_l_angle(ind_new, lambda_[id[i]])<point_l_angle(point,lambda_[id[i]]))){
				population_.replace(id[i], new Solution(ind_new));
			}*/
		}
		
	}
//end
	
	
	double innerproduct(double[] vec1, double[] vec2) {
		double sum = 0;
		for (int i = 0; i < vec1.length; i++)
			sum += vec1[i] * vec2[i];
		return sum;
	}

	double norm_vector(Vector<Double> x) {
		double sum = 0.0;
		for (int i = 0; i < (int) x.size(); i++)
			sum = sum + x.get(i) * x.get(i);
		return Math.sqrt(sum);
	}

	double norm_vector(double[] z) {
		double sum = 0.0;
		for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
			sum = sum + z[i] * z[i];
		return Math.sqrt(sum);
	}

	/**
	 * Evaluate the fitness function by the decomposition method
	 * 
	 * @param individual:
	 *            current solution
	 * @param lambda:
	 *            : weight vector
	 * @return
	 */
	double fitnessFunction(Solution individual, double[] lambda) {
		double fitness;
		fitness = 0.0;

		if (functionType_.equals("_TCHE1")) {
			double maxFun = -1.0e+30;

			for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
				double diff = Math.abs(individual.getObjective(n) - z_[n]);

				double feval;
				if (lambda[n] == 0) {
					feval = 0.0001 * diff;
				} else {
					feval = diff * lambda[n];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for

			fitness = maxFun;
		} else if (functionType_.equals("_TCHE2")) {
			double maxFun = -1.0e+30;

			for (int i = 0; i < problem_.getNumberOfObjectives(); i++) {
				double diff = Math.abs(individual.getObjective(i) - z_[i]);

				double feval;
				if (lambda[i] == 0) {
					feval = diff / 0.000001;
				} else {
					feval = diff / lambda[i];
				}
				if (feval > maxFun) {
					maxFun = feval;
				}
			} // for
			fitness = maxFun;
		} else if (functionType_.equals("_PBI"))// if
		{
			double theta; // penalty parameter
			theta = 5.0;

			// normalize the weight vector (line segment)
			double nd = norm_vector(lambda);
			for (int i = 0; i < problem_.getNumberOfObjectives(); i++)
				lambda[i] = lambda[i] / nd;

			double[] realA = new double[problem_.getNumberOfObjectives()];
			double[] realB = new double[problem_.getNumberOfObjectives()];

			// difference between current point and reference point
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++)
				realA[n] = (individual.getObjective(n) - z_[n]);

			// distance along the line segment
			double d1 = Math.abs(innerproduct(realA, lambda));

			// distance to the line segment
			for (int n = 0; n < problem_.getNumberOfObjectives(); n++)
				realB[n] = (individual.getObjective(n) - (z_[n] + d1 * lambda[n]));
			double d2 = norm_vector(realB);

			fitness = d1 + theta * d2;
		} else {
			System.out.println("MOEAD.fitnessFunction: unknown type " + functionType_);
			System.exit(-1);
		}
		return fitness;
	} // fitnessEvaluation
	
	/**
	 * aibo 2016.11.21 at home To calculate the distance the solutions to the
	 * subproblems vector Using trigonometric function method
	 * 
	 * @param individual
	 * @param lambda
	 * @return
	 */
	double p_l_angle(Solution individual, double[] lambda) {
		double sinAngle = 0.0;

		// vector module of the vector from point z_ to point solution
		double sz_vm = 0;
		// vector module of the vector lambda( subproblem)
		double lambda_vm = 0;

		// inner product of the sz and lambda
		double sz_lambda_innerP = 0;

		double diff;
		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
			diff = Math.abs(individual.getObjective(n) - z_[n]);
			sz_lambda_innerP = sz_lambda_innerP + diff * lambda[n];
			sz_vm = sz_vm + diff * diff;
			lambda_vm = lambda_vm + lambda[n] * lambda[n];
		}
		sinAngle = Math.sqrt(1 - sz_lambda_innerP*sz_lambda_innerP / (sz_vm * lambda_vm));
		return sinAngle;
	}
	
	/**
	 * aibo 2016.11.21 at home To calculate the distance the solutions to the
	 * subproblems vector Using trigonometric function method
	 * 
	 * @param individual
	 * @param lambda
	 * @return
	 */
	double point_l_angle(double[] point, double[] lambda) {
		double sinAngle = 0.0;

		// vector module of the vector from point z_ to point solution
		double sz_vm = 0;
		// vector module of the vector lambda( subproblem)
		double lambda_vm = 0;

		// inner product of the sz and lambda
		double sz_lambda_innerP = 0;

		double diff;
		for (int n = 0; n < problem_.getNumberOfObjectives(); n++) {
	//		diff = Math.abs(individual.getObjective(n) - z_[n]);
			sz_lambda_innerP = sz_lambda_innerP + point[n] * lambda[n];
			sz_vm = sz_vm + point[n] * point[n];
			lambda_vm = lambda_vm + lambda[n] * lambda[n];
		}
		sinAngle = Math.sqrt(1 - sz_lambda_innerP*sz_lambda_innerP / (sz_vm * lambda_vm));
		return sinAngle;
	}
} // MOEAD
