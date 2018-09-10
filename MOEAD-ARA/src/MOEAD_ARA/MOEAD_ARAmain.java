
package MOEAD_ARA;

import jmetal.core.Algorithm;
import jmetal.core.Operator;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.operators.crossover.CrossoverFactory;
import jmetal.operators.mutation.MutationFactory;
import jmetal.operators.selection.SelectionFactory;
import jmetal.qualityIndicator.QualityIndicator;
import jmetal.qualityIndicator.QualityIndicator_self;
import jmetal.util.Configuration;
import jmetal.util.JMException;

import java.io.IOException;
import java.util.HashMap;
import java.util.logging.FileHandler;
import java.util.logging.Logger;

public class MOEAD_ARAmain {
	public static Logger logger_; // Logger object
	public static FileHandler fileHandler_; // FileHandler object
	
	public static String[] alg = {"MOEAD_ARA"};

	public static String[] probs= {"MOP1"};

	public static void main(String[] args) throws JMException, SecurityException, IOException, ClassNotFoundException {
		Problem problem; // The problem to solve
		Algorithm algorithm; // The algorithm to use
		Operator crossover; // Crossover operator
		Operator mutation; // Mutation operator

		HashMap parameters; // Operator parameters
		String date = "test";
		int runsum = 5;

		String[] pps = probs;
		String curDir = System.getProperty("user.dir");
		String DATAFILE = curDir + "/" + "ExperimentalData";

		for (int pp = 0; pp < pps.length; pp++) {
			problem = Util_ARA.getProblem(pps[pp]);

			for (int al = 0; al < alg.length; al++) {

				String algName = alg[al];
				String PFFILE = DATAFILE + "/" + date + "/" + "PF" + "/";
				String PSFILE = DATAFILE + "/" + date + "/" + "PS" + "/";

				algorithm = Util_ARA.getAlgorithm(problem, alg[al]);
				algorithm.setInputParameter("algName", algName);
			
				algorithm.setInputParameter("populationSize", 300);
				algorithm.setInputParameter("maxEvaluations", 300000);
				
				algorithm.setInputParameter("dataDirectory", curDir + "/DirectionVector");

				

					algorithm.setInputParameter("finalSize", 300); // used by//
																	// MOEAD_DRA
					algorithm.setInputParameter("T", 20);
					algorithm.setInputParameter("delta", 0.9);
					algorithm.setInputParameter("nr", 2);

					// Crossover operator
					parameters = new HashMap();
					parameters.put("CR", 1.0);
					parameters.put("F", 0.5);
					crossover = CrossoverFactory.getCrossoverOperator("DifferentialEvolutionCrossover", parameters);

					// Mutation operator
					parameters = new HashMap();
					parameters.put("probability", 1.0 / problem.getNumberOfVariables());
					parameters.put("distributionIndex", 20.0);
					mutation = MutationFactory.getMutationOperator("PolynomialMutation", parameters);

					algorithm.addOperator("crossover", crossover);
					algorithm.addOperator("mutation", mutation);

				

				for (int iter = 0; iter < runsum; iter++) {
					// Execute the Algorithm
		//			long initTime = System.currentTimeMillis();
					Util_ARA.initPrint(iter, date, algName, problem);
					SolutionSet population = algorithm.execute();
		//			long estimatedTime = System.currentTimeMillis() - initTime;

					String PFPATH = PFFILE + "/" + "pf_" + algName + "_" + problem.getName() + "_run" + iter + "_"+ "_end.txt";
					population.printObjectivesToFile(PFPATH);
					String PSPATH = PSFILE + "/" + "ps_" + algName + "_" + problem.getName() + "_run" + iter + "_"+ "_end.txt";
					population.printVariablesToFile(PSPATH);
				}
			}
		}
	} // main

} // MOEAD_main
