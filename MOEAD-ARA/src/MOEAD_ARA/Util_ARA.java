package MOEAD_ARA;

import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;

import jmetal.core.Algorithm;
import jmetal.core.Problem;
import jmetal.core.SolutionSet;
import jmetal.problems.ZDT.*;
import jmetal.problems.MOP.*;
import jmetal.problems.WFG.*;
import jmetal.problems.LZ09.*;
import jmetal.problems.DTLZ.*;
import jmetal.problems.cec2009Competition.*;
import jmetal.qualityIndicator.QualityIndicator;

public class Util_ARA {
	/*
	 * aibo psPrintpath 在程序运行中打印进化过程中的ps 的路径 pfPrintpath 在程序运行中打印进化过程中的pf 的路径
	 * runnum 算法独立运行的哪一次 date 运行的日期
	 */
	public static String psPrintpath;
	public static String pfPrintpath;
	public static String METRICSFILE;
	public static int runnum;
	public static String date;
	public static String algName;
	public static String problemName;
	public static Problem problem;
	public static QualityIndicator indicator;

	public static void initPrint(int num, String dateStr, String algStr, Problem prob) {
		date = dateStr;
		String curDir = System.getProperty("user.dir");
		String DATAFILE = curDir + "/" + "ExperimentalData" + "/" + date;
		pfPrintpath = DATAFILE + "/" + "PF" + "/";
		psPrintpath = DATAFILE + "/" + "PS" + "/";
		METRICSFILE = DATAFILE + "/" + "METRICS" + "/";
		runnum = num;
		algName = algStr;
		problemName = prob.getName();
		problem = prob;
	}

	public static void processPrint(int evaluations, SolutionSet population) {
	//	if (runnum == 1) {
			boolean isprint = false;
			if ("MOP6".equals(problemName) || "MOP7".equals(problemName)||
					"UF8".equals(problemName) || "UF9".equals(problemName)
					|| "UF10".equals(problemName)) {
				if (evaluations % 6000 == 0) {
					isprint = true;
				}
			}else {
				if (evaluations % 3000 == 0) {
					isprint = true;
				}
			}
			
			String curDir = System.getProperty("user.dir");
			if (isprint) {
				String PFPATH = pfPrintpath + "/pos/" + "pf_" + algName + "_" + problemName + "_run" + runnum + "_e"
						+ evaluations + ".txt";
				population.printObjectivesToFile(PFPATH);
				String PsPATH =pfPrintpath + "/pos/" + "ps_" + algName + "_" + problemName + "_run" + runnum + "_e"
						+ evaluations + ".txt";
				population.printVariablesToFile(PsPATH);

				String METRICSPATH = METRICSFILE + "/pos/" + "metr_" + algName + "_" + problem.getName() + "_run" + runnum
						+ ".txt";
				printVectorToFile(METRICSPATH, new Double[] { evaluations + 0.0,
						indicator.getIGD(population), indicator.getHypervolume(population) }, true);
			}
	//	}
	}

	
	/*public static void processPrint(int evaluations, SolutionSet population,String alg_Name,String proName) {
				boolean isprint = false;
			
				if (evaluations % 6000 == 0) {
					isprint = true;
				}
				String curDir = System.getProperty("user.dir");
				if (isprint) {
					String PFPATH = curDir+"/expdata/pf" + "/" + "pf_" + alg_Name + "_" + proName + "_run" + runnum + "_e"
							+ evaluations + ".txt";
					population.printObjectivesToFile(PFPATH);
					String PsPATH = curDir+"/expdata/ps" + "/" + "ps_" + alg_Name + "_" + proName + "_run" + runnum + "_e"
							+ evaluations + ".txt";
				//	population.printObjectivesToFile(PsPATH);
				
				}
		}
	*/
	
	public static void initProcessPrint() {
		//if (runnum == 1) {
			String curDir = System.getProperty("user.dir");
			String paretoFrontFile = curDir + "/" + "trueParetoFront/" + problemName + ".dat";
			indicator = new QualityIndicator(problem, paretoFrontFile);

			String METRICSPATH = curDir+"/expdata/metri" + "/" + "metr_" + algName + "_" + problem.getName() + "_run" + runnum
					+ ".txt";
			printVectorToFile(METRICSPATH, new String[] { "evaluations", "IGD", "HV" }, false);
	//	}
	}

	public static double getAverage(double[] d) {
		double sum = 0;
		for (int i = 0; i < d.length; i++) {
			sum = sum+d[i];
		}
		return (sum / d.length);
	}
	
	 public static double getStandardDevition(double[] d){
	        double sum = 0;
	        double mean=getAverage(d);
	        for(int i = 0;i < d.length;i++){
	            sum += (d[i] -mean) * (d[i] -mean);
	        }
	        return (Math.sqrt(sum / (d.length)));
	    }
	 
	 
	  public static void printResult_mean_sd(String path, double[][] mean,double[][] sd,String[] probName,String[] algName,boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path,append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);
		      for(int i=0;i<mean.length;i++){
		    	  bw.write(probName[i]);
		    	  for(int j=0;j<mean[i].length;j++){
		    		  bw.write("  "+mean[i][j]+"("+sd[i][j]+")");
		    	  }
		    	  bw.newLine();
		      }
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }
	  
	  public static void print_iter_result(String path, double[][] mean,String[] probName,String[] algName,boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path,append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);
		      for(int i=0;i<mean.length;i++){
		    	  for(int j=0;j<mean[i].length;j++){
		    		  bw.write("  "+mean[i][j]);
		    	  }
		    	  bw.newLine();
		      }
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }
	  public static void printResult_pos(String path, double[][] mean,String[] probName,String[] algName,boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path,append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);
		      for(int i=0;i<mean.length;i++){
		    	  bw.write(probName[i]);
		    	  for(int j=0;j<mean[i].length;j++){
		    		  bw.write("  "+mean[i][j]);
		    	  }
		    	  bw.newLine();
		      }
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }

	  public static void printResult_init(String path, String[] algName,boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path,append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);
		      for(int i=0;i<algName.length;i++){
		    	  bw.write("      "+algName[i]+"      ");
		      }
		      bw.newLine();
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }
	  
	  public static void printResult_posinit(String path, String[] algName,boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path,append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);
		      for(int i=0;i<algName.length;i++){
		    	  bw.write("      "+algName[i]+"      ");
		      }
		      bw.newLine();
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }
	  
	  public static void printResult_igd(String path, double[] igd,String probName,String[] algName,boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path,append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);
		      bw.write(probName);
		      bw.newLine();
		      for(int i=0;i<igd.length;i++){
		    	
		    	  bw.write(igd[i]+"");
		    	  bw.newLine();
		      }
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }
	  
	  public static void printResultName(String path, double[][] mean,double[][] sd,String[] probName,String[] algName,boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path,append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);
		      for(int j=0;j<mean[0].length;j++){
		    	  bw.write("                           "+algName[j].toString());
		      }
		      bw.newLine();
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }

	  
	  public static void printVectorToFile(String path, Object[] vector, boolean append) {
		    try {
		      FileOutputStream fos = new FileOutputStream(path, append);
		      OutputStreamWriter osw = new OutputStreamWriter(fos);
		      BufferedWriter bw = new BufferedWriter(osw);

		      if (vector.length > 0)
		      {
		        for (int j = 0; j < vector.length; j++)
		          bw.write(vector[j].toString() + "  ");
		       bw.write("\n");
		        bw.newLine();
		      }
		      bw.close();
		    } catch (IOException e) {
		      e.printStackTrace();
		    }
		  }

	  public static Algorithm getAlgorithm(Problem p, String algName) {
		  
			return new MOEAD_ARA(p);
			
		}
	  
	  
	  public static Problem getProblem(String problemName) throws ClassNotFoundException {
			Problem p;
			if ("MOP1".equals(problemName)) {
				p = new MOP1("Real");
			} else if ("MOP2".equals(problemName)) {
				p = new MOP2("Real");
			} else if ("MOP3".equals(problemName)) {
				p = new MOP3("Real");
			} else if ("MOP4".equals(problemName)) {
				p = new MOP4("Real");
			} else if ("MOP5".equals(problemName)) {
				p = new MOP5("Real");
			} else if ("MOP6".equals(problemName)) {
				p = new MOP6("Real");
			} else if ("MOP7".equals(problemName)) {
				p = new MOP7("Real");
			}else if ("UF1".equals(problemName)) {
				p = new UF1("Real");
			} else if ("UF2".equals(problemName)) {
				p = new UF2("Real");
			} else if ("UF3".equals(problemName)) {
				p = new UF3("Real");
			} else if ("UF4".equals(problemName)) {
				p = new UF4("Real");
			} else if ("UF5".equals(problemName)) {
				p = new UF5("Real");
			} else if ("UF6".equals(problemName)) {
				p = new UF6("Real");
			} else if ("UF7".equals(problemName)) {
				p = new UF7("Real");
			} else if ("UF8".equals(problemName)) {
				p = new UF8("Real");
			} else if ("UF9".equals(problemName)) {
				p = new UF9("Real");
			} else if ("UF10".equals(problemName)) {
				p = new UF10("Real");
			}else {
				//p = new OKA2("Real");
				p=null;
			}
			return p;
		}

	  
}
