package evaluation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.StringTokenizer;
import java.util.concurrent.locks.ReentrantReadWriteLock.ReadLock;

import utils.MatrixUtilities;

/**
 * @author Giuseppe Manco
 * 
 *         To change the template for this generated type comment go to Window -
 *         Preferences - Java - Code Generation - Code and Comments
 */
public class ClusteringEvaluator {

	
	double[][] ct;
	int nrows;
	int ncols;
	double n;
	double Z;
	double snr, snc;
	
	private ClusteringEvaluator(double[][] mat, boolean summaries) {
		if (mat != null) {
			n = 0;
			if (summaries) {
				nrows = mat.length - 1;
				ncols = mat[0].length - 1;
			} else {
				nrows = mat.length;
				ncols = mat[0].length;
			}
			ct = new double[nrows + 1][ncols + 1];

			for (int i = 0; i < nrows; i++) {
				for (int j = 0; j < ncols; j++) {
					ct[i][j] = mat[i][j];
					n += ct[i][j];
					Z += ct[i][j] * ct[i][j];
					ct[i][ncols] += ct[i][j];
					ct[nrows][j] += ct[i][j];
					ct[nrows][ncols] += ct[i][j];
				}
			}
			snr = 0;
			for (int i = 0; i < nrows; i++)
				snr += ct[i][ncols] * ct[i][ncols];

			snc = 0;
			for (int i = 0; i < ncols; i++)
				snc += ct[nrows][i] * ct[nrows][i];

		}
	}

	public double getRandMeasure() {
		double a, d, M;
		a = getA();
		d = getD();
		M = getM();

		return (a + d) / M;
	}
	
	public double getCorrectedRandMeasure() {
		 double num, den;
		 num = (Z/2.0 - n/2.0)- (2.0/(n*(n-1)))*(snc -n)/2.0*(snr-n)/2.0;
		 den =1.0/2.0*((snc -n)/2.0 + (snr-n)/2.0)- (2.0/(n*(n-1)))*(snc -n)/2.0*(snr-n)/2.0;
		 return num/den;
	}

	/**
	 * @return
	 */
	private double getM() {
		return n * (n - 1) / 2;
	}

	/**
	 * @return
	 */
	private double getA() {
		return 1.0 / 2 * Z - n / 2;
	}

	/**
	 * @return
	 */
	private double getB() {
		return 1.0 / 2 * (snc - Z);
	}

	private double getC() {
		return 1.0 / 2 * (snr - Z);
	}

	/**
	 * @return
	 */
	private double getD() {
		// return getM() -getA() -getB() -getC();
		return (n * n + Z - (snr + snc)) / 2;
	}

	/**
	 * @return
	 */
	public double getGammaMeasure() {
		double a = getA();
		double b = getB();
		double c = getC();
		// double d = getD();
		double M = getM();
		double m1 = a + b, m2 = a + c;
		double gamma = (M * a - m1 * m2)
				/ Math.sqrt(m1 * m2 * (M - m1) * (M - m2));
		return gamma;
	}

	/**
	 * @return
	 */
	public double getFowlMeasure() {
		double a = getA();
		double b = getB();
		double c = getC();

		return a / (Math.sqrt((a + c) * (a + b)));
	}

	/**
	 * @return
	 */
	public double getJaccMeasure() {
		double a = getA();
		double b = getB();
		double c = getC();
		return a / (a + b + c);
	}

	public double getFMeasure() {
		double f[] = new double[ncols];
		double prec, rec, find;

		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++) {
				prec = ct[i][j] / ct[i][ncols];
				rec = ct[i][j] / ct[nrows][j];
				if (prec != 0 && rec != 0)
					find = 2.0 * prec * rec / (prec + rec);
				else
					find = 0;
				if (find > f[j])
					f[j] = find;
			}
		find = 0;
		for (int j = 0; j < f.length; j++)
			find += f[j] * ct[nrows][j] / ct[nrows][ncols];
		return find;
	}
	
	public double getMutualInfo(){
		//calculate the MI (unadjusted)
		double MI=0.0;
		for (int i=1; i < nrows; i++){
		    for (int j=1; j < ncols; j++){
		        if (ct[i][j]>0) 
		        		MI=MI+(ct[i][j]/n)*Math.log((double)ct[i][j]*n/(ct[i][ncols]*ct[nrows][j]));
		    }
		}
		MI=MI/n;
		return MI;
	}

	
	public  double getNMI(){
		double nmi=0;
		double I=getMutualInfo();
		
		// Compute entropy
		double H_rows=0;
		for(int i=0;i<nrows;i++){
			double p_x_i=(double)ct[i][ncols]/n;
			H_rows+=-p_x_i*Math.log(p_x_i);	
		}
		
		double H_cols=0;
		for(int j=0;j<ncols;j++){
			double p_y_j=(double)ct[nrows][j]/n;
			H_cols+=-p_y_j*Math.log(p_y_j);	
		}
		nmi=I/((H_rows+H_cols)/2);
		return nmi;
	}
	
	/**
	 * @return
	 */
	public double getError() {
		int[] cm = new int[nrows];

		for (int i = 0; i < nrows; i++) {
			cm[i] = -1;
			double max = 0;
			for (int j = 0; j < ncols; j++)
				if (ct[i][j] > max) {
					cm[i] = j;
					max = ct[i][j];
				}
		}
		double err = 0;
		for (int i = 0; i < nrows; i++)
			for (int j = 0; j < ncols; j++)
				if (j != cm[i])
					err += ct[i][j];
		err /= ct[nrows][ncols];

		return err;
	}

	public double getWeigthedError() {
		double err = getError();
		if (nrows > ncols)
			err *= nrows * 1.0 / ncols;
		else
			err *= ncols * 1.0 / nrows;

		return err;
	}
	
	
	/*public double getMutualInfo(){
			double mi = 0, p_xy,p_x,p_y;
			for (int i = 0; i < nrows; i++) {
				for (int j = 0; j < ncols; j++) {
					p_xy = ct[i][j]/n;
					p_x = ct[i][ncols]/n;
					p_y = ct[ncols][j]/n;
					mi+= p_xy*Math.log(p_xy/(p_x*p_y));
				}
			}
			return mi;
	}
		 
	 public double getNMInfo(){
		    double nmi = getMutualInfo();
		    
		    double h_x = 0, h_y = 0, p_x,p_y;
		    
		    for (int i = 0; i < nrows; i++) {
		      p_x = ct[i][ncols]/n;
		      h_x += p_x*Math.log(p_x);
		    }
		    for (int j = 0; j < ncols; j++) {
		       p_y = ct[nrows][j]/n;
		        h_y += p_y*Math.log(p_y);
		    }
		    return (nmi/(h_x + h_y));
	 }*/
	
	public static void evaluate(int n_real_clusters,  HashMap<Integer, Integer> truth,HashMap<Integer, Integer> model){
		evaluate(n_real_clusters,n_real_clusters ,truth, model);
	}
	
	public static void evaluate(int n_real_clusters, int n_predicted_clusters, HashMap<Integer, Integer> truth,HashMap<Integer, Integer> model){
		double[][] a = new double[n_real_clusters][n_predicted_clusters];
		for(Integer element:truth.keySet()){
			Integer realclass=truth.get(element);
			Integer predictedClass=model.get(element);
			if(realclass==null || predictedClass==null)
				throw new RuntimeException("Exception on element\t"+element);
			a[realclass-1][predictedClass-1]++;
		}
		ClusteringEvaluator e = new ClusteringEvaluator(a, false);
		System.out.println(" - - - Clustering Evaluation - - - - ");
		System.out.println("F-Index:\t\t" + e.getFMeasure());
		System.out.println("Jaccard Statistics:\t" + e.getJaccMeasure());
		System.out.println("Rand Statistics:\t" + e.getRandMeasure());
		System.out.println("ARI Statistics:\t" + e.getCorrectedRandMeasure());
		System.out.println("Fowlkes Statistics:\t" + e.getFowlMeasure());
		System.out.println("Gamma Statistics:\t" + e.getGammaMeasure());
		System.out.println("Error:\t\t\t" + e.getError());
		System.out.println("Mutual:\t\t\t" + e.getMutualInfo());
		System.out.println("NMI:\t\t\t" + e.getNMI());

		System.out.println("Confusion matrix");
		MatrixUtilities.print(a);
		
		
		System.out.println();
		System.out.println(" - - - - - - - END - - - - - - - ");
	}
	

	public static void evaluate(ArrayList<Integer> truth, ArrayList<Integer> model){
		int real_offset = Collections.min(truth);
		int n_real_clusters = Collections.max(truth) - real_offset +1;
		int pred_offset = Collections.min(model);
		int n_predicted_clusters = Collections.max(model) - pred_offset +1;
		
		if (truth.size() != model.size())
			return;
		
		double[][] a = new double[n_real_clusters][n_predicted_clusters];
		for(int i= 0; i < truth.size(); i++)
			a[truth.get(i)-real_offset][model.get(i)-pred_offset]++;
		
		ClusteringEvaluator e = new ClusteringEvaluator(a, false);
		System.out.println(" - - - Clustering Evaluation - - - - ");
		System.out.println("F-Index:\t" + e.getFMeasure());
		System.out.println("Jaccard:\t" + e.getJaccMeasure());
		System.out.println("Rand:\t" + e.getRandMeasure());
		System.out.println("ARI:\t" + e.getCorrectedRandMeasure());
		System.out.println("Fowlkes:\t" + e.getFowlMeasure());
		System.out.println("Gamma:\t" + e.getGammaMeasure());
		System.out.println("Error:\t" + e.getError());
		System.out.println("MI:\t" + e.getMutualInfo());
		System.out.println("NMI:\t" + e.getNMI());

		System.out.println("Confusion matrix");
		MatrixUtilities.print(a);
		
		
		System.out.println();
		System.out.println(" - - - - - - - END - - - - - - - ");
	}
	
	
	public static void main(String[] args) throws Exception {
		System.out.println(" - - - - Clustering Evaluation - - - - ");
		if(args.length==0){
			printUsage();
			return;	
		}
		
		int nclusters=0;
		String model="";
		String groundTruth="";
		
		for (int i = 0; i < args.length; i++) {
			if (args[i].equalsIgnoreCase("--help")) {
				printUsage();
				return;
			}
			
			if (args[i].equals("-c")) {
				nclusters = Integer.parseInt(args[i + 1]);
				i++;
			}
			
			
			if (args[i].equals("-m")) {
				model = args[i + 1];
				i++;
			}
			
			
			if (args[i].equals("-g")) {
				groundTruth = args[i + 1];
				i++;
			}
			
		
			
		}// for each args
	
		
		BufferedReader in_ground = new BufferedReader(new FileReader(groundTruth));
		
		System.out.println("Reading ground truth from \t"+groundTruth);
		
		ArrayList<Integer> t = new ArrayList<Integer>();
		for (String line = in_ground.readLine(); line != null; line = in_ground.readLine()) {
			t.add(Integer.parseInt(line));
		}			
		in_ground.close();
		
		BufferedReader in_model = new BufferedReader(new FileReader(model));
		System.out.println("Reading model from \t"+model);

		ArrayList<Integer> p = new ArrayList<Integer>();
		for (String line = in_model.readLine(); line != null; line = in_model.readLine()) {
			p.add(Integer.parseInt(line));
		}			
		in_model.close();

		evaluate(t,p);
	}//main
	

	private static void printUsage() {
		System.out.println("-c <nClusters> -m <clusterModel> -g <groundTruth> ");
	}//printUsage
	
	
}
