package survivalFactorization;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.List;

import utils.Multinomial;

import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import cern.colt.matrix.tint.impl.SparseIntMatrix2D;

import data.CascadeData;
import data.CascadeEvent;
import data.WordOccurrence;

public class GibbsSamplerState {

    static double MIN_DELTA=1.0;
    
	protected int[] M_k;
	protected int[] M_v;
	protected int[] Z;
	public double [] Gamma_k;
	public double [] tilde_Gamma_k;
	
	protected SparseIntMatrix2D Y; //n_cascades x n_users

	int n_nodes;
	int n_words;
	int n_features;
	int n_cascades;
	
	/*
	 * N_k_u_v[u](u,v) = \sum_{c: u \in V^+(c), v \prec_c u} z_{c,k}* y^c_{u,v}
	 */
	SparseDoubleMatrix2D[] N_k_u_v; // n_features x n_users x n_users x 
	
	/*
	 * \sum_v n^k_{u,v}
	 */
	double[][] n_k_u_pre; //  n_users x n_features
	
	/*
     * \sum_v n^k_{v,u}
     */
	double [][] n_k_u_post; //  n_users x n_features 
	
	SparseDoubleMatrix2D[] L_k_u_v; // n_features x n_users x n_users
	int[][] N_k_w; // n_words x n_features
	int[][] C_k_w; // n_words x n_features

	public GibbsSamplerState(int n_users, int n_cascades, int n_words,
			int n_features) {
		this.n_nodes = n_users;
		this.n_words = n_words;
		this.n_features = n_features;
		this.n_cascades = n_cascades;
		
		/*
		 * Initialize counter
		 */
		resetCounters();
		
		/*
		 * Initialize structures for latent factors
		 */
		this.Z = new int[n_cascades];
		this.Y = new SparseIntMatrix2D(n_cascades, n_users);
	}

	/*
	 * Given a new status of latent variable,
	 * this method updates the counters of the gibbs sampler state
	 */
	protected void update(CascadeData data, int[] z_new,
	        SparseIntMatrix2D y_new,double[][]F_curr,Counters counters) {
		
	    // first of all: reset counters
	    resetCounters();
	    
	    //update the status of the latent variable
		this.Z = z_new;
		this.Y = y_new;

		//loop on cascade
		for (int c = 0; c < n_cascades; c++) {
		    
		    // factor associated with cascade c
			int k_c = Z[c];
			
			M_k[k_c]++;
			
			// current cascade (activation times)
			List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);

			// number of events in this cascade
			int n_events_cascade = cascadeEvents.size();

			// loop on all the events, starting from the second
			if(n_events_cascade>1)
    			for (int e = 2; e < n_events_cascade; e++) {
    				// current event
    			    CascadeEvent curr=cascadeEvents.get(e);
    			    
    				int u = curr.node;
    				double t_u = curr.timestamp;
    				
    				// id influencer
    				int v = (int) (Y.get(c, u));
    
    				// time of activation of the influencer
    				double t_v = data.getActivationTimestamp(v,c);
    
    				// time gap between activations
    				double delta_uv = t_u - t_v;
    				if(delta_uv==0){
    				    //System.err.println("Delta is zero");
    				    delta_uv=MIN_DELTA;
    				}
    				// update counters
    				N_k_u_v[k_c].set(u, v, N_k_u_v[k_c].get(u, v) + 1);
    				L_k_u_v[k_c].set(u, v, L_k_u_v[k_c].get(u, v) + Math.log(delta_uv));
    				
    				//update pre and post
    				n_k_u_pre[u][k_c]++;
    				n_k_u_post[v][k_c]++;
    				

    				M_v[v]++;
    			}
			
			// process the content of the cascade
			List<WordOccurrence> cascadeContent = data.getCascadeContent(c);

			for(WordOccurrence wo:cascadeContent){
			    int w=wo.word;
			    int n_w=wo.cnt;
			    N_k_w[w][k_c] += n_events_cascade - 1;
                C_k_w[w][k_c] += n_w;
			}//for each word in the cascade

			
			//update Gamma_k and tilde_Gamma_k
			Gamma_k[k_c]+=F_curr[c][k_c]*counters.A_c_k[c][k_c];
			tilde_Gamma_k[k_c]+=F_curr[c][k_c]*counters.tilde_A_c_k[c][k_c];
			
		}//for each cascade
		
	}//update
	
	protected void updateGamma(CascadeData data, double[][] F_curr,
			Counters counters) {
	    Gamma_k=new double[n_features];
	    tilde_Gamma_k=new double[n_features];
	    // loop on cascade
		for (int c = 0; c < n_cascades; c++) {

			// factor associated with cascade c
			int k_c = Z[c];

			// update Gamma_k and tilde_Gamma_k
			Gamma_k[k_c] += F_curr[c][k_c] * counters.A_c_k[c][k_c];
			tilde_Gamma_k[k_c] += F_curr[c][k_c] * counters.tilde_A_c_k[c][k_c];

		}// for each cascade

	}// updateGamma

	public void randomInitZ(double[] p) {
		Multinomial m = new Multinomial(p);
		for (int c = 1; c < n_cascades; c++)
			Z[c] = m.sample();
	}//randomInitZ

	/*
	 * Allocate new counters
	 */
	private void resetCounters() {
		this.N_k_u_v = new SparseDoubleMatrix2D[n_features];
		this.L_k_u_v = new SparseDoubleMatrix2D[n_features];
		this.n_k_u_pre=new double[n_nodes][n_features];
        this.n_k_u_post=new double[n_nodes][n_features];
		for (int k = 0; k < n_features; k++) {
			this.N_k_u_v[k] = new SparseDoubleMatrix2D(n_nodes, n_nodes);
			this.L_k_u_v[k] = new SparseDoubleMatrix2D(n_nodes, n_nodes);
		}
		this.N_k_w = new int[n_words][n_features];
		this.C_k_w = new int[n_words][n_features];
		this.M_k = new int[n_features];
		this.M_v = new int[n_nodes];
		this.Gamma_k=new double[n_features];
		this.tilde_Gamma_k=new double[n_features];
	}//resetCounters

	/*
	 * Prints a summary of the status of hidden variables on a file
	 */
	public void printSummaryStatus(CascadeData data,String file) throws Exception{
	    PrintWriter pw=new PrintWriter(new FileWriter(file));
	    pw.println(" SUMMARY OF THE CURRENT STATE");
	    pw.println();
	    printZ(Z, data,pw);
	    printMk(pw);
	    printY(Y, data,pw);
	    printMv(pw);
	    pw.println();
	    pw.println("===========================================================");
	    pw.flush();
	    pw.close();
	}
	 
    private void printZ(int Z[],CascadeData data,PrintWriter pw){
        pw.println("****** Z *********");
        int n_cascades=data.n_cascades;
        pw.println("CascadeId \t Z");
        for(int c=0;c<n_cascades;c++){
           pw.println(""+c+"\t"+Z[c]);
        }//for each c   
        pw.println("********************");
    }//printZ
    
    
    private void printY(SparseIntMatrix2D Y, CascadeData data,PrintWriter pw){
        
        int n_cascades=data.n_cascades;
        pw.println("****************** Y ******************************************");
        pw.println("Cascade \t Node \t Timestamp \t NodeInfluencer \t TimestampInfluencer \t Delta");

        for(int c=0;c<n_cascades;c++){
            List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);
            int n_events_cascade=cascadeEvents.size();
        
            for(int e=0;e<n_events_cascade;e++){
                CascadeEvent ce=cascadeEvents.get(e);
                int u=ce.node;
                double t_u=ce.timestamp;
                if(e>0){
                    int v=Y.get(c,u);
                    double t_v=data.getActivationTimestamp(v, c);
                    double delta=t_u-t_v;
                    pw.format("%d \t %d \t %.3f \t %d \t %.3f \t %.3f %n", c,u,t_u,v,t_v,delta);
                }
            
            }// for each event
        }//for each c
        pw.println("****************************************************************");
        pw.println();
    }//printY
    
    private void printMv(PrintWriter pw){
        pw.println();
        pw.println("***********************");
        pw.println("NodeId\t CntInfluencer");
        for(int u=0;u<M_v.length;u++)
            pw.println(""+u+"\t"+M_v[u]);
        pw.println("***********************");
    }
    
    
    private void printMk(PrintWriter pw){
        pw.println();
        pw.println("***********************");
        pw.println("Feature\t Cnt");
        for(int k=0;k<M_k.length;k++)
            pw.println(""+k+"\t"+M_k[k]);
        pw.println("***********************");
    }
    
    
    
    
}//GibbsSamplerState
