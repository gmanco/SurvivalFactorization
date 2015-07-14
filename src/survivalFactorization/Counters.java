package survivalFactorization;

import java.util.List;

import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import data.CascadeData;
import data.CascadeEvent;


public class Counters {
	
    // factor-wise counters
    public double[] S_k;
    public double[] A_k;
    public double[] Phi_k;

	// cascade based counter
    public double[][] S_c_k;
    public double[][] A_c_k;
    public double[][] tilde_S_c_k;
    public double[][] tilde_A_c_k;
	//double[][] log_S_c_k;
	//double[][] log_A_c_k;

	
	//observation based counters
    public SparseDoubleMatrix2D[] S_c_u_k;
    public SparseDoubleMatrix2D[] A_c_u_k;
    public SparseDoubleMatrix2D[] tilde_S_c_u_k;
    public SparseDoubleMatrix2D[] tilde_A_c_u_k;
	
	//properties of the data
	int n_cascades;
	int n_nodes;
	int n_features;
	
	public Counters(int n_cascades,int n_nodes,int n_features){
		this.n_cascades = n_cascades;
		this.n_nodes = n_nodes;
		this.n_features = n_features;
		init();
	}
	
	private void init(){
	    S_k=new double[n_features];
	    A_k=new double[n_features];
	    Phi_k=new double[n_features];
	    
	    S_c_k = new double[n_cascades][n_features];
        A_c_k = new double[n_cascades][n_features];
        
        tilde_S_c_k = new double[n_cascades][n_features];
        tilde_A_c_k = new double[n_cascades][n_features];

       // log_S_c_k = new double[n_cascades][n_features];
       // log_A_c_k = new double[n_cascades][n_features];

        
        S_c_u_k = new SparseDoubleMatrix2D[n_cascades];
        A_c_u_k = new SparseDoubleMatrix2D[n_cascades];
        
        tilde_S_c_u_k = new SparseDoubleMatrix2D[n_cascades];
        tilde_A_c_u_k = new SparseDoubleMatrix2D[n_cascades];
        
        for (int c = 0; c < n_cascades; c++){
            S_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);
            A_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);
            tilde_S_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);
            tilde_A_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);            
        }//for each cascade
	}//init
	
	/*
	 * Update counters 
	 */
	public void update(CascadeData data, Model model){
		//reset counters
	    init();
		
	    double A[][]=model.getA();
	    double S[][]=model.getS();
	    double Phi[][]=model.getPhi();
	    
	    int n_nodes=model.n_nodes;
	    int n_cascades=data.n_cascades;
	    int n_features=model.n_features;
	    int n_words=model.n_words;
	    
	    // compute factor wise counters
	    for(int k=0;k<n_features;k++){
    	    for(int u=0;u<n_nodes;u++){
    	            A_k[k]+=A[u][k];
    	            S_k[k]+=S[u][k];
    	    }
    	    for(int w=0;w<n_words;w++)
    	        Phi_k[k]+=Phi[w][k];        
	    }// for each k
	    
	    // compute cascade and event based counter
	    for (int c = 0; c < n_cascades; c++){
            List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);
            int n_events_cascade=cascadeEvents.size();
                  
            // data structure for recursive definitions
            double cumulative_A_v_k[]=new double[n_features];
            double cumulative_tilde_A_v_k[]=new double[n_features];
            double cumulative_S_v_k[]=new double[n_features];
            double cumulative_tilde_S_v_k[]=new double[n_features];
            
            // loop over events
            for(int e=0;e<n_events_cascade;e++){
                CascadeEvent ce=cascadeEvents.get(e);
                int u=ce.node;
                double t=ce.timestamp;
				for (int k = 0; k < n_features; k++) {
					S_c_k[c][k] += S[u][k];
					// log_S_c_k[c][k]+=Math.log(S[u][k]);
					tilde_S_c_k[c][k] += S[u][k] * t;
					// log_A_c_k[c][k]+=Math.log(A[u][k]);
					A_c_k[c][k] += A[u][k];
					tilde_A_c_k[c][k] += A[u][k] * t;

					cumulative_S_v_k[k] += S[u][k];
					cumulative_tilde_S_v_k[k] += S[u][k] * t;
					cumulative_A_v_k[k] += A[u][k];
					cumulative_tilde_A_v_k[k] += A[u][k] * t;

					S_c_u_k[c].set(u, k, cumulative_S_v_k[k]);
					tilde_S_c_u_k[c].set(u, k, cumulative_tilde_S_v_k[k]);
					A_c_u_k[c].set(u, k, cumulative_A_v_k[k]);
					tilde_A_c_u_k[c].set(u, k, cumulative_tilde_A_v_k[k]);

				}// for each k
			}//
	    
	    }//for each cascade
	    
	    
	}//update

	/*
	 * Update counters 
	 */
	public void updateA(CascadeData data, Model model){
		//reset counters
        A_k=new double[n_features];        
        A_c_k = new double[n_cascades][n_features];
        tilde_A_c_k = new double[n_cascades][n_features];
        A_c_u_k = new SparseDoubleMatrix2D[n_cascades];
        tilde_A_c_u_k = new SparseDoubleMatrix2D[n_cascades];        
        for (int c = 0; c < n_cascades; c++){
            A_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);
            tilde_A_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);            
        }//for each cascade
        
		
	    double A[][]=model.getA();
	    
	    int n_nodes=model.n_nodes;
	    int n_cascades=data.n_cascades;
	    int n_features=model.n_features;
	    
	    // compute factor wise counters
	    for(int k=0;k<n_features;k++)
    	    for(int u=0;u<n_nodes;u++){
    	            A_k[k]+=A[u][k];
    	    }
	    
	    // compute cascade and event based counter
	    for (int c = 0; c < n_cascades; c++){
            List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);
            int n_events_cascade=cascadeEvents.size();
                  
            // data structure for recursive definitions
            double cumulative_A_v_k[]=new double[n_features];
            double cumulative_tilde_A_v_k[]=new double[n_features];
            
            // loop over events
            for(int e=0;e<n_events_cascade;e++){
                CascadeEvent ce=cascadeEvents.get(e);
                int u=ce.node;
                double t=ce.timestamp;
				for (int k = 0; k < n_features; k++) {
					// log_A_c_k[c][k]+=Math.log(A[u][k]);
					A_c_k[c][k] += A[u][k];
					tilde_A_c_k[c][k] += A[u][k] * t;

					cumulative_A_v_k[k] += A[u][k];
					cumulative_tilde_A_v_k[k] += A[u][k] * t;

					A_c_u_k[c].set(u, k, cumulative_A_v_k[k]);
					tilde_A_c_u_k[c].set(u, k, cumulative_tilde_A_v_k[k]);

				}// for each k
			}//
	    
	    }//for each cascade
	    
	    
	}//updateA
	
	/*
	 * Update counters 
	 */
	public void updateS(CascadeData data, Model model){
	    S_k=new double[n_features];
        S_c_k = new double[n_cascades][n_features];        
        tilde_S_c_k = new double[n_cascades][n_features];
        S_c_u_k = new SparseDoubleMatrix2D[n_cascades];        
        tilde_S_c_u_k = new SparseDoubleMatrix2D[n_cascades];
        for (int c = 0; c < n_cascades; c++){
            S_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);
            tilde_S_c_u_k[c] = new SparseDoubleMatrix2D(n_nodes,n_features);
        }//for each cascade
		
	    double S[][]=model.getS();
	    
	    int n_nodes=model.n_nodes;
	    int n_cascades=data.n_cascades;
	    int n_features=model.n_features;
	    
	    // compute factor wise counters
	    for(int k=0;k<n_features;k++)
    	    for(int u=0;u<n_nodes;u++){
    	            S_k[k]+=S[u][k];
    	    }
	    
	    // compute cascade and event based counter
	    for (int c = 0; c < n_cascades; c++){
            List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);
            int n_events_cascade=cascadeEvents.size();
                  
            // data structure for recursive definitions
            double cumulative_S_v_k[]=new double[n_features];
            double cumulative_tilde_S_v_k[]=new double[n_features];
            
            // loop over events
            for(int e=0;e<n_events_cascade;e++){
                CascadeEvent ce=cascadeEvents.get(e);
                int u=ce.node;
                double t=ce.timestamp;
				for (int k = 0; k < n_features; k++) {
					S_c_k[c][k] += S[u][k];
					// log_S_c_k[c][k]+=Math.log(S[u][k]);
					tilde_S_c_k[c][k] += S[u][k] * t;

					cumulative_S_v_k[k] += S[u][k];
					cumulative_tilde_S_v_k[k] += S[u][k] * t;

					S_c_u_k[c].set(u, k, cumulative_S_v_k[k]);
					tilde_S_c_u_k[c].set(u, k, cumulative_tilde_S_v_k[k]);
				}// for each k
			}//
	    
	    }//for each cascade
	    
	    
	}//update
}//Counters
