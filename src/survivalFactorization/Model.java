package survivalFactorization;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import data.CascadeData;
import data.CascadeEvent;
import data.WordOccurrence;

public class Model implements Serializable{
	
    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    int n_nodes;
    int n_words;
    int n_features;
    int n_cascades;
    
    HyperParameters hyperParams;
    
    private double[][] A; // n_nodes x n_features
    private double[][] S;  // n_nodes x n_features
    private double[][] Phi; // n_words x n_features
    private double[][] F; // n_cascade x n_features

	public Model() {
		this.n_nodes = -1;
		this.n_words = -1;
		this.n_features = -1;
	}

	public Model(Model m) {
		this.n_nodes = m.n_nodes;
		this.n_words = m.n_words;
		this.n_features = m.n_features;
		this.hyperParams = m.hyperParams;
		this.A = m.A;
		this.S = m.S;
		this.Phi = m.Phi;
		this.F = m.F;
	}

	public Model(int n_users, int n_words, int n_features, HyperParameters h) {
		this.n_nodes = n_users;
		this.n_words = n_words;
		this.n_features = n_features;
		this.hyperParams = h;
	}

	public void init(Randoms rng) {
		double a = hyperParams.getA();
		double b = hyperParams.getB();
		double[] C = hyperParams.getC();
		double[] D = hyperParams.getD();

		final int N = this.n_nodes;
		final int K = this.n_features;
		final int W = this.n_words;

		A = new double[N][K];
		S = new double[N][K];
		Phi = new double[W][K];

		for (int k = 0; k < K; k++) {
			for (int n = 0; n < N; n++) {
				A[n][k] = rng.nextGamma(a, b);
				S[n][k] = rng.nextGamma(a, b);
			}
			for (int w = 0; w < W; w++)
				Phi[w][k] = rng.nextGamma(C[w], D[w]);
		}

	}

	

	public double computeLLk(CascadeData data) {
		double llk = 0;
		//compute S_k[k]=\sum_u S[u][k]
		double S_k[]=new double[n_features];
		for(int u=0;u<n_nodes;u++){
		    for(int k=0;k<n_features;k++){
		        S_k[k]+=S[u][k];
		    }
		}//for each node
		
		//compute Phi_k[k]=\sum_w Phi[w][k]
		double Phi_k[]=new double[n_features];
		for(int w=0;w<n_words;w++){
		    for(int k=0;k<n_features;k++){
		        Phi_k[k]+=Phi[w][k];
		    }
		}// for each word
		    
		
		
		
		for (int c = 0; c < n_cascades; c++){
			List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);
			List<WordOccurrence> cascadeContent=data.getCascadeContent(c);
		    int n_events_cascade=cascadeEvents.size();
		
		    double F_c[]=computeF(cascadeContent);
		    
			   
		    double A_c_k[]=new double[n_features];
		    double tilde_A_c_k[]=new double[n_features];
		    double S_c_k[]=new double[n_features];
		    
		    //first scan: compute counters
		    for(int e=0;e<n_events_cascade;e++){
		        CascadeEvent ce=cascadeEvents.get(e);
		        int u=ce.node;
		        double t=ce.timestamp;
		        for(int k=0;k<n_features;k++){
	                if(e>0){ //skip first activation
    		            S_c_k[k]+=S[u][k];
	                }
	                A_c_k[k]+=A[u][k];
	                tilde_A_c_k[k]+=A[u][k]*t;
	            }
		    }//first scan
		   
		    
		    double first_component[]=new double[n_features];
		    for(int k=0;k<n_features;k++){
		        first_component[k]=-F_c[k]*(S_k[k]-S_c_k[k])*(data.t_max*A_c_k[k]-tilde_A_c_k[k]);
		    }
		    
		    double secondComponent[]=new double[n_features];
		   
	        double sum_log_S_active[]=new double[n_features];
            double sum_log_cumulative_A[]=new double[n_features];
		    double cumulative_A_v_k[]=new double[n_features];
		    double cumulative_tilde_A_v_k[]=new double[n_features];
		    //second scan
		    for(int e=0;e<n_events_cascade;e++){
                CascadeEvent ce=cascadeEvents.get(e);
                int u=ce.node;
                double t_u=ce.timestamp;
		    
                //note the order of these updates: it matters!
                for(int k=0;k<n_features;k++){
                    secondComponent[k]+=S[u][k]*(t_u*cumulative_A_v_k[k]-cumulative_tilde_A_v_k[k]); 
                    
                    if(Math.log(S[u][k])>0)
                        sum_log_S_active[k]+=Math.log(S[u][k]); // this goes into the third component
                    else
                        throw new RuntimeException();
                    
                    if(e>1) { //this goes into the third component
                        if(cumulative_A_v_k[k]>0)
                            sum_log_cumulative_A[k]+=Math.log(cumulative_A_v_k[k]);
                        else
                            throw new RuntimeException();
                    }
                    
                    //update cumulatives
                    cumulative_A_v_k[k]+=A[u][k];
                    cumulative_tilde_A_v_k[k]+=A[u][k]*t_u;
                  
                }//for each k
              
		    }// for each event
		    
		    //multiply for F_c_k
		    for(int k=0;k<n_features;k++){
		        secondComponent[k]=-F_c[k]*secondComponent[k];
		    }
		    
		    double third_component[]=new double[n_features];
		    for(int k=0;k<n_features;k++){
		        third_component[k]+=sum_log_cumulative_A[k]+sum_log_S_active[k]+(n_events_cascade-1)*F_c[k];
		    }// for each k
		    
		    //compute llk on the content
		    double content_llk[]=new double[n_features];
		    for(WordOccurrence wo:cascadeContent){
		        int w=wo.word;
		        int n_w_c=wo.cnt;
		        for(int k=0;k<n_features;k++){
		            if(Phi[w][k]>0)
		                content_llk[k]+=n_w_c*Math.log(Phi[w][k]);
		            else
		                throw new RuntimeException();
		        }//for each k
		    }// for each word
		    
		    for(int k=0;k<n_features;k++){
		        content_llk[k]-=(cascadeContent.size())*Phi_k[k];
		    }//for each k
		    
		    //finally.. 
			for(int k=0;k<n_features;k++){
			    llk+=first_component[k]+secondComponent[k]+third_component[k]+content_llk[k];
			}//for each k
			
			
		}//for each cascade
		
		return llk;
	}//computeLLk
	
	
	

    public double[] computeF(List<WordOccurrence> W_c) {
       
        double[] F = new double[n_features];
        for (int k = 0; k < n_features; k++)
            F[k] = 1;
       
        for(WordOccurrence wo:W_c){
            int w=wo.word;
            for (int k = 0; k < n_features; k++)
                F[k] *= Phi[w][k];
        }
        return F;
    }//computeF
	

	public double[][] computeFAllCascades(CascadeData data) {
		double[][] F_curr = new double[n_cascades][n_features];
		for (int c = 0; c < n_cascades; c++) {
			List<WordOccurrence>W_c = data.getCascadeContent(c);
			F_curr[c] = computeF(W_c);
		}
		return F_curr;
	}//computeFAllCascades

    public void setN_nodes(int n_nodes) {
        this.n_nodes = n_nodes;
    }

    public void setN_words(int n_words) {
        this.n_words = n_words;
    }

    public void setN_features(int n_features) {
        this.n_features = n_features;
    }

    public void setN_cascades(int n_cascades) {
        this.n_cascades = n_cascades;
    }

    public void setHyperParams(HyperParameters hyperParams) {
        this.hyperParams = hyperParams;
    }

    public void setA(double[][] a) {
        A = a;
    }

    public void setS(double[][] s) {
        S = s;
    }

    public void setPhi(double[][] phi) {
        Phi = phi;
    }

    public void setF(double[][] f) {
        F = f;
    }

}//Model
