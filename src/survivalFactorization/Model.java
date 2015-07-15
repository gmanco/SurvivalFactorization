package survivalFactorization;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;

import utils.MatrixUtilities;
import utils.Randoms;

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
    
    HyperParameters hyperParams;
    
    private double[][] A; // n_nodes x n_features
    private double[][] S;  // n_nodes x n_features
    private double[][] Phi; // n_words x n_features


	/*
	 * This is the deep copy constructor
	 */
	public Model(Model m) {
		this.n_nodes = m.n_nodes;
		this.n_words = m.n_words;
		this.n_features = m.n_features;
		this.hyperParams = m.hyperParams;

        final int N = this.n_nodes;
        final int K = this.n_features;
        final int W = this.n_words;
        
        this.A = new double[N][K];
        MatrixUtilities.copy(m.A, this.A);
        
        this.S = new double[N][K];
        MatrixUtilities.copy(m.S, this.S);
        
        this.Phi = new double[W][K];
        MatrixUtilities.copy(m.Phi, this.Phi);
        
       this.hyperParams=new HyperParameters(m.hyperParams);

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
				A[n][k] = rng.nextGamma(a, 1/b);
				S[n][k] = rng.nextGamma(a, 1/b);
			}
			for (int w = 0; w < W; w++)
				Phi[w][k] = rng.nextGamma(C[w], 1/D[w]);
		}

	}//init

	

	public double computeLLk(CascadeData data, Counters counters) {
		double llk=0.0;
		int n_cascades=data.n_cascades;
		for (int c = 0; c < n_cascades; c++){
			List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);
			List<WordOccurrence> cascadeContent=data.getCascadeContent(c);
		    int n_events_cascade=cascadeEvents.size();
		
		    double F_c[]=counters.F_curr[c];
		      
		    
		    double first_component[]=new double[n_features];
		    for(int k=0;k<n_features;k++){
		        first_component[k]=-F_c[k]*(counters.S_k[k]-counters.S_c_k[c][k])*(data.t_max*counters.A_c_k[c][k]-counters.tilde_A_c_k[c][k]);
		    }
		    
		    double[]secondComponent=new double[n_features];
		    double[]thirdComponent=new double[n_features];
		   
		    //scan over events
		    for(int e=0;e<n_events_cascade;e++){
                CascadeEvent ce=cascadeEvents.get(e);
                int u=ce.node;
                double t_u=ce.timestamp;

                for(int k=0;k<n_features;k++){
                    secondComponent[k]+=S[u][k]*(t_u*counters.A_c_u_k[c].get(u, k)-counters.tilde_A_c_u_k[c].get(u,k));
                   
                    if(e>0){
                        double a=counters.A_c_u_k[c].get(u,k)-A[u][k];
                        if(a>0)
                            thirdComponent[k]+=Math.log(counters.A_c_u_k[c].get(u,k)-A[u][k]);
                        else
                            throw new RuntimeException(""+counters.A_c_u_k[c].get(u,k)+"\t"+A[u][k]);
                    }
                }//for each k
              
		    }// for each event
		    
		    
		    for(int k=0;k<n_features;k++){
		        secondComponent[k]=-F_c[k]*secondComponent[k];
		        if(F_c[k]>0){		        
		            thirdComponent[k]+=counters.log_S_c_k[c][k] + (n_events_cascade-1)*Math.log(F_c[k]);
		        }
		        else throw new RuntimeException();
		    }
		    
	   
		    //compute llk on the content
		    double content_llk[]=new double[n_features];
		    int cascadeContentSize=0;
		    for(WordOccurrence wo:cascadeContent){
		        int w=wo.word;
		        int n_w_c=wo.cnt;
		        cascadeContentSize+=n_w_c;
		        for(int k=0;k<n_features;k++){
		            if(Phi[w][k]>0)
		                content_llk[k]+=n_w_c*Math.log(Phi[w][k]);
		            else
		                throw new RuntimeException();
		        }//for each k
		    }// for each word
		    
		    for(int k=0;k<n_features;k++){
		        content_llk[k]-=(cascadeContentSize)*counters.Phi_k[k];
		    }//for each k
		    
		    //finally.. 
			for(int k=0;k<n_features;k++){
			    llk+=first_component[k]+secondComponent[k]+thirdComponent[k]+content_llk[k];
			}//for each k
			
			
		}//for each cascade
		
		return llk;
	}//computeLLk
	
	/*
	 * Compute F given a list of word occurrences
	 */
    public double[] computeF(List<WordOccurrence> W_c) {
        double[] F_c = new double[n_features];
        Arrays.fill(F_c,1.0);
       
        for(WordOccurrence wo:W_c){
            int w=wo.word;
            for (int k = 0; k < n_features; k++)
                F_c[k] *= Phi[w][k];
            for (int k = 0; k < n_features; k++)
                if (F_c[k] == 0)
                	F_c[k] = GibbsSampler.DEFAULT_SHAPE/GibbsSampler.DEFAULT_RATE_PHI;
            	
        }
        return F_c;
    }//computeF
	
    /*
     * Compute F for all the cascades
     */
	public double[][] computeFAllCascades(CascadeData data) {
	    int n_cascades=data.n_cascades;
	    double[][] F_curr = new double[n_cascades][n_features];
		for (int c = 0; c < n_cascades; c++) {
			F_curr[c] = computeF(data.getCascadeContent(c));
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


    public void setHyperParams(HyperParameters hyperParams) {
        this.hyperParams = hyperParams;
    }

    public void setA(double[][] a) {
        this.A = a;
    }
    

    public double[][] getA() {
        return A;
    }

    public double[][] getS() {
        return S;
    }

    public double[][] getPhi() {
        return Phi;
    }

   
    public void setS(double[][] s) {
        this.S = s;
    }

    public void setPhi(double[][] phi) {
        this.Phi = phi;
    }



}//Model
