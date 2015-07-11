package it.cnr.adalab.surivalfactorization;

import java.io.Serializable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;

import data.CascadeData;
import data.WordOccurrence;

public class Model implements Serializable{
	
    int n_nodes;
    int n_words;
    int n_features;
    int n_cascades;
    
    HyperParameters hyperParams;
    
    private double[][] A;
    private double[][] S;
    private double[][] Phi;
    private double[][] F;

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

		for (int c = 1; c < n_cascades; c++){
			double[] curr_llk = new double[n_features];
			
			//double[][] cascade = data.getCascadeEvents(c);
			
			for (int k = 1; k < n_features; k++){
				
				
			}
		}
		
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
		for (int c = 1; c < n_cascades; c++) {
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
