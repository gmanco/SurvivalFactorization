package it.cnr.adalab.surivalfactorization;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

public class Model {
	int n_users;
	int n_words;
	int n_features;
	int n_cascades;
	HyperParameters hyperParams;
	double[][] A;
	double[][] S;
	double[][] Phi;
	double[][] F;

	public Model() {
		this.n_users = -1;
		this.n_words = -1;
		this.n_features = -1;
	}

	public Model(Model m) {
		this.n_users = m.n_users;
		this.n_words = m.n_words;
		this.n_features = m.n_features;
		this.hyperParams = m.hyperParams;
		this.A = m.A;
		this.S = m.S;
		this.Phi = m.Phi;
		this.F = m.F;
	}

	public Model(int n_users, int n_words, int n_features, HyperParameters h) {
		this.n_users = n_users;
		this.n_words = n_words;
		this.n_features = n_features;
		this.hyperParams = h;
	}

	public void init(Randoms rng) {
		double a = hyperParams.getA();
		double b = hyperParams.getB();
		double[] C = hyperParams.getC();
		double[] D = hyperParams.getD();

		final int N = this.n_users;
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

	public double[] computeF(Set<Integer> W_c) {
		Iterator<Integer> it = W_c.iterator();
		double[] F = new double[n_features];
		for (int i = 0; i < n_features; i++)
			F[i] = 1;
		while (it.hasNext()) {
			int w = it.next();
			for (int k = 0; k < n_features; k++)
				F[k] *= Phi[w][k];
		}
		return F;
	}

	public double computeLLk(CascadeData data) {
		double llk = 0;

		for (int c = 1; c < n_cascades; c++){
			double[] curr_llk = new double[n_features];
			
			double[][] cascade = data.getCascadeEvents(c);
			
			for (int k = 1; k < n_features; k++){
				
				
			}
		}
		
		return llk;
	}

	public double[][] computeFAllCascades(CascadeData data) {
		double[][] F_curr = new double[n_cascades][n_features];
		for (int c = 1; c < n_cascades; c++) {
			HashMap<Integer, Double> W_c = data.getContent(c);
			F_curr[c] = computeF(W_c.keySet());
		}
		return F_curr;
	}

}
