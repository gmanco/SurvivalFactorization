package it.cnr.adalab.surivalfactorization;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

public class GibbsSamplerState {

	protected int[] M_k;
	protected int[] M_v;
	protected int[] Z;
	protected SparseDoubleMatrix2D Y;

	int n_users;
	int n_words;
	int n_features;
	int n_cascades;
	SparseDoubleMatrix2D[] N_k_u_v; // n_users x n_users x n_features
	SparseDoubleMatrix2D[] L_k_u_v; // n_users x n_users x n_features
	int[][] N_k_w; // n_words x n_features
	int[][] C_k_w; // n_words x n_features

	public GibbsSamplerState(int n_users, int n_cascades, int n_words,
			int n_features) {
		this.n_users = n_users;
		this.n_words = n_words;
		this.n_features = n_features;
		this.n_cascades = n_cascades;
		this.N_k_u_v = new SparseDoubleMatrix2D[n_features];
		this.L_k_u_v = new SparseDoubleMatrix2D[n_features];
		for (int c = 1; c < n_features; c++) {
			this.N_k_u_v[c] = new SparseDoubleMatrix2D(n_users, n_users);
			this.L_k_u_v[c] = new SparseDoubleMatrix2D(n_users, n_users);
		}
		this.N_k_w = new int[n_words][n_features];
		this.C_k_w = new int[n_words][n_features];
		this.M_k = new int[n_features];
		this.M_v = new int[n_features];
		this.Z = new int[n_cascades];
		this.Y = new SparseDoubleMatrix2D(n_cascades, n_users);
	}

	protected void update(CascadeData data, int[] z_new,
			SparseDoubleMatrix2D y_new) {
		resetCounters();
		this.Z = z_new;
		this.Y = y_new;

		for (int c = 1; c < n_cascades; c++) {
			int k = Z[c];
			M_k[k]++;

			// process events

			// current cascade (activation times)
			double[][] T_c = data.getCascadeEvents(c);

			// number of events in this cascade
			int n_events_cascade = T_c.length;

			// loop on all the events, starting from the second
			for (int e = 2; e < n_events_cascade; e++) {
				// current event
				int u = (int)T_c[e][1];
				double t_u = T_c[e][2];

				// id influencer
				int v = (int) (Y.get(c, u));

				// time of activation of the influencer
				double t_v = data.getCascadeContent().getColumn(c).get(v);

				// time gap between activations
				double delta_uv = t_u - t_v;

				// update counter
				N_k_u_v[k].set(u, v, N_k_u_v[k].get(u, v) + 1);
				L_k_u_v[k].set(u, v, L_k_u_v[k].get(u, v) + Math.log(delta_uv));
				M_v[v]++;

			}

			// loop on the content of the cascade

			HashMap<Integer, Double> W_c = data.getContent(c);

			// process content
			Iterator<Entry<Integer, Double>> it = W_c.entrySet().iterator();
			while (it.hasNext()) {
				// current word and count occurrence
				Entry<Integer, Double> word = it.next();
				int w = word.getKey();
				double n_w = word.getValue();
				N_k_w[w][k] += n_events_cascade - 1;
				C_k_w[w][k] += n_w;
			}

		}
	}

	public void randomInitZ(double[] p) {
		Multinomial m = new Multinomial(p);
		for (int c = 1; c < n_cascades; c++)
			Z[c] = m.sample();
	}

	public void resetCounters() {
		this.N_k_u_v = new SparseDoubleMatrix2D[n_features];
		this.L_k_u_v = new SparseDoubleMatrix2D[n_features];
		for (int c = 1; c < n_features; c++) {
			this.N_k_u_v[c] = new SparseDoubleMatrix2D(n_users, n_users);
			this.L_k_u_v[c] = new SparseDoubleMatrix2D(n_users, n_users);
		}
		this.N_k_w = new int[n_words][n_features];
		this.C_k_w = new int[n_words][n_features];
		this.M_k = new int[n_features];
		this.M_v = new int[n_features];
	}

}
