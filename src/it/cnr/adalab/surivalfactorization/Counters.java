package it.cnr.adalab.surivalfactorization;

public class Counters {
	double[] S_k;
	double[] A_k;
	SparseDoubleMatrix2D S_c_k;
	SparseDoubleMatrix2D A_c_k;
	SparseDoubleMatrix2D tilde_S_c_k;
	SparseDoubleMatrix2D tilde_A_c_k;
	SparseDoubleMatrix2D log_S_c_k;
	SparseDoubleMatrix2D log_A_c_k;

	
	
	SparseDoubleMatrix2D[] S_c_u_k;
	SparseDoubleMatrix2D[] A_c_u_k;
	SparseDoubleMatrix2D[] tilde_S_c_u_k;
	SparseDoubleMatrix2D[] tilde_A_c_u_k;
	int n_cascades;
	int n_users;
	int n_features;
	
	public Counters(int n_cascades,int n_users,int n_features){
		this.n_cascades = n_cascades;
		this.n_users = n_users;
		this.n_features = n_features;
		
		S_c_k = new SparseDoubleMatrix2D(n_cascades,n_features);
		A_c_k = new SparseDoubleMatrix2D(n_cascades,n_features);
		
		tilde_S_c_k = new SparseDoubleMatrix2D(n_cascades,n_features);
		tilde_S_c_k = new SparseDoubleMatrix2D(n_cascades,n_features);

		log_S_c_k = new SparseDoubleMatrix2D(n_cascades,n_features);
		log_A_c_k = new SparseDoubleMatrix2D(n_cascades,n_features);

		
		
		S_c_u_k = new SparseDoubleMatrix2D[n_cascades];
		A_c_u_k = new SparseDoubleMatrix2D[n_cascades];
		
		tilde_S_c_u_k = new SparseDoubleMatrix2D[n_cascades];
		tilde_A_c_u_k = new SparseDoubleMatrix2D[n_cascades];
		
		for (int c = 0; c < n_cascades; c++){
			S_c_u_k[c] = new SparseDoubleMatrix2D(n_users,n_features);
			A_c_u_k[c] = new SparseDoubleMatrix2D(n_users,n_features);
			tilde_S_c_u_k[c] = new SparseDoubleMatrix2D(n_users,n_features);
			tilde_A_c_u_k[c] = new SparseDoubleMatrix2D(n_users,n_features);			
		}
		
		
		
	}
	
	
	public void update(CascadeData data, Model model){
		
	}

}
