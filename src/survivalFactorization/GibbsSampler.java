package survivalFactorization;

import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;
import data.CascadeData;

public class GibbsSampler {

	double a;
	double b;
	double[] C; // n_words
	double[] D; // n_words
	double[] Alpha; // n_features
	double[] Beta; // n_users
	SamplerSettings settings;

	public GibbsSampler(SamplerSettings settings) {
		this.settings = settings;
	}

	public Model[] runInference(CascadeData data, int n_features) {

		int n_nodes = data.getNNodes();
		int n_words = data.getNWords();
		int n_cascades = data.getNCascades();
		int n_iterations = settings.n_iterations;
		int burnin = settings.burnin;

		// seed the seed
		int seed = settings.seed;
		Randoms rnd = new Randoms(seed);

		Model[] models = new Model[n_iterations - burnin];

		inferHyperParams(data, n_features);

		// init paramers
		Model model = new Model(n_nodes, n_words, n_features, new HyperParameters(a,b,C,D));
		model.init(rnd);

		long tot_time = 0;

		GibbsSamplerState curr_state = new GibbsSamplerState(n_nodes,
				n_cascades, n_words, n_features);

		double[] p = (new Dirichlet(Alpha)).nextDistribution();

		curr_state.randomInitZ(p);
		
		Counters counters=new Counters(n_cascades, n_nodes, n_features);
		counters.update(data, model);
		// Gibbs sampling
		//
		// Z is n_cascades x 1
		// Y is n_cascades x n_users
		//
		for (int epoch = 0; epoch < n_iterations; epoch++) {
			long t = System.currentTimeMillis();
			
			
			
			GibbsSamplerState next_state = sampleNextState(model, data,
					curr_state,counters);
			curr_state = next_state;

			double[][] F_curr = computeF(model, data);

			double[][] A_new = sampleA(model, data, curr_state, F_curr,counters);
			model.setA(A_new);

			double[][] S_new = sampleS(model, data, curr_state, F_curr,counters);
			model.setS(S_new);

			double[][] Phi_new = samplePhi(model, data, curr_state, F_curr,counters);
			model.setPhi(Phi_new);

			//now it's time to update counters
            counters.update(data, model);
            
			if (epoch > burnin) {
			    Model curr_model = new Model(model);
				models[epoch - burnin] = curr_model;
			}

			// update time
			long iteration_time = (System.currentTimeMillis() - t) / 1000;
			tot_time += iteration_time;

			// print out information about iteration
			if (epoch % settings.llk_interval == 0 && settings.compute_llk == 1) {
				double llk = model.computeLLk(data,counters);
				System.out
						.format("Iteration %d completed [elapsed time: %.0fs (llk: %-.2f)].\n",
								epoch, iteration_time, llk);
			} else
				System.out.format(
						"Iteration %d completed [elapsed time: %.0fs].\n",
						epoch, iteration_time);

		}
		// save(output,'Models');

		System.out
				.format("**************************************************\n");
		System.out.format(
				"*      SAMPLING COMPLETED [elapsed time: %.0fs]     *\n",
				tot_time);
		System.out
				.format("**************************************************\n");

		return models;
	}

	double[][] computeF(Model model, CascadeData data) {
		return model.computeFAllCascades(data);
	}

	private GibbsSamplerState sampleNextState(Model model, CascadeData data,
			GibbsSamplerState curr_state,Counters counters) {

		GibbsSamplerState next_state = new GibbsSamplerState(data.n_nodes,
				data.n_cascades, data.n_words, model.n_features);

		// Counters
		int[] M_k = curr_state.M_k;
		int[] M_v = curr_state.M_v;
		int[] Z = curr_state.Z;
		SparseDoubleMatrix2D Y = curr_state.Y;

		SparseDoubleMatrix2D Y_new = sampleY(model, data, Y, Z, M_v,counters);
		int[] Z_new = sampleZ(model, data, Y_new, Z, M_k,counters);

		next_state.update(data, Z_new, Y_new);

		return next_state;

	}//sampleNextState

	protected void inferHyperParams(CascadeData data, int n_features) {
		
		
		inferPriorRates(data);

		inferPriorWords(data, n_features);

		// these are default values used in topic models
		Alpha = new double[n_features];
		Beta = new double[data.n_nodes];
		for (int k = 1; k < n_features; k++)
			Alpha[k] = 50 / n_features;
		for (int u = 1; u < data.n_nodes; u++)
			Beta[u] = 200 / data.n_nodes;

	}

	private void inferPriorWords(CascadeData data, int n_features) {
		// TODO Auto-generated method stub

	}

	private void inferPriorRates(CascadeData data) {
		// TODO Auto-generated method stub

	}

	private double[][] sampleA(Model model, CascadeData data,
			GibbsSamplerState curr_state, double[][] F_curr,Counters counters) {
		double[][] A_new = new double[data.n_nodes][model.n_features];

		return A_new;
	}

	private double[][] sampleS(Model model, CascadeData data,
			GibbsSamplerState curr_state, double[][] F_curr,Counters counters) {
		double[][] S_new = new double[data.n_nodes][model.n_features];

		return S_new;

	}

	private double[][] samplePhi(Model model, CascadeData data,
			GibbsSamplerState curr_state, double[][] F_curr,Counters counters) {
		double[][] Phi_new = new double[data.n_words][model.n_features];

		return Phi_new;

	}

	private int[] sampleZ(Model model, CascadeData data,
			SparseDoubleMatrix2D y, int[] z, int[] m_k,Counters counters) {
		// TODO Auto-generated method stub
	    
	    
		return null;
	}//sampleZ

	private SparseDoubleMatrix2D sampleY(Model model, CascadeData data,
			SparseDoubleMatrix2D y, int[] z, int[] m_v,Counters counters) {
		// TODO Auto-generated method stub
		return null;
	}

}
