package evaluation;

import java.util.ArrayList;
import java.util.HashSet;

import survivalFactorizationEM.SurvivalFactorizationEM_Model;
import data.CascadeData;
import data.CascadeEvent;

public class CascadeCompletion {

	// Data c
	// chiamiamo c_a la parte iniziale
	// e c_b la parte finale
	// (e.g. c_a contiene il 30/40/50/70%)
	// assumiamo T_c il tempo massimo della cascata
	// per ogni u in c_a:
	// per ogni v non in c_a:
	// calcoliamo t_v = min(t_u,v)
	// dove t_u,v ~ t_u + exp(A_u * S_v)
	// se t_u < T_c allora abbiamo 1
	// altrimenti 0
	// il file risultante è
	// cascade , user , predicted_class, true_class
	// dove true_class = 1 se v è in c_b
	// 0 altrimenti

	public static void evaluate(SurvivalFactorizationEM_Model model,
			CascadeData cascades, double coeff, double maxDelay) {

		final int n_nodes = cascades.n_nodes;

		final HashSet<Integer> uninfectedNodes = new HashSet<Integer>(n_nodes);
		final HashSet<Integer> infectedNodes = new HashSet<Integer>(n_nodes);
		final Integer[] nodeIDs = new Integer[n_nodes];

		for (int i = 0; i < n_nodes; ++i) {
			final Integer n = Integer.valueOf(i);

			uninfectedNodes.add(n);
			nodeIDs[i] = n;
		}

		for (int c = 0; c < cascades.n_cascades; ++c) {
			final ArrayList<CascadeEvent> cascade = (ArrayList<CascadeEvent>) cascades
					.getCascadeEvents(c);

			final int nEvents = cascade.size();
			int i;

			for (i = 0; i < nEvents * coeff; ++i) {
				final CascadeEvent e = cascade.get(i);
				final Integer id = nodeIDs[e.node];

				uninfectedNodes.remove(id);
				infectedNodes.add(id);
			}

			for (i = 0; i < nEvents * coeff; ++i) {
				final CascadeEvent e = cascade.get(i);

				for (final Integer uninfectedNode : uninfectedNodes) {
					final double delay = scalarProduct(model.getA()[e.node],
							model.getS()[uninfectedNode]);

					System.out.println(delay);

					// TODO
				}
			}

			uninfectedNodes.addAll(infectedNodes);
			infectedNodes.clear();
		}
	}

	public static double scalarProduct(double[] v, double[] w) {
		if (v.length != w.length)
			throw new RuntimeException("Invalid scalar product"
					+ " of vectors with different size");

		double sum = 0;

		for (int k = 0; k < v.length; ++k)
			sum += v[k] * w[k];

		return sum;
	}
}