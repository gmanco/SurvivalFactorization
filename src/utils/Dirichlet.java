package utils;

/**
 * Various useful functions related to Dirichlet distributions.
 *
 * @author Andrew McCallum and David Mimno
 */

public class Dirichlet {

	double magnitude = 1;
	double[] partition;

	Randoms random = null;

	/** A dirichlet parameterized with a single vector of positive reals */
	public Dirichlet(double[] p) {
		magnitude = 0;
		partition = new double[p.length];

		// Add up the total
		for (final double element : p)
			magnitude += element;

		for (int i = 0; i < p.length; i++)
			partition[i] = p[i] / magnitude;
	}

	/**
	 * A symmetric Dirichlet with alpha_i = 1.0 and <code>size</code> dimensions
	 */
	public Dirichlet(int size) {
		this(size, 1.0);
	}

	/**
	 * A symmetric dirichlet: E(X_i) = E(X_j) for all i, j
	 *
	 * @param n
	 *            The number of dimensions
	 * @param alpha
	 *            The parameter for each dimension
	 */
	public Dirichlet(int size, double alpha) {
		magnitude = size * alpha;

		partition = new double[size];

		partition[0] = 1.0 / size;
		for (int i = 1; i < size; i++)
			partition[i] = partition[0];
	}

	private void initRandom() {
		if (random == null)
			random = new Randoms();
	}

	public double[] nextDistribution() {
		final double distribution[] = new double[partition.length];
		initRandom();

		// For each dimension, draw a sample from Gamma(mp_i, 1)
		double sum = 0;

		for (int i = 0; i < distribution.length; i++) {
			distribution[i] = random.nextGamma(partition[i] * magnitude, 1);
			if (distribution[i] <= 0)
				distribution[i] = 0.0001;
			sum += distribution[i];
		}

		// Normalize
		for (int i = 0; i < distribution.length; i++)
			distribution[i] /= sum;

		return distribution;
	}
}