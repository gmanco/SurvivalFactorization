package utils;

public class ExpLogUtils {

	public static void expNormalization(double[] exponents) {
		logExpNormalization(exponents);

		for (int i = 0; i < exponents.length; ++i)
			exponents[i] = Math.exp(exponents[i]);
	}

	public static void logExpNormalization(double[] exponents) {
		final double div = logSumOfExps(exponents);

		for (int i = 0; i < exponents.length; ++i)
			exponents[i] -= div;
	}

	public static void logExpNormalization(double[] exponents,
			double logSumOfExps) {

		for (int i = 0; i < exponents.length; ++i)
			exponents[i] -= logSumOfExps;
	}

	public static double logSum(double[] v) {
		int iMax = 0;
		double max = v[0];

		for (int i = 1; i < v.length; ++i)
			if (v[i] > max) {
				max = v[i];
				iMax = i;
			}

		double sum = 0;
		final double logMax = Math.log(max);

		for (int i = 0; i < iMax; ++i)
			sum += Math.exp(Math.log(v[i]) - logMax);

		for (int i = iMax + 1; i < v.length; ++i)
			sum += Math.exp(Math.log(v[i]) - logMax);

		return logMax + Math.log1p(sum);

	}

	public static double logSumOfExps(double[] exponents) {
		int iMax = 0;
		double max = exponents[0];

		for (int i = 1; i < exponents.length; ++i)
			if (exponents[i] > max) {
				max = exponents[i];
				iMax = i;
			}

		double sum = 0;

		for (int i = 0; i < iMax; ++i)
			sum += Math.exp(exponents[i] - max);

		for (int i = iMax + 1; i < exponents.length; ++i)
			sum += Math.exp(exponents[i] - max);

		return max + Math.log1p(sum);
	}
}