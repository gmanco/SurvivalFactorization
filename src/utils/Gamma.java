package utils;

public class Gamma {

	private static final double[] P = { 676.5203681218851, -1259.1392167224028,
		771.32342877765313, -176.61502916214059, 12.507343278686905,
		-0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7 };

	public static double cumulativeDistribution(double x, int shape, double rate) {
		double sum = 1;
		double factorial = 1;
		int i = 0;

		for (; i < 2; ++i)
			sum -= Math.pow(rate * x, i) * Math.exp(-x * rate);

		if (shape > 0)
			sum -= Math.exp(-x * rate);

		if (shape > 1)
			sum -= rate * x * Math.exp(-x * rate);

		factorial = 1;
		for (i = 2; i < shape; ++i) {
			factorial *= i;
			sum -= Math.pow(rate * x, i) * Math.exp(-x * rate) / factorial;
		}

		return sum;
	}

	public static double densityDistribution(double x, double shape, double rate) {
		if (x < 0 || shape <= 0 || rate <= 0)
			throw new RuntimeException();

		return Math.pow(rate, shape) * Math.pow(x, shape - 1)
				* Math.exp(-x * rate) / function(shape);
	}

	public static int factorial(int z) {
		if (z < 0)
			throw new RuntimeException();

		if (z == 0 || z == 1)
			return 1;

		int prod = 1;

		for (int i = 2; i <= z; ++i)
			prod *= i;

		return prod;
	}

	public static double function(double z) {
		if (z < 0.5)
			return Math.PI / (Math.sin(Math.PI * z) * function(1 - z));
		else {
			z -= 1;
			double x = 0.99999999999980993;

			for (int i = 0; i < P.length; ++i)
				x += P[i] / (z + i + 1);

			final double t = z + P.length - 0.5;

			return Math.sqrt(2 * Math.PI) * Math.pow(t, z + 0.5) * Math.exp(-t)
					* x;
		}
	}

	public static void main(String[] args) {
		for (int i = 5; i <= 10; ++i) {
			System.out.println("Gamma(" + i + ") = " + function(i));
			System.out.println(i - 1 + "! = " + factorial(i - 1));
		}

		System.out.println("Gamma(" + 0.0078125 + ") = " + function(0.0078125));
		System.out.println("Gamma(" + -1 + ") = " + function(-1));
		System.out.println("CumGamma(" + -1 + ") = " + function(-1));
	}
}