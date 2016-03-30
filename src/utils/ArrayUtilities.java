package utils;

import java.util.Arrays;
import java.util.LinkedList;

public class ArrayUtilities {

	public static double SMALL = 1e-6;

	/**
	 *
	 * @param a
	 *            first array
	 * @param b
	 *            second array
	 * @param threshold
	 * @return true if |a - b|/|a| < threshold false otherwise
	 */
	public static boolean arrayConverged(double a[], double b[],
			double threshold) {

		double a2sum = 0;
		double difference2 = 0;

		for (int i = 0; i < a.length; i++) {
			a2sum += a[i] * a[i];
			difference2 += (a[i] - b[i]) * (a[i] - b[i]);
		}
		/**
		 * sqrt ( (a-b)^2/a^2 )< threshold
		 */
		// System.out.println("A2sum::: "+a2sum);
		if (Math.sqrt(difference2 / a2sum) <= threshold)
			return true;
		else
			return false;
	}// arrayConverged

	public static boolean arrayConverged(float a[], float b[], double threshold) {

		double a2sum = 0;
		double difference2 = 0;

		for (int i = 0; i < a.length; i++) {
			a2sum += a[i] * a[i];
			difference2 += (a[i] - b[i]) * (a[i] - b[i]);
		}
		/**
		 * sqrt ( (a-b)^2/a^2 )< threshold
		 */
		// System.out.println("A2sum::: "+a2sum);
		if (Math.sqrt(difference2 / a2sum) < threshold)
			return true;
		else
			return false;
	}// arrayConverged

	public static String arrayToString(double m[]) {
		final StringBuffer sb = new StringBuffer();
		for (final double element : m)
			sb.append("" + element + "\t");
		return sb.toString();
	}

	public static double avg(double[] array) {
		double avg = 0.0;
		for (final double d : array)
			avg += d;
		return avg / array.length;
	}

	// swaps array elements i and j
	public static void exch(int[] a, int i, int j) {
		final int swap = a[i];
		a[i] = a[j];
		a[j] = swap;
	}

	// swaps array elements i and j
	public static void exch(short[] a, int i, int j) {
		final short swap = a[i];
		a[i] = a[j];
		a[j] = swap;
	}

	public static void exp(double[] a) {
		for (int i = 0; i < a.length; i++) {
			if (Double.isInfinite(a[i]))
				System.err.println("@Warning: the " + i
						+ "-th value of the array is infinite");
			else if (Double.isNaN(a[i]))
				System.err.println("@Warning: the " + i
						+ "-th value of the array is NaN");
			a[i] = Math.exp(a[i]);
			// if(a[i]==0.0)a[i]=SMALL;
		}
	}// exp

	/**
	 *
	 * @param a
	 * @return res[0]=index of the max value res[1]=max value;
	 */
	public static double[] findMax(double[] a) {
		double max = Double.MIN_VALUE;
		int index = -1;
		for (int i = 0; i < a.length; i++)
			if (a[i] > max) {
				max = a[i];
				index = i;
			}
		return new double[] { index, max };
	}// findMax

	public static int[] findMax(int[] a) {
		int max = Integer.MIN_VALUE;
		int index = -1;
		for (int i = 0; i < a.length; i++)
			if (a[i] > max) {
				max = a[i];
				index = i;
			}
		return new int[] { index, max };
	}// findMax

	/**
	 *
	 * @param a
	 * @return res[0]=index of the min value res[1]=min value;
	 */
	public static double[] findMin(double[] a) {
		double max = Double.MAX_VALUE;
		int index = -1;
		for (int i = 0; i < a.length; i++)
			if (a[i] < max) {
				max = a[i];
				index = i;
			}
		return new double[] { index, max };
	}// findMin

	public static double HarmonicMean(double[] array) {

		final double num = array.length;
		double den = 0.0;
		for (final double d : array)
			den += 1 / d;
		return num / den;

	}

	public static void log(double[] a) {
		for (int i = 0; i < a.length; i++)
			if (a[i] <= 0)
				throw new RuntimeException(
						"Cannot compute the log of a value less than zero");
			else
				a[i] = Math.log(a[i]);
	}// log

	public static void main(String[] args) {
		final short a[] = { 1, 2, 3, 4, 5 };
		shuffle(a);
		print(a);
	}

	public static double median(double[] array) {
		Arrays.sort(array);
		return array[array.length / 2];
	}

	public static void normalize(double[] a) {
		double sum = 0.0;
		for (final double element : a)
			sum += element;
		if (Double.isNaN(sum))
			throw new RuntimeException(
					"Cannot normalize a vector with NaN elements ");
		if (Double.isInfinite(sum))
			throw new RuntimeException(
					"Cannot normalize a vector. The sum of the elements is infinite ");
		if (sum != 0)
			for (int i = 0; i < a.length; i++)
				a[i] /= sum;
	}// normalize

	public static void normalize(double[] a, double sum) {
		if (Double.isNaN(sum))
			throw new RuntimeException(
					"Cannot normalize a vector with NaN elements ");
		if (Double.isInfinite(sum))
			throw new RuntimeException(
					"Cannot normalize a vector. The sum of the elements is infinite ");
		if (sum != 0)
			for (int i = 0; i < a.length; i++)
				a[i] /= sum;
	}// normalize

	public static void normalize(Double[] a, double sum) {
		if (Double.isNaN(sum))
			throw new RuntimeException(
					"Cannot normalize a vector with NaN elements ");
		if (Double.isInfinite(sum))
			throw new RuntimeException(
					"Cannot normalize a vector. The sum of the elements is infinite ");
		if (sum != 0)
			for (int i = 0; i < a.length; i++)
				a[i] /= sum;
	}// normalize

	public static void print(double[] m) {
		// System.out.println();
		for (final double element : m)
			System.out.printf("" + element + "\t");
		System.out.println();
	}// print

	public static void print(int m[]) {
		System.out.println();
		for (final int element : m)
			System.out.printf("" + element + "\t");
		System.out.println();
	}// print

	public static void print(short[] m) {
		System.out.println();
		for (final short element : m)
			System.out.printf("" + element + "\t");
		System.out.println();
	}// print

	public static void printforExcel(double[] m) {
		// System.out.println();
		for (final double element : m) {
			final String s = (element + "").replace('.', ',');
			System.out.printf("" + s + "\n");
		}
		System.out.println();
	}// print

	/**
	 * Retrive the quartile value from an array .
	 *
	 * @param values
	 *            THe array of data
	 * @param lowerPercent
	 *            The percent cut off. For the lower quartile use 25, for the
	 *            upper-quartile use 75
	 * @return
	 */
	public static double quartile(double[] values, double lowerPercent) {

		if (values == null || values.length == 0)
			throw new IllegalArgumentException(
					"The data array either is null or does not contain any data.");

		// Rank order the values
		final double[] v = new double[values.length];
		System.arraycopy(values, 0, v, 0, values.length);
		Arrays.sort(v);

		final int n = (int) Math.round(v.length * lowerPercent / 100);

		return v[n];

	}

	public static void shadeVector(double a[], double step) {
		for (int i = 0; i < a.length; i++) {
			System.out.print("" + (i + 1) + " ");
			final double a_i = a[i];
			double tempCumulative = 0;

			while (tempCumulative < a_i) {
				System.out.print(":");
				tempCumulative += step;
			}
			System.out.print("\n");
		}
	}

	// take as input an array of strings and rearrange them in random order
	public static void shuffle(int[] a) {
		final int N = a.length;
		for (int i = 0; i < N; i++) {
			final int r = i + (int) (Math.random() * (N - i)); // between i and
			// N-1
			exch(a, i, r);
		}
	}

	// take as input an array of strings and rearrange them in random order
	public static void shuffle(short[] a) {
		final int N = a.length;
		for (int i = 0; i < N; i++) {
			final int r = i + (int) (Math.random() * (N - i)); // between i and
			// N-1
			exch(a, i, r);
		}
	}

	public static double sum(double[] a) {
		double sum = 0.0;
		for (final double element : a)
			sum += element;
		return sum;
	}// sum

	public static double sum(int[] a) {
		double sum = 0.0;
		for (final int element : a)
			sum += element;
		return sum;
	}// sum

	public static void toDoubleArray(double[] values,
			LinkedList<Double> valueList) {
		int i = 0;
		for (final Double double1 : valueList)
			values[i++] = double1;
	}// toDoubleArray

	public static void toIntArray(int[] indexes, LinkedList<Integer> indexList) {
		int i = 0;
		for (final Integer integer : indexList)
			indexes[i++] = integer;
	}// toDoubleArray

}// ArrayUtilities