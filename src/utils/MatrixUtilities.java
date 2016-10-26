package utils;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;

public class MatrixUtilities {

	public static void copy(double[][] src, double[][] dest) {
		for (int i = 0; i < src.length; i++)
			System.arraycopy(src[i], 0, dest[i], 0, dest[i].length);
	}

	public static void copy(int[][] src, int[][] dest) {
		for (int i = 0; i < src.length; i++)
			System.arraycopy(src[i], 0, dest[i], 0, dest[i].length);
	}

	public static void fill(double[][] a, double val) {
		for (int i = 0; i < a.length; i++)
			for (int j = 0; j < a[0].length; j++)
				a[i][j] = val;
	}

	public static void print(double[][] m) {
		System.out.println();
		for (final double[] element : m) {
			for (final double element2 : element) {
				System.out.printf("%1$f", element2);
				System.out.print(" ");
			}
			System.out.println();
		}
		System.out.println();
	}

	public static void print(double[][] m, String filePath) throws IOException {
		System.out.println("Writing matrix in: " + filePath);

		final PrintWriter pw = new PrintWriter(new FileWriter(filePath));

		final DecimalFormat df = new DecimalFormat("#");
		df.setMaximumFractionDigits(32);
		final DecimalFormatSymbols dfs = df.getDecimalFormatSymbols();
		dfs.setDecimalSeparator('.');
		df.setDecimalFormatSymbols(dfs);

		try {
			for (final double[] row : m) {
				int i, n;

				for (i = 0, n = row.length - 1; i < n; ++i)
					pw.print(df.format(row[i]) + "\t");

				pw.println(df.format(row[i]));
			}
		} catch (final Exception ex) {
			System.out.flush();
			ex.printStackTrace();
			System.err.flush();
		} finally {
			pw.close();
		}

		System.out.println("... done!");
	}

	public static void print(int[][] m) {
		System.out.println();
		for (final int[] element : m) {
			for (final int element2 : element) {
				System.out.print(element2);
				System.out.print(" ");
			}
			System.out.println();
		}
		System.out.println();
	}

	public static double[][] product(double[][] a, double[][] b) {
		final double prod[][] = new double[a.length][b[0].length];
		for (int i = 0; i < prod.length; i++)
			for (int j = 0; j < prod[0].length; j++)
				for (int k = 0; k < a[0].length; k++)
					prod[i][j] += a[i][k] * b[k][j];
		return prod;
	}

	public static double[][] transpose(double[][] a) {
		final double[][] tmp = new double[a[0].length][a.length];

		for (int j = 0; j < a[0].length; j++)
			for (int i = 0; i < a.length; i++)
				tmp[j][i] = a[i][j];

		return tmp;
	}

	public static int[][] transpose(int[][] a) {
		final int[][] tmp = new int[a[0].length][a.length];

		for (int j = 0; j < a[0].length; j++)
			for (int i = 0; i < a.length; i++)
				tmp[j][i] = a[i][j];

		return tmp;
	}
}