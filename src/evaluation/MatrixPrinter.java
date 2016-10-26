package evaluation;

import java.io.File;
import java.util.TreeSet;

import survivalFactorizationEM.SurvivalFactorizationEM_Model;
import utils.MatrixUtilities;

public class MatrixPrinter {

	private static class Entry implements Comparable<Entry> {
		int index;
		double value;

		public Entry(int index, double value) {
			this.index = index;
			this.value = value;
		}

		// DESC order
		@Override
		public int compareTo(Entry o) {
			return Double.compare(o.value, value);
		}
	}

	@SuppressWarnings("unused")
	private static double[][] columnNormalization(double[][] m) {
		final int rows = m.length;

		if (rows <= 1)
			return null;

		final int cols = m[0].length;

		if (cols < 1)
			return null;

		final double[] maxs = new double[cols];

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				if (m[i][j] > maxs[j])
					maxs[j] = m[i][j];

		final double[][] r = new double[rows][cols];

		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j)
				r[i][j] = m[i][j] / maxs[j];

		return r;
	}

	@SuppressWarnings("unchecked")
	public static double[][] heatMap(double[][] m) {
		if (m == null)
			return null;

		final int rows = m.length;

		if (rows <= 1)
			return null;

		final int cols = m[0].length;

		if (cols < 1)
			return null;

		final TreeSet<Entry>[] counts = new TreeSet[cols];

		for (int j = 0; j < cols; ++j)
			counts[j] = new TreeSet<Entry>();

		for (int i = 0; i < rows; ++i) {
			final int imax = iMax(m[i]);
			counts[imax].add(new Entry(i, m[i][imax]));
		}

		final double[][] r = new double[rows][cols];
		int i = 0;
		int j = 0;

		for (final TreeSet<Entry> set : counts) {
			for (final Entry e : set)
				if (m[e.index][j] != 0)
					System.arraycopy(m[e.index], 0, r[i++], 0, cols);

			++j;
		}

		return r;
	}

	private static int iMax(double[] v) {
		int imax = 0;
		double max = v[0];

		for (int i = 1, n = v.length; i < n; ++i)
			if (max < v[i]) {
				imax = i;
				max = v[i];
			}

		return imax;
	}

	public static void main(String[] args) throws Exception {
		final String modelFolder = args[0];
		final String outputFolder = args[1];

		final File[] files = new File(modelFolder).listFiles();

		for (final File f : files)
			if (f.getName().endsWith(".model")) {
				System.out.print("Loading model" + f.getName() + "...");
				final SurvivalFactorizationEM_Model model = SurvivalFactorizationEM_Model
						.readFromFile(f.getAbsolutePath());
				System.out.println("... done!");

				System.out.println("Matrix A");
				MatrixUtilities.print(heatMap(model.A), outputFolder
						+ File.separator + f.getName() + ".sorted.A");

				System.out.println("Matrix S");
				MatrixUtilities.print(heatMap(model.S), outputFolder
						+ File.separator + f.getName() + ".sorted.S");

				if (model.Phi != null) {
					System.out.println("Matrix Phi");
					MatrixUtilities.print(heatMap(model.Phi), outputFolder
							+ File.separator + f.getName() + ".sorted.phi");
				}

				System.out.println("\r\n... done!");
				System.out.println("\r\n");
			}
	}
}