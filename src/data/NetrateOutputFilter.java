package data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.StringTokenizer;

public class NetrateOutputFilter {

	public static void main(String[] args) throws Exception {
		final String in = args[0];
		final String out = args[1];

		final BufferedReader br = new BufferedReader(new FileReader(in));
		final PrintWriter pw = new PrintWriter(new FileWriter(out));

		pw.println("Prediction\tClass");

		String line;

		int negatives = 0;
		int total = 0;
		int nulls = 0;
		int positives = 0;

		while ((line = br.readLine()) != null) {
			++total;

			final StringTokenizer st = new StringTokenizer(line);

			String s1 = null;
			if (st.hasMoreTokens())
				s1 = st.nextToken();

			String s2 = null;
			if (st.hasMoreTokens())
				s2 = st.nextToken();

			String s3 = null;
			if (st.hasMoreTokens())
				s3 = st.nextToken();

			String s4 = null;
			if (st.hasMoreTokens())
				s4 = st.nextToken();

			if (s1 != null && s2 != null && s3 != null && s4 != null
					&& s1.length() > 0 && s2.length() > 0 && s3.length() > 0
					&& s4.length() > 0) {

				final double d3 = Double.parseDouble(s3) + 1e-200;
				if (Double.isNaN(d3)) {
					++nulls;
					continue;
				} else if (d3 <= 0) {
					++negatives;
					continue;
				} else
					++positives;

				pw.println(d3 + "\t" + s4);

			} else
				++nulls;
		}

		br.close();
		pw.close();

		System.out.println("Total: " + total);
		System.out.println("Nulls: " + nulls);
		System.out.println("Positives: " + positives);
		System.out.println("Negatives: " + negatives);
	}
}