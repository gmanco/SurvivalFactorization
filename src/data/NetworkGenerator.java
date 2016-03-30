package data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.HashSet;
import java.util.StringTokenizer;

public class NetworkGenerator {

	private static boolean contains(int u1, int u2,
			HashMap<Integer, HashSet<Integer>> positives) {

		final HashSet<Integer> set = positives.get(u1);

		if (set == null)
			return false;

		return set.contains(u2);
	}

	public static void main(String[] args) throws Exception {
		final String linkFile = args[0];
		final String outputFile = args[1];

		final BufferedReader br = new BufferedReader(new FileReader(linkFile));
		String line;

		final HashSet<Integer> distinct = new HashSet<Integer>();
		final HashMap<Integer, HashSet<Integer>> positives = new HashMap<Integer, HashSet<Integer>>(
				1024);

		br.readLine(); // skip first line

		while ((line = br.readLine()) != null) {
			final StringTokenizer st = new StringTokenizer(line);

			final int u1 = Integer.parseInt(st.nextToken());
			final int u2 = Integer.parseInt(st.nextToken());

			HashSet<Integer> set = positives.get(u1);

			if (set == null) {
				set = new HashSet<Integer>();
				positives.put(u1, set);
			}

			set.add(u2);

			distinct.add(u1);
			distinct.add(u2);
		}

		br.close();

		final PrintWriter pw = new PrintWriter(new FileWriter(outputFile));

		for (final int u1 : distinct)
			for (final int u2 : distinct)
				if (u1 != u2) // skip auto-loop
					if (contains(u1, u2, positives))
						pw.println(u1 + "\t" + u2 + "\t" + 1);
					else
						pw.println(u1 + "\t" + u2 + "\t" + 0);

		pw.close();
	}
}