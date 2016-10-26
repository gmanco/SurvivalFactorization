package data;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;

public class TrainingTestSplitter {

	private static boolean canBeMoved(int limit, int c,
			HashMap<Integer, Set<Integer>> node_cascades,
			HashMap<Integer, Set<Integer>> cascade_nodes) {

		final Set<Integer> nodes = cascade_nodes.get(c);

		if (nodes.size() <= limit)
			return false;

		for (final int node : nodes)
			if (node_cascades.get(node).size() <= limit)
				return false;

		return true;
	}

	public static void main(String[] args) throws Exception {
		final String file_events = args[0];
		final double splitCoeff = Double.parseDouble(args[1]);
		final String trainingFile = args[2];
		final String testFile = args[3];
		final int minCascadeLength = Integer.parseInt(args[4], 10);
		final long seed = Long.parseLong(args[5], 10);

		split(file_events, splitCoeff, trainingFile, testFile,
				minCascadeLength, seed);
	}

	private static void mapping(CascadeData cascades,
			HashMap<Integer, Set<Integer>> node_cascades,
			HashMap<Integer, Set<Integer>> cascade_nodes) {

		for (int c = 0; c < cascades.n_cascades; ++c)
			for (final CascadeEvent e : cascades.getCascadeEvents(c)) {
				Set<Integer> set = node_cascades.get(e.node);

				if (set == null) {
					set = new TreeSet<Integer>();
					node_cascades.put(e.node, set);
				}

				set.add(c);

				set = cascade_nodes.get(c);

				if (set == null) {
					set = new TreeSet<Integer>();
					cascade_nodes.put(c, set);
				}

				set.add(e.node);
			}
	}

	private static void move(int c, CascadeData cascades,
			List<List<CascadeEvent>> test,
			HashMap<Integer, Set<Integer>> node_cascades,
			HashMap<Integer, Set<Integer>> cascade_nodes) {

		test.add(cascades.getCascadeEvents(c));

		final Set<Integer> nodes = cascade_nodes.remove(c);

		for (final int node : nodes)
			node_cascades.get(node).remove(c);
	}

	private static void print(List<List<CascadeEvent>> cascades, String file)
			throws Exception {

		final PrintWriter pw = new PrintWriter(new FileWriter(file));
		pw.println("User\tCascade\tTimestamp");

		int c = 0;
		for (final List<CascadeEvent> cascade : cascades) {
			for (final CascadeEvent e : cascade)
				pw.println(e.node + 1 + "\t" + (c + 1) + "\t" + e.timestamp);

			++c;
		}

		pw.close();
	}

	private static void shuffle(int[] v, Random r) {
		for (int i = v.length; i > 1; --i) {
			final int otherIndex = r.nextInt(i);
			final int tmp = v[i - 1];
			v[i - 1] = v[otherIndex];
			v[otherIndex] = tmp;
		}
	}

	public static void split(String file_events, double splitCoeff,
			String trainingFile, String testFile, int minCascadeLength,
			long seed) throws Exception {

		final Random r = new Random(seed);
		final CascadeData cascades = new CascadeData(file_events, null);
		final int n_nodes = cascades.n_nodes;
		final int n_cascades = cascades.n_cascades;

		final HashMap<Integer, Set<Integer>> node_cascades = new HashMap<Integer, Set<Integer>>(
				n_nodes);
		final HashMap<Integer, Set<Integer>> cascade_nodes = new HashMap<Integer, Set<Integer>>(
				n_nodes);

		mapping(cascades, node_cascades, cascade_nodes);

		final List<List<CascadeEvent>> training = new LinkedList<List<CascadeEvent>>();
		final List<List<CascadeEvent>> test = new LinkedList<List<CascadeEvent>>();

		final int[] cascadeIndexes = new int[n_cascades];

		for (int i = 0; i < n_cascades; ++i)
			cascadeIndexes[i] = i;

		shuffle(cascadeIndexes, r);

		for (final int c : cascadeIndexes)
			if (test.size() < n_cascades * splitCoeff
					&& canBeMoved(minCascadeLength, c, node_cascades,
							cascade_nodes))
				move(c, cascades, test, node_cascades, cascade_nodes);
			else
				training.add(cascades.getCascadeEvents(c));

		print(training, trainingFile);
		print(test, testFile);
	}
}