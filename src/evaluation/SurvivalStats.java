package evaluation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.StringTokenizer;

import survivalFactorizationEM.SurvivalFactorizationEM_Model;
import utils.Node;
import data.CascadeData;
import data.CascadeEvent;

public class SurvivalStats {

	private static double expInfluence(SurvivalFactorizationEM_Model model,
			int node, int k) {

		return model.A[node][k];
	}

	public static void main(String[] args) throws Exception {
		final String file_events = "resources/datasets/memeTracker/v1/cleaned_activations";
		final String file_content = "resources/datasets/memeTracker/v1/cleaned_content";
		final String modelFile = "resources/datasets/memeTracker/models/v1/meme_new_8f.model";
		final String topicFile = "resources/datasets/memeTracker/models/v1/meme_new_8f.clusters.full";

		final boolean isExponential = true;
		final String influencerFile = "resources/datasets/memeTracker/models/v1/cascadeInflunecers.txt";
		final String delayFile = "resources/datasets/memeTracker/models/v1/delay_avg.txt";
		final String topInflFile = "resources/datasets/memeTracker/models/v1/topInfluencers.txt";

		final SurvivalStats stats = new SurvivalStats(file_events,
				file_content, modelFile, topicFile);
		stats.printOutput(isExponential, influencerFile, delayFile, topInflFile);
	}

	private static double rayleighInfluence(
			SurvivalFactorizationEM_Model model, CascadeEvent prev,
			CascadeEvent curr, int k) {

		return model.A[prev.node][k] * (curr.timestamp - prev.timestamp);
	}

	private static int[] topicAssignment(CascadeData cascadeData,
			String topicFile) throws Exception {
		final int[] topicAssignment = new int[cascadeData.n_cascades];

		final BufferedReader br = new BufferedReader(new FileReader(topicFile));
		String line;
		int c = 0;

		while ((line = br.readLine()) != null) {
			line = line.trim();

			if (line.length() == 0)
				continue;

			final StringTokenizer st = new StringTokenizer(line);
			st.nextToken();
			topicAssignment[c++] = Integer.parseInt(st.nextToken(), 10) - 1;
		}

		br.close();
		return topicAssignment;
	}

	private final CascadeData cascadeData;
	private final SurvivalFactorizationEM_Model model;
	private final int[] topicAssignment;

	public SurvivalStats(String file_events, String file_content,
			String modelFile, String topicFile) throws Exception {

		cascadeData = new CascadeData(file_events, file_content);
		model = SurvivalFactorizationEM_Model.readFromFile(modelFile);
		topicAssignment = topicAssignment(cascadeData, topicFile);
	}

	private int getTopicAssignment(int c) {
		return topicAssignment[c];
	}

	private void predictedInfluencer(PrintWriter pw, boolean isExponential) {

		pw.println("item\tuser\tassignedTopic\tpredInfluencer");

		for (int c = 0; c < cascadeData.getNCascades(); ++c) {
			CascadeEvent prevEvent = null;
			final int k = getTopicAssignment(c);

			int bestInfluencer = -1;
			double influence = -1;

			for (final CascadeEvent currentEvent : cascadeData
					.getCascadeEvents(c)) {

				if (prevEvent == null)
					pw.println(c + 1 + "\t" + currentEvent.node + "\t"
							+ (k + 1) + "\t" + -1);
				else {
					final double value = isExponential ? expInfluence(model,
							prevEvent.node, k) : rayleighInfluence(model,
							prevEvent, currentEvent, k);

					if (influence < value) {
						bestInfluencer = prevEvent.node;
						influence = value;
					}

					pw.println(c + 1 + "\t" + currentEvent.node + "\t"
							+ (k + 1) + "\t" + bestInfluencer);
				}

				prevEvent = currentEvent;
			}
		}
	}

	private void printAvgDelay(PrintWriter pw) {
		final double[] delays = new double[model.nFactors];
		final int[] counters = new int[model.nFactors];

		for (int c = 0; c < cascadeData.getNCascades(); ++c) {
			CascadeEvent prevEvent = null;
			final int k = getTopicAssignment(c);

			for (final CascadeEvent currentEvent : cascadeData
					.getCascadeEvents(c)) {

				if (prevEvent != null) {
					delays[k] += currentEvent.timestamp - prevEvent.timestamp;
					++counters[k];
				}

				prevEvent = currentEvent;
			}
		}

		for (int i = 0; i < model.nFactors; ++i)
			pw.println("Topic:\t" + i + "\tdelay (avg):\t" + delays[i]
					/ counters[i]);
	}

	public void printOutput(boolean isExponential, String influencerFile,
			String delayFile, String topInflFile) throws Exception {

		PrintWriter pw = new PrintWriter(new FileWriter(influencerFile));
		predictedInfluencer(pw, isExponential);
		pw.close();

		pw = new PrintWriter(new FileWriter(delayFile));
		printAvgDelay(pw);
		pw.close();

		pw = new PrintWriter(new FileWriter(topInflFile));
		printTopInfluencers(pw);
		pw.close();

		printTopNInfluencersForTopicK(100, 2,
				"resources/datasets/memeTracker/v1/dictionary_hostDictionary_inverted");
	}

	private void printTopInfluencers(PrintWriter pw) {
		final double[][] A = model.A;

		for (int k = 0; k < A[0].length; ++k) {
			int topUser = -1;
			double bestVal = -1;

			for (int u = 0; u < A.length; ++u)
				if (A[u][k] > bestVal) {
					bestVal = A[u][k];
					topUser = u + 1;
				}

			pw.println("topic:\t" + k + "\tbest influencer:\t" + topUser);
		}
	}

	private void printTopNInfluencersForTopicK(int n, int k, String mapFile)
			throws Exception {

		final double[][] A = model.A;
		final ArrayList<Node> list = new ArrayList<Node>(A.length);

		for (int u = 0; u < A.length; ++u)
			list.add(new Node(u, A[u][k]));

		Collections.sort(list);

		final HashMap<Integer, String> map = new HashMap<Integer, String>(
				A.length);

		final BufferedReader br = new BufferedReader(new FileReader(mapFile));
		String line;

		br.readLine(); // Skip first line

		while ((line = br.readLine()) != null) {
			line = line.trim();

			if (line.length() == 0)
				continue;

			final StringTokenizer st = new StringTokenizer(line);

			map.put(Integer.parseInt(st.nextToken(), 10) - 1, st.nextToken());
		}

		br.close();

		for (int i = list.size() - 1; i >= 0 && i > list.size() - 1 - n; --i) {
			final Node node = list.get(i);
			System.out.println(k + "\t" + map.get(node.getId()) + "\t"
					+ node.getValue());
		}
	}
}