package data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.StringTokenizer;
import java.util.TreeSet;

import cern.colt.matrix.tdouble.impl.SparseDoubleMatrix2D;

public class CascadeData {

	public static void main(String[] args) throws Exception {
		final String file_events = "/Users/barbieri/Dropbox/shared ICAR/SurvivalFactorization/exp/meme_tracker/debug/cleaned_debug_activations";
		final String file_content = "/Users/barbieri/Dropbox/shared ICAR/SurvivalFactorization/exp/meme_tracker/debug/cleaned_debug_content";

		final CascadeData data = new CascadeData(file_events, file_content);
		data.getInfo();

	}

	/*
	 * Provides fast access to t_u(c)
	 */
	protected SparseDoubleMatrix2D cascadeEvents;

	/*
	 * For each cascade it contains an ArrayList of cascade events (sorted)
	 */
	protected List<CascadeEvent>[] cascadeEventsSrt;

	protected Set<Integer>[] cascadesForNode;
	protected List<WordOccurrence>[] cascadeContent;
	protected Set<Integer>[] cascadesForWord;

	protected int[] lengthCascade;

	/*
	 * Set of nodes ids
	 */
	protected Set<Integer> nodeSet = new TreeSet<Integer>();

	/*
	 * Set of cascadesIds
	 */
	protected Set<Integer> cascadeSet = new TreeSet<Integer>();

	/*
	 * Set of word ids
	 */
	protected Set<Integer> wordSet = new TreeSet<Integer>();
	// properties
	public int n_nodes;
	public int n_cascades;
	public int n_words = 0;

	public double t_max;

	public CascadeData(String file_events, String file_content) {
		// here we read the file
		// file_cascade is in the format (u,c,t_u(c))
		// file_content is in the format (w,c,n_{w,c})
		processEventFile(file_events);

		if (file_content != null)
			processContentFile(file_content);
	}

	/*
	 * Return the activation timestamp for the pair (node,cascade). Returns -1
	 * if missing
	 */
	public double getActivationTimestamp(int nodeId, int cascadeId) {
		Double t = cascadeEvents.get(nodeId, cascadeId);
		if (t == null)
			t = -1.0;
		return t;
	}// getTimestamp

	public List<WordOccurrence> getCascadeContent(int c) {
		if (cascadeContent == null || cascadeContent[c] == null)
			return new ArrayList<WordOccurrence>();
		return cascadeContent[c];
	}// WordOccurrence

	public List<CascadeEvent> getCascadeEvents(int cascadeId) {
		List<CascadeEvent> events = cascadeEventsSrt[cascadeId];
		if (events == null)
			events = new ArrayList<CascadeEvent>();
		return events;
	}// getSrtEventsForCascade

	public Set<Integer> getCascadeIds() {
		return cascadeSet;
	}

	/*
	 * Returns a set with cascade ids on which the node is active
	 */
	public Set<Integer> getCascadeIdsForNode(int node) {
		Set<Integer> ris = cascadesForNode[node];
		if (ris == null)
			ris = new HashSet<Integer>();
		return ris;
	}// getCascadeIdsForNode

	public Set<Integer> getCascadeIdsForWord(int word) {
		if (cascadesForWord == null || cascadesForWord[word] == null)
			return new HashSet<Integer>();

		return cascadesForWord[word];
	}// getCascadeIdsForWord

	public void getInfo() {
		System.out.println("*** Statistics cascade data ****");
		System.out.format("Number of nodes: %d\n", n_nodes);
		System.out.format("Number of cascades %d\n", n_cascades);
		System.out.format("Number of words %d\n", n_words);
		System.out.format("T_max %f\n", t_max);
		System.out.println("*********************************");
	}// getInfo

	public int getLenghtOfCascadeContent(int c) {
		if (lengthCascade == null)
			return 0;
		return lengthCascade[c];
	}

	public int getNCascades() {
		return n_cascades;
	}

	public int getNNodes() {
		return n_nodes;
	}

	public Set<Integer> getNodeIds() {
		return nodeSet;
	}

	public int getNWords() {
		return n_words;
	}

	public double getTMax() {
		return t_max;
	}

	public Set<Integer> getWordIds() {
		return wordSet;
	}

	/*
	 * Read Content file
	 */
	@SuppressWarnings("unchecked")
	private void processContentFile(String file_content) {
		System.out.print("Reading cascade content file from " + file_content
				+ " ...");
		final long tic = System.currentTimeMillis();
		try {
			wordSet = new TreeSet<Integer>();

			// read dimensions
			BufferedReader br = new BufferedReader(new FileReader(file_content));
			String line = br.readLine();
			String tokens[];
			line = br.readLine();// skip header

			while (line != null) {
				tokens = line.split("\t");
				final int word = Integer.parseInt(tokens[0]) - 1;
				wordSet.add(word);
				line = br.readLine();
			}
			br.close();

			n_words = wordSet.size();

			cascadeContent = new ArrayList[n_cascades];
			cascadesForWord = new HashSet[n_words];
			lengthCascade = new int[n_cascades];
			br = new BufferedReader(new FileReader(file_content));
			line = br.readLine();
			line = br.readLine();// skip header

			while (line != null) {
				tokens = line.split("\t");
				final int word = Integer.parseInt(tokens[0]) - 1;
				final int cascadeId = Integer.parseInt(tokens[1]) - 1;
				final int cnt = Integer.parseInt(tokens[2]);

				if (cascadeContent[cascadeId] == null)
					cascadeContent[cascadeId] = new ArrayList<WordOccurrence>();
				cascadeContent[cascadeId].add(new WordOccurrence(word, cnt));

				if (cascadesForWord[word] == null)
					cascadesForWord[word] = new HashSet<Integer>();
				cascadesForWord[word].add(cascadeId);
				lengthCascade[cascadeId] += cnt;

				line = br.readLine();
			}
			br.close();

		} catch (final Exception e) {
			e.printStackTrace();
		}
		final long time = (System.currentTimeMillis() - tic) / 1000;
		System.out.print(" Done (" + time + " secs)\n");

	}// processContentFile

	/*
	 * Read event file
	 */
	@SuppressWarnings("unchecked")
	private void processEventFile(String file_events) {

		System.out.print("Reading event file from " + file_events + " ... ");
		final long tic = System.currentTimeMillis();
		try {

			nodeSet = new TreeSet<Integer>();
			cascadeSet = new TreeSet<Integer>();
			t_max = -1;

			// read dimensions
			BufferedReader br = new BufferedReader(new FileReader(file_events));
			String line = br.readLine(); // skip header

			while ((line = br.readLine()) != null) {
				final StringTokenizer st = new StringTokenizer(line);
				final int nodeId = Integer.parseInt(st.nextToken()) - 1;
				final int cascadeId = Integer.parseInt(st.nextToken()) - 1;
				nodeSet.add(nodeId);
				cascadeSet.add(cascadeId);
			}

			br.close();

			// init variables
			n_nodes = nodeSet.size();
			n_cascades = cascadeSet.size();
			cascadeEvents = new SparseDoubleMatrix2D(n_nodes, n_cascades);
			cascadeEventsSrt = new ArrayList[n_cascades];
			cascadesForNode = new HashSet[n_nodes];

			// read cascades
			br = new BufferedReader(new FileReader(file_events));
			line = br.readLine(); // skip header

			while ((line = br.readLine()) != null) {
				final StringTokenizer st = new StringTokenizer(line);
				final int nodeId = Integer.parseInt(st.nextToken()) - 1;
				final int cascadeId = Integer.parseInt(st.nextToken()) - 1;
				final double timestamp = Double.parseDouble(st.nextToken());

				// set the timestamp
				cascadeEvents.set(nodeId, cascadeId, timestamp);

				// add event to CascadeEventsSrt
				if (cascadeEventsSrt[cascadeId] == null)
					cascadeEventsSrt[cascadeId] = new ArrayList<CascadeEvent>();
				cascadeEventsSrt[cascadeId].add(new CascadeEvent(nodeId,
						timestamp));

				if (cascadesForNode[nodeId] == null)
					cascadesForNode[nodeId] = new HashSet<Integer>();
				cascadesForNode[nodeId].add(cascadeId);

				// check t_max
				if (t_max < timestamp)
					t_max = timestamp;
			}

			br.close();

			// sort events within cascades
			for (final int cascadeId : cascadeSet)
				Collections.sort(cascadeEventsSrt[cascadeId]);

		} catch (final Exception e) {
			e.printStackTrace();
		}
		final long time = (System.currentTimeMillis() - tic) / 1000;
		System.out.print(" Done (" + time + " secs)\n");
	}// processEventFile
}