package data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;
import java.util.StringTokenizer;
import java.util.TreeMap;

public class CascadeFilter {

	private static HashMap<Integer, Integer> buildCascadeMapping(
			TreeMap<Integer, Integer> counts, String mapFile) throws Exception {

		final HashMap<Integer, Integer> mapping = new HashMap<Integer, Integer>(
				counts.size());

		final PrintWriter pw = new PrintWriter(new FileWriter(mapFile));
		pw.println("originalId" + "\t" + "newId");

		int i = 0;
		for (final Map.Entry<Integer, Integer> e : counts.entrySet())
			if (e.getValue() > 1) {
				mapping.put(e.getKey(), ++i);
				pw.println(e.getKey() + "\t" + i);
			}

		pw.close();

		return mapping;
	}

	private static HashMap<Integer, Integer> buildUserMapping(String linkFile,
			TreeMap<Integer, Integer> cascadeCounts, String userMapFile)
			throws Exception {

		final BufferedReader br = new BufferedReader(new FileReader(linkFile));
		String line;

		final HashMap<Integer, Integer> mapping = new HashMap<Integer, Integer>(
				1024);
		int nextIndex = 0;

		final PrintWriter pw = new PrintWriter(new FileWriter(userMapFile));
		pw.println("originalId" + "\t" + "newId");

		br.readLine(); // skip first line

		while ((line = br.readLine()) != null) {
			final StringTokenizer st = new StringTokenizer(line);

			st.nextToken(); // skip EventId
			final int userId = Integer.parseInt(st.nextToken());
			final int c = Integer.parseInt(st.nextToken());

			final Integer i = cascadeCounts.get(c);

			if (i != null && i > 1 && !mapping.containsKey(userId)) {
				mapping.put(userId, ++nextIndex);
				pw.println(userId + "\t" + nextIndex);
			}
		}

		br.close();
		pw.close();

		return mapping;
	}

	public static void main(String[] args) throws Exception {
		final String linkFile = args[0];
		final String outputFile = args[1];
		final String cascadeMapFile = args[2];
		final String userMapFile = args[3];

		BufferedReader br = new BufferedReader(new FileReader(linkFile));
		String line;

		final TreeMap<Integer, Integer> cascadeCounts = new TreeMap<Integer, Integer>();

		br.readLine(); // skip first line

		while ((line = br.readLine()) != null) {
			final StringTokenizer st = new StringTokenizer(line);

			st.nextToken(); // skip EventId
			st.nextToken(); // skip UserId
			final int c = Integer.parseInt(st.nextToken());

			final Integer i = cascadeCounts.get(c);

			if (i == null)
				cascadeCounts.put(c, 1);
			else
				cascadeCounts.put(c, i + 1);
		}

		br.close();

		final HashMap<Integer, Integer> cascadeMapping = buildCascadeMapping(
				cascadeCounts, cascadeMapFile);
		final HashMap<Integer, Integer> userMapping = buildUserMapping(
				linkFile, cascadeCounts, userMapFile);

		br = new BufferedReader(new FileReader(linkFile));
		final PrintWriter pw = new PrintWriter(new FileWriter(outputFile));

		StringTokenizer st = new StringTokenizer(br.readLine());
		st.nextToken(); // skip EventId
		pw.print(st.nextToken() + "\t");
		pw.print(st.nextToken() + "\t");
		pw.print(st.nextToken() + "\t");
		pw.println();

		while ((line = br.readLine()) != null) {
			st = new StringTokenizer(line);

			st.nextToken(); // skip EventId
			final int userId = Integer.parseInt(st.nextToken());
			final int c = Integer.parseInt(st.nextToken());
			final String time = st.nextToken();

			final Integer i = cascadeCounts.get(c);

			if (i != null && i > 1)
				pw.println(userMapping.get(userId) + "\t"
						+ cascadeMapping.get(c) + "\t" + time);
		}

		br.close();
		pw.close();
	}
}