package data;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Random;
import java.util.StringTokenizer;

public class NetrateInputGenerator {

	public static void main(String[] args) throws Exception {
		final String activationFile = args[0];
		final String networkFile = args[1];
		final String newActivationFile = args[2];
		final String newNetworkFile = args[3];
		final int newNCascades = Integer.parseInt(args[4], 10);
		final int limit = Integer.parseInt(args[5], 10);
		final long seed = Long.parseLong(args[6], 10);

		final Random rnd = new Random(seed);

		final CascadeData cascadeData = new CascadeData(activationFile, null);
		final int nCascades = cascadeData.n_cascades;

		final int[] selectedCascades = selectedCascades(nCascades,
				newNCascades, rnd);

		final HashMap<Integer, Integer> userMap = new HashMap<Integer, Integer>(
				cascadeData.n_nodes);

		int newCascadeIndex = 1;
		int newUserIndex = 1;

		PrintWriter pw = new PrintWriter(new FileWriter(newActivationFile));

		pw.println("UserId\tItemId\tTimestamp");

		for (final int cascadeIndex : selectedCascades) {
			final List<CascadeEvent> eventsCurrCascade = cascadeData
					.getCascadeEvents(cascadeIndex);

			if (eventsCurrCascade.size() > limit)
				continue;

			if (eventsCurrCascade.size() <= 1)
				continue;

			for (final CascadeEvent event : eventsCurrCascade) {
				final int user = event.node;
				final double t = event.timestamp;

				Integer newUser = userMap.get(user);

				if (newUser == null) {
					newUser = newUserIndex++;
					userMap.put(user, newUser);
				}

				pw.println(newUser + "\t" + newCascadeIndex + "\t"
						+ String.format("%.3f", t));
			}

			++newCascadeIndex;

			if (newCascadeIndex > newNCascades)
				break;
		}

		pw.close();

		final BufferedReader br = new BufferedReader(
				new FileReader(networkFile));
		pw = new PrintWriter(new FileWriter(newNetworkFile));
		String line;

		pw.println(br.readLine()); // skip header

		final HashMap<Integer, Integer> badUsers = new HashMap<Integer, Integer>(
				cascadeData.n_nodes);

		while ((line = br.readLine()) != null) {
			final StringTokenizer st = new StringTokenizer(line);

			final int from = Integer.parseInt(st.nextToken(), 10);
			final int to = Integer.parseInt(st.nextToken(), 10);

			Integer newFrom = userMap.get(from);
			Integer newTo = userMap.get(to);

			if (newFrom == null && newTo == null)
				continue;

			if (newFrom == null) {
				newFrom = badUsers.get(from);

				if (newFrom == null) {
					newFrom = newUserIndex++;
					badUsers.put(from, newFrom);
				}
			}

			if (newTo == null) {
				newTo = badUsers.get(to);

				if (newTo == null) {
					newTo = newUserIndex++;
					badUsers.put(to, newTo);
				}
			}

			pw.println(newFrom + "\t" + newTo);
		}

		br.close();
		pw.close();

		System.out.println("#Cascades: " + (newCascadeIndex - 1));
		System.out.println("#Users: " + (badUsers.size() + userMap.size()));
		System.out.println("\t#ActiveUsers: " + userMap.size());
		System.out.println("\t#NonActiveUsers: " + badUsers.size());
	}

	private static int[] selectedCascades(int nCascades, int newNCascades,
			Random rnd) {
		final int[] shuffle = new int[nCascades];

		for (int i = 0; i < nCascades; ++i)
			shuffle[i] = i;

		shuffle(shuffle, rnd);

		return shuffle;
	}

	private static void shuffle(int[] list, Random rnd) {
		for (int i = list.length; i > 1; --i)
			swap(list, i - 1, rnd.nextInt(i));
	}

	private static void swap(int[] arr, int i, int j) {
		final int tmp = arr[i];
		arr[i] = arr[j];
		arr[j] = tmp;
	}
}