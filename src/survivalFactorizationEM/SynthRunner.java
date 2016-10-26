package survivalFactorizationEM;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.Arrays;

import data.CascadeData;
import evaluation.ClusteringEvaluator;

public class SynthRunner {

	private static double computeLinkPrediction(
			SurvivalFactorizationEM_Model model, int source, int destination) {

		double score = 0.0;

		for (int k = 0; k < model.nFactors; k++)
			score += model.pi[k] * model.S[source][k] * model.A[destination][k];

		return score;
	}

	private static void evaluateClusters(File[] files) throws Exception {
		int counter = 0;

		for (final File f : files) {
			final String name = f.getName();

			if (name.endsWith(".cascades_assignments")) {
				final int index = name.indexOf(".", name.indexOf("k") + 1);
				final String radix = f.getParent() + File.separator
						+ name.substring(0, index);
				final String groundTruth = radix + ".clusters";
				final String outputFile = radix + ".clusterComparison";

				// -m <clusterModel> -g <groundTruth>
				final String[] ceArgs = { "-m", f.getAbsolutePath(), "-g",
						groundTruth };

				final PrintStream stdOut = System.out;
				final PrintStream ps = new PrintStream(new FileOutputStream(
						outputFile));
				System.setOut(ps);
				ClusteringEvaluator.main(ceArgs);
				System.out.flush();
				ps.close();
				System.setOut(stdOut);

				System.out.println("DONE: " + counter++);
			}
		}
	}

	private static void generateModels(File[] files) throws Exception {
		for (final File f : files) {
			final String name = f.getName();

			if (name.endsWith(".cascades")) {
				final int index = name.indexOf("k") + 1;
				final int nTopics = Integer.parseInt(
						name.substring(index, name.indexOf(".", index)), 10);

				final CascadeData cascadeData = new CascadeData(
						f.getAbsolutePath(), null);
				cascadeData.getInfo();

				for (final int k : new int[] { 2, 4, 8, 16, 32, 64, 128 }) {
					System.out.println("\r\n\r\nRunning inference (" + k
							+ " factors), file:" + name + "...\r\n\r\n");

					final SurvivalFactorizationEM_LearnerOPT inf = new SurvivalFactorizationEM_LearnerOPT();

					final SurvivalFactorizationEM_Model model = inf.build(
							cascadeData, k, 1000, f.getAbsolutePath()
							+ "_assignments");

					System.out.println("Done.");

					System.out.println("\r\n\r\nSaving model...");

					model.store(f.getAbsolutePath() + "_" + nTopics + "_" + k
							+ "f.model");

					System.out.println("Done.");
				}
			}
		}
	}

	private static void generatePreds(File[] files) throws Exception {
		int counter = 1;

		for (final File f : files) {
			final String name = f.getName();

			if (name.endsWith(".model")) {
				final int index = name.indexOf(".", name.indexOf("k") + 1);

				final String testFile = f.getParent() + File.separator
						+ name.substring(0, index) + ".network.2_hops";
				final String outputFile = f.getParent() + File.separator + name
						+ ".pred";

				final SurvivalFactorizationEM_Model model = SurvivalFactorizationEM_Model
						.readFromFile(f.getAbsolutePath());

				final PrintWriter pw = new PrintWriter(new FileWriter(
						outputFile));
				final BufferedReader br = new BufferedReader(new FileReader(
						testFile));
				String line = br.readLine();

				pw.println("Prediction\tClass");

				String tokens[];

				while ((line = br.readLine()) != null) {
					tokens = line.split("\t");
					final int source = Integer.parseInt(tokens[0]) - 1;
					final int destination = Integer.parseInt(tokens[1]) - 1;
					final String label = tokens[2];
					final double score = computeLinkPrediction(model, source,
							destination);

					pw.println("" + score + "\t" + label);
				}

				pw.close();
				br.close();
				System.out.println("DONE: " + counter++);
			}
		}
	}

	public static void main(String[] args) throws Exception {
		final String dataFolder = args[0];
		final int operation = Integer.parseInt(args[1], 10);

		final File folder = new File(dataFolder);
		final File[] files = folder.listFiles();

		Arrays.sort(files);

		if (operation == 0)
			generateModels(files);
		else if (operation == 1)
			generatePreds(files);
		else if (operation == 2)
			evaluateClusters(files);
	}
}