package survivalFactorizationEM;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.LinkedList;
import java.util.Properties;
import java.util.StringTokenizer;

public class FullTestNetworkReconstruction {

	private static double computeLinkPrediction(
			SurvivalFactorizationEM_Model model, int source, int destination) {
		double score = 0.0;
		for (int k = 0; k < model.nFactors; k++)
			score += model.pi[k] * model.S[source][k] * model.A[destination][k];
		return score;
	}

	public static void main(String[] args) throws Exception {
		if (args == null || args.length != 1) {
			printUsage();
			return;
		}

		final String conf = args[0];

		final Properties prop = new Properties();
		prop.load(new FileInputStream(conf));

		final String modelFolder = prop.getProperty("model_folder");
		final String[] modelFiles = parseArray(prop.getProperty("model_files"));
		final String testFile = prop.getProperty("test_file");
		final String outputFolder = prop.getProperty("output_folder");
		final String[] outputFiles = parseArray(prop
				.getProperty("output_files"));

		if (modelFiles.length != outputFiles.length) {
			System.err.println("model_files and"
					+ " output_files must have same length");
			return;
		}

		for (int i = 0; i < modelFiles.length; ++i) {
			final File modelPath = new File(modelFolder, modelFiles[i]);
			final File outputFile = new File(outputFolder, outputFiles[i]);

			System.out.print("Loading model" + modelPath.getName() + "...");

			final SurvivalFactorizationEM_Model model = SurvivalFactorizationEM_Model
					.readFromFile(modelPath.getAbsolutePath());

			System.out.println("Done");

			final PrintWriter pw = new PrintWriter(new FileWriter(outputFile));
			final BufferedReader br = new BufferedReader(new FileReader(
					testFile));
			String line = br.readLine();

			pw.println("Prediction\tClass");

			String tokens[];
			while (line != null) {
				tokens = line.split("\t");
				final int source = Integer.parseInt(tokens[0]) - 1;
				final int destination = Integer.parseInt(tokens[1]) - 1;
				final String label = tokens[2];
				final double score = computeLinkPrediction(model, source,
						destination);

				pw.println("" + score + "\t" + label);

				line = br.readLine();
			}

			pw.flush();
			pw.close();
			br.close();
			System.out.println("DONE");
		}
	}

	private static String[] parseArray(String s) {
		final LinkedList<String> l = new LinkedList<String>();
		final StringTokenizer st = new StringTokenizer(s, ";");

		while (st.hasMoreTokens())
			l.add(st.nextToken());

		final String[] v = new String[l.size()];

		int i = 0;
		for (final String s1 : l)
			v[i++] = s1;

		return v;
	}

	private static void printUsage() {
		System.out.println("usage: <configurationFile>");
		System.out.println("example: \"resources/"
				+ "Weibo_dpu_FullTestNetworkReconstruction_run.properties\"");
	}
}
