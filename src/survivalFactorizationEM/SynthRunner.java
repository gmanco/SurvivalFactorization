package survivalFactorizationEM;

import java.io.File;
import java.util.Arrays;

import data.CascadeData;

public class SynthRunner {

	private static void generateModels(String dataFolder) throws Exception {
		final File folder = new File(dataFolder);
		final File[] files = folder.listFiles();

		Arrays.sort(files);

		for (final File f : files) {
			final String name = f.getName();
			if (name.endsWith(".cascades")) {
				final int index = name.indexOf("k") + 1;
				final int nTopics = Integer.parseInt(
						name.substring(index, name.indexOf(".", index)), 10);

				final CascadeData cascadeData = new CascadeData(
						f.getAbsolutePath(), null);
				cascadeData.getInfo();

				System.out.println("\r\n\r\nRunning inference (" + nTopics
						+ " factors), file:" + name + "...\r\n\r\n");

				final SurvivalFactorizationEM_LearnerOPT inf = new SurvivalFactorizationEM_LearnerOPT();

				final SurvivalFactorizationEM_Model model = inf.build(
						cascadeData, nTopics, 1000, f.getAbsolutePath()
						+ "_assignments");

				System.out.println("Done.");

				System.out.println("\r\n\r\nSaving model...");

				model.store(f.getAbsolutePath() + "_" + nTopics + "f.model");

				System.out.println("Done.");
			}
		}
	}

	private static void generatePreds(String dataFolder) {
		// TODO Auto-generated method stub

	}

	public static void main(String[] args) throws Exception {
		final String dataFolder = args[0];
		final String flag = args[1];

		if (flag.equalsIgnoreCase("t") || flag.equalsIgnoreCase("true"))
			generateModels(dataFolder);
		else
			generatePreds(dataFolder);
	}
}