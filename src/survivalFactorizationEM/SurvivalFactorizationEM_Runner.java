package survivalFactorizationEM;

import java.io.FileInputStream;
import java.util.LinkedList;
import java.util.Properties;
import java.util.StringTokenizer;

import data.CascadeData;

public class SurvivalFactorizationEM_Runner {

	public static void main(String[] args) throws Exception {

		// args =new String[]{"resources/twitter_run.properties"};

		System.out.println("*** Survival Factorization EM ***");

		if (args.length == 0) {
			printUsage();

			return;
		}

		String file_events = null;
		String file_content = null;

		int[] nFactors = { SurvivalFactorizationEM_Configuration.DEFAULT_N_FACTORS };
		int nMaxIterations = SurvivalFactorizationEM_Configuration.DEFAULT_N_ITERATIONS;
		String outputFile = SurvivalFactorizationEM_Configuration.DEFAULT_OUTPUT;
		String assignmentFile = SurvivalFactorizationEM_Configuration.DEFAULT_ASSIGNMENT_FILE;

		System.out.print("\r\n\r\nReading parameters...");

		final String conf = args[0];

		final Properties prop = new Properties();
		prop.load(new FileInputStream(conf));

		if (!prop.containsKey("event_file")
				&& !prop.containsKey("content_file")) {

			System.out.println("Cascade Data must be specified. "
					+ "Please specify either events or contents");

			return;
		}

		file_events = prop.getProperty("event_file");
		file_content = prop.getProperty("content_file");

		if (prop.containsKey("n_factors"))
			nFactors = parseArray(prop.getProperty("n_factors"));

		if (prop.containsKey("max_iterations"))
			nMaxIterations = Integer.parseInt(prop
					.getProperty("max_iterations"));

		if (prop.containsKey("output"))
			outputFile = prop.getProperty("output");

		if (prop.containsKey("assignment_file"))
			assignmentFile = prop.getProperty("assignment_file");

		System.out.println("Done.");

		System.out.println("\r\n\r\nLoading data...");

		final CascadeData cascadeData = new CascadeData(file_events,
				file_content);

		System.out.println("Done.");

		cascadeData.getInfo();

		final long time = System.currentTimeMillis();

		for (final int k : nFactors) {
			System.out.println("\r\n\r\nRunning inference (" + k
					+ " factors)...\r\n\r\n");

			final SurvivalFactorizationEM_LearnerOPT inf = new SurvivalFactorizationEM_LearnerOPT();

			final SurvivalFactorizationEM_Model model = inf.build(cascadeData,
					k, nMaxIterations, assignmentFile);

			System.out.println("Done.");

			System.out.println("\r\n\r\nSaving model...");

			model.store(outputFile + "_" + k + "f.model");

			System.out.println("Done.");
		}

		System.out.println("\r\n\r\nAll models were built.");
		System.out.println("Elapsed time: "
				+ time(System.currentTimeMillis() - time));
	}

	private static int[] parseArray(String s) {
		final LinkedList<String> l = new LinkedList<String>();
		final StringTokenizer st = new StringTokenizer(s, ";");

		while (st.hasMoreTokens())
			l.add(st.nextToken());

		final int[] v = new int[l.size()];

		int i = 0;
		for (final String s1 : l)
			v[i++] = Integer.parseInt(s1, 10);

		return v;
	}

	private static void printUsage() {
		System.out.println("usage: <configurationFile>");
		System.out.println("example: \"resources/twitter_run.properties\"");
	}

	private static String time(long l) {
		final int millis = (int) (l % 1000);
		l /= 1000;
		final int secs = (int) (l % 60);
		l /= 60;
		final int mins = (int) (l % 60);
		l /= 60;
		final int hours = (int) (l % 24);
		l /= 24;
		final int days = (int) l;

		final String sDays = days == 1 ? "1 day" : days + " days";
		final String sHours = hours == 1 ? "1 hour" : hours + " hours";
		final String sMins = mins == 1 ? "1 minute" : mins + " minutes";
		final String sSecs = secs == 1 ? "1 second" : secs + " seconds";
		final String sMillis = millis == 1 ? "1 millisecond" : millis
				+ " milliseconds";

		final String comma = ", ";

		return sDays + comma + sHours + comma + sMins + comma + sSecs + comma
				+ sMillis;
	}
}
