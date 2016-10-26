package survivalFactorizationEM;

public class SuperRunner {

	public static void main(String[] args) throws Exception {
		SurvivalFactorizationEM_Runner
				.main(new String[] { "resources/twitter3_s1_run.properties" });
		SurvivalFactorizationEM_Runner
				.main(new String[] { "resources/twitter3_s2_run.properties" });
	}
}