package survivalFactorizationEM;

import utils.Randoms;

public class SurvivalFactorizationEM_Configuration {

    
	public static Randoms randomGen=new Randoms(131);
    
    public static double eps=0.000001;
    
    public static int DEFAULT_N_ITERATIONS=100;
    public static final int DEFAULT_N_FACTORS = 2;
    public static int DEFAULT_SAVE_STEP=10;
    public static String DEFAULT_OUTPUT="survFactorization.model";
    
}
