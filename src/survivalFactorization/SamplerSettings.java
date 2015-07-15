package survivalFactorization;

public class SamplerSettings {

	public int seed=13;
	public int n_iterations=100;
	public int burnin=10;
	public int llk_interval=1;
	public boolean compute_llk=true;
    public String curr_state_log="resources/curr_state_info";
	
	public static SamplerSettings getDefaultSettings(){
	    SamplerSettings s=new SamplerSettings();
	    return s;
	}
	
	public SamplerSettings(){}
	
	public SamplerSettings(int seed, int n_iterations, int burnin,
            int llk_interval, boolean compute_llk) {
        this.seed = seed;
        this.n_iterations = n_iterations;
        this.burnin = burnin;
        this.llk_interval = llk_interval;
        this.compute_llk = compute_llk;
    }

	
	
	
}//SamplerSettings
