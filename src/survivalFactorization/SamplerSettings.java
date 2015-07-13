package survivalFactorization;

public class SamplerSettings {

	public int seed=13;
	public int n_iterations=100;
	public int burnin=10;
	public int llk_interval=5;
	public boolean compute_llk=true;
    
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
