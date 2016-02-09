package survivalFactorizationEM;

import java.util.List;

import data.CascadeData;
import data.CascadeEvent;

public class SurvivalFactorizationEM_ModelCounters {
	double [][][] A_c_u_k;
	double [][][] tilde_A_c_u_k;
	double [][] A_c_k;
	double [][] tilde_A_c_k;
	
	double [][][] R_c_u_k;
	
	double [] S_k;
	double [][] S_c_k;
	double [][] tilde_S_c_k;
	
	double [][] L_c_k;
	
    public int nVertices;
    public int nFactors;
    public int nCascades;

	public SurvivalFactorizationEM_ModelCounters(int n, int k, int c){
		this.nCascades = c; 
		this.nCascades = k; 
		this.nCascades = n; 
		
		resetCounters();	 
	}

	public void update(CascadeData cascadeData, SurvivalFactorizationEM_Model model) {
		resetCounters();
        List<CascadeEvent> eventsCurrCascade;
		
	       for (int c = 0; c < cascadeData.getNCascades(); c++) {
	            
	    	   eventsCurrCascade = cascadeData.getCascadeEvents(c);
	            for (CascadeEvent currentEvent : eventsCurrCascade) {
		            for(int k=0;k<nFactors;k++){
		            		
		            
		            }
	            }
	       }
		
	}

	private void resetCounters() {
		A_c_u_k = new double[nCascades][nCascades][nCascades];
		tilde_A_c_u_k = new double[nCascades][nCascades][nCascades];
		A_c_k = new double[nCascades][nCascades];
		tilde_A_c_k = new double[nCascades][nCascades];
		
		R_c_u_k = new double[nCascades][nCascades][nCascades];
		
		S_k = new double[nCascades];
		S_c_k = new double[nCascades][nCascades];		
		tilde_S_c_k = new double[nCascades][nCascades];
		L_c_k = new double[nCascades][nCascades];		 
	}
}
