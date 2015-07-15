package survivalFactorization;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Set;

import utils.Dirichlet;
import utils.Multinomial;
import utils.Randoms;
import utils.Weka_Utils;
import cern.colt.matrix.tint.impl.SparseIntMatrix2D;
import data.CascadeData;
import data.CascadeEvent;
import data.WordOccurrence;


public class GibbsSampler {
    static double DEFAULT_SHAPE=1.0;
    static double DEFAULT_SCALE=Math.pow(10,-15);
	double a;
	double b;
	double[] C; // n_words
	double[] D; // n_words
	double[] Alpha; // n_features
	double[] Beta; // n_users
	SamplerSettings settings;
	Randoms randomGenerator;
	
	// FIXME: for debug only, remove after checking		
	PrintWriter pwA;
	PrintWriter pwS;
	//
	
	public GibbsSampler(SamplerSettings settings) {
		this.settings = settings;
		this.randomGenerator= new Randoms(settings.seed);
		
		// FIXME: for debug only, remove after checking		
		try {
			pwA=new PrintWriter(new FileWriter("current_info_A.txt"));
			pwS=new PrintWriter(new FileWriter("current_info_S.txt"));

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		//

	}

	public Model[] runInference(CascadeData data, int n_features) {
	    System.out.println();
	    System.out.println("**** RUNNING INFERENCE ****");
		
	    int n_nodes = data.getNNodes();
		int n_words = data.getNWords();
		int n_cascades = data.getNCascades();
		int n_iterations = settings.n_iterations;
		int burnin = settings.burnin;

	
		

		Model[] models = new Model[n_iterations - burnin];

		inferHyperParams(data, n_features);

		// init paramers
		Model model = new Model(n_nodes, n_words, n_features, new HyperParameters(a,b,C,D));
		model.init(randomGenerator);

		long tot_time = 0;

		GibbsSamplerState curr_state = new GibbsSamplerState(n_nodes,
				n_cascades, n_words, n_features);

		double[] p = (new Dirichlet(Alpha)).nextDistribution();

		curr_state.randomInitZ(p);
		
		Counters counters=new Counters(n_cascades, n_nodes, n_features);
		counters.update(data, model,curr_state);
		
		
		// Gibbs sampling
		//
		// Z is n_cascades x 1
		// Y is n_cascades x n_users
		//
		for (int epoch = 0; epoch < n_iterations; epoch++) {
			long t = System.currentTimeMillis();
			
			GibbsSamplerState next_state = sampleNextState(model, data,
					curr_state,counters);
			curr_state = next_state;
			
			counters.update(data, model,curr_state);
			
					
			long tic=System.currentTimeMillis();
	     //   System.out.print("Sampling A...");
			double[][] A_new = sampleA(model, data, curr_state,counters);
			model.setA(A_new);
			
			counters.update(data, model,curr_state);
			
			
	        long toc=(System.currentTimeMillis()-tic)/1000;
	    //    System.out.println(" Done ("+toc+" secs)");


	        tic=System.currentTimeMillis();
        //    System.out.print("Sampling S...");
			double[][] S_new = sampleS(model, data, curr_state,counters);
			model.setS(S_new);
			
			counters.update(data, model,curr_state);

			toc=(System.currentTimeMillis()-tic)/1000;
        //    System.out.println(" Done ("+toc+" secs)");

            tic=System.currentTimeMillis();
        //    System.out.print("Sampling Phi...");
			double[][] Phi_new = samplePhi(model, data, curr_state,counters);
			model.setPhi(Phi_new);
			toc=(System.currentTimeMillis()-tic)/1000;
        //    System.out.println(" Done ("+toc+" secs)");
			
          
            
            // update counters
            counters.update(data, model,curr_state);
            
			if (epoch > burnin) {
			    Model curr_model = new Model(model);
				models[epoch - burnin] = curr_model;
			}

			// update time
			long iteration_time = (System.currentTimeMillis() - t) / 1000;
			tot_time += iteration_time;

			// print out information about iteration
			if (epoch % settings.llk_interval == 0 && settings.compute_llk) {
				double llk = model.computeLLk(data,counters);
				System.out
						.format("Iteration %d completed [elapsed time: %d s (llk: %.5f )].\n",
								epoch, iteration_time, llk);
				
				 /*
                 * Print current state
                 */
				if(settings.curr_state_log!=null)
                    try {
                        curr_state.printSummaryStatus(data,settings.curr_state_log);
                    } catch (Exception e) {
                        e.printStackTrace();
                    }
                
			} else
				System.out.format(
						"Iteration %d completed [elapsed time: %d s].\n",
						epoch, iteration_time);

		}
		// save(output,'Models');

		System.out
				.format("**************************************************\n");
		System.out.format(
				"*      SAMPLING COMPLETED [elapsed time: %.0fs]     *\n",
				tot_time);
		System.out
				.format("**************************************************\n");

		return models;
	}

	double[][] computeF(Model model, CascadeData data) {
		return model.computeFAllCascades(data);
	}

	private GibbsSamplerState sampleNextState(Model model, CascadeData data,
			GibbsSamplerState curr_state,Counters counters) {

		GibbsSamplerState next_state = new GibbsSamplerState(data.n_nodes,
				data.n_cascades, data.n_words, model.n_features);

		// Counters
		int[] M_k = curr_state.M_k;
		int[] M_v = curr_state.M_v;
		int[] Z = curr_state.Z;
		SparseIntMatrix2D Y = curr_state.Y;

		/*
		 * SAMPLE Y
		 */
		long tic=System.currentTimeMillis();
	//	System.out.print("Sampling Y...");
		SparseIntMatrix2D Y_new = sampleY(model, data, Y, Z, M_v);
		long toc=(System.currentTimeMillis()-tic)/1000;
	//	System.out.println(" Done ("+toc+" secs)");
		
		/*
		 * SAMPLE Z
		 */
		tic=System.currentTimeMillis();
	//	System.out.print("Sampling Z...");
		int[] Z_new = sampleZ(model, data, Y_new, Z, M_k,counters);
		toc=(System.currentTimeMillis()-tic)/1000;
	//    System.out.println(" Done ("+toc+" secs)");
		
	    /*
	     * Update the state
	     */
	    next_state.update(data, Z_new, Y_new);

		return next_state;

	}//sampleNextState

	/*
	 * This method initializes hyper-parameters
	 */
	protected void inferHyperParams(CascadeData data, int n_features) {
		
		inferPriorRates(data);

		inferPriorWords(data);

		// these are default values used in topic models
		Alpha = new double[n_features];
		Beta = new double[data.n_nodes];
		for (int k = 0; k < n_features; k++)
			Alpha[k] = (double)50 / n_features;
		
		for (int u = 0; u < data.n_nodes; u++)
			Beta[u] = (double)200 / data.n_nodes;
		
	}//inferHyperParams

	private void inferPriorWords(CascadeData data) {
		this.C=new double[data.n_words];
		Arrays.fill(C,DEFAULT_SCALE);
		this.D=new double[data.n_words];
        Arrays.fill(D,DEFAULT_SCALE);
	}

	private void inferPriorRates(CascadeData data) {
	    this.a=DEFAULT_SHAPE;
        this.b=DEFAULT_SCALE;
	}

	private double[][] sampleA(Model model, CascadeData data,
			GibbsSamplerState curr_state,Counters counters) {
		double[][] A_new = new double[data.n_nodes][model.n_features];
		
		int n_nodes=data.n_nodes;
		int n_features=model.n_features;
		
		
		// loop on the nodes
		for(int u=0;u<n_nodes;u++){
		    
		    double shape_u[]=new double[n_features];
		    double rate_u[]=new double[n_features];
		    Arrays.fill(shape_u, model.hyperParams.a);
		    Arrays.fill(rate_u, model.hyperParams.b);
		        
		    /*
		     * Compute shape
		     */
		    for(int k=0;k<n_features;k++){
		        shape_u[k]+=curr_state.n_k_u_post[u][k]; 
		    }
		    
		    /*
		     * Compute rate
		     */
		    Set<Integer> cascades_u=data.getCascadeIdsForNode(u);
		    for(int c:cascades_u){
		        int k_c=curr_state.Z[c];
		        double t_u=data.getActivationTimestamp(u, c);
		       

		        
		        double contribute_cascade=0.0;
		        contribute_cascade+=counters.F_curr[c][k_c]*(
		        					counters.tilde_S_c_k[c][k_c]
                                - counters.tilde_S_c_u_k[c].get(u,k_c)
		                           + counters.S_c_u_k[c].get(u,k_c)*t_u
		                           +counters.S_k[k_c]*data.t_max
		                           -counters.S_c_k[c][k_c]*data.t_max
		                           -counters.S_k[k_c]*t_u
		                           );
		                         
		        
		        rate_u[k_c]+=contribute_cascade;
		    }//for each cascade on which the user is active
		    
		    
		   // sample from gamma
		   for(int k=0;k<n_features;k++){
		       A_new[u][k]=randomGenerator.nextGamma(shape_u[k],1.0/rate_u[k]);
		       
				// FIXME: for debug only, remove after checking		
		       pwA.println("("+u+","+k+")\t"+shape_u[k]+"\t\t"+rate_u[k]+"\t\t"+(shape_u[k]/rate_u[k])+"\t\t"+A_new[u][k]);
		       //
		   }
		}// for each node
		
		// FIXME: for debug only, remove after checking		
	       pwA.println();
	       pwA.println();
	       //	
	       
	     return A_new;
	}//sampleA

	private double[][] sampleS(Model model, CascadeData data,
			GibbsSamplerState curr_state,Counters counters) {
		
	    double[][] S_new = new double[data.n_nodes][model.n_features];

		int n_nodes=data.n_nodes;
        int n_features=model.n_features;
        
        double shape_u[]=null;
        double rate_u[]=null;
 
        Set<Integer> cascades_u=null;
        
        // loop on the nodes
        for(int u=0;u<n_nodes;u++){
            
            shape_u=new double[n_features];
            rate_u=new double[n_features];
            
            /*
             * Compute rate
             */
            for(int k=0;k<n_features;k++){
                shape_u[k]=curr_state.n_k_u_pre[u][k]+model.hyperParams.a;
                rate_u[k]=counters.Gamma_k[k]*data.t_max-counters.tilde_Gamma_k[k] + model.hyperParams.b;
            }
            
          
            /*
             * Compute scale
             */
            cascades_u=data.getCascadeIdsForNode(u);
            for(int c:cascades_u){
               
                int k_c=curr_state.Z[c];
                double t_u=data.getActivationTimestamp(u, c);
                
                double contribute_cascade=0.0;
            
                contribute_cascade=counters.F_curr[c][k_c]*(
                                    t_u*counters.A_c_u_k[c].get(u, k_c)
                                    - counters.tilde_A_c_u_k[c].get(u, k_c)
                                    - data.t_max*counters.A_c_k[c][k_c]
                                    + counters.tilde_A_c_k[c][k_c]
                                    );
               
                                   
                rate_u[k_c]+=contribute_cascade;
            }//for each cascade on which the user is active
            
            
           // sample from gamma
           for(int k=0;k<n_features;k++){
               S_new[u][k]=randomGenerator.nextGamma(shape_u[k],1.0/rate_u[k]);
				// FIXME: for debug only, remove after checking		
		       pwS.println("("+u+","+k+")\t"+shape_u[k]+"\t\t"+rate_u[k]+"\t\t"+(shape_u[k]/rate_u[k])+"\t\t"+S_new[u][k]);
		       //

           }
           
        }// for each node

			// FIXME: for debug only, remove after checking		
	       pwS.println();
	       pwS.println();
	       //	

        
		return S_new;
	    
	}//sampleS

    private double[][] samplePhi(Model model, CascadeData data,
            GibbsSamplerState curr_state,Counters counters) {
        
        
        double[][] Phi_new = new double[data.n_words][model.n_features];

        int n_features=model.n_features;
        
        
        Set<Integer> cascades_w;
        
        double F_curr[][]=counters.F_curr;
        double Phi[][]=model.getPhi();
        
        for(int w=0;w<data.n_words;w++){
            
            
            double shape_w[]=new double[n_features];
            double rate_w[]=new double[n_features];
            Arrays.fill(rate_w, model.hyperParams.D[w]);
            
            /*
             * Compute rate
             */
            for(int k=0;k<n_features;k++){
                shape_w[k]=curr_state.N_w_k[w][k]+curr_state.C_w_k[w][k]+model.hyperParams.C[w];
            }
            
            cascades_w=data.getCascadeIdsForWord(w);
            
            for(int c:cascades_w){
                
                int k_c=curr_state.Z[c];
                
                //update F_curr[c]
               F_curr[c][k_c]=F_curr[c][k_c]/Phi[w][k_c]; 
              
               int lenghtContent_c=data.getLenghtOfCascadeContent(c);
               
                  
               double contributeA= F_curr[c][k_c]*( counters.A_c_k[c][k_c]*counters.tilde_S_c_k[c][k_c] 
                                                 -counters.A_prod_tilde_S_c_k[c][k_c]
                                                 - counters.tilde_A_c_k[c][k_c]*counters.S_c_k[c][k_c]
                                                 + counters.tilde_A_prod_S_c_k[c][k_c]
                                               ) ;
              double countributeB= F_curr[c][k_c]*(counters.S_k[k_c]-counters.S_c_k[c][k_c])*
                                                       (data.t_max*counters.A_c_k[c][k_c]-counters.tilde_A_c_k[c][k_c]);
                   
                   
              rate_w[k_c]+=contributeA+countributeB+lenghtContent_c; 
                   
                
            }
            
            // sample from gamma
           for(int k=0;k<n_features;k++){
               Phi_new[w][k]=randomGenerator.nextGamma(shape_w[k],1.0/rate_w[k]);
           }
           
           //now update F_curr
           for(int c:cascades_w){      
               int k_c=curr_state.Z[c];
               F_curr[c][k_c]=F_curr[c][k_c]*Phi_new[w][k_c]; 
           }
           
        }//for each word
            
         return Phi_new;
        
    }//samplePhi


	private int[] sampleZ(Model model, CascadeData data,
	        SparseIntMatrix2D Y, int[] z, int[] m_k,Counters counters) {
	  
	    int Z_new[]=new int[data.n_cascades];
	    
	    int n_features=model.n_features;
	    int n_cascades=data.n_cascades;
	           
        Multinomial multinomial;
        
	    for (int c = 0; c < n_cascades; c++){
	        
	        double logProbEvents[]=computeLogProbEvents(model,data,c,Y,counters);
	        double logProbContent[]=computeLogProbContent(model,data,c,counters);
	        
	        //update m_k by removing the old assignment
	        int k_old=z[c];
	        m_k[k_old]--;
	        if(m_k[k_old]<0)
	            throw new RuntimeException(""+c+" "+k_old);
	        
	        double logPrior[]=new double[n_features];
	        for(int k=0;k<n_features;k++){
	            logPrior[k]=Math.log(m_k[k]+Alpha[k]);
	        }
	        
	        double logProbCascade[]=new double[n_features];
	        for(int k=0;k<n_features;k++){
	            logProbCascade[k]=logProbEvents[k]+logProbContent[k]+logPrior[k];
	        }
	        
	        double theta[]=Weka_Utils.logs2probs(logProbCascade);
	        
	        multinomial = new Multinomial(theta);
	        int k_new=multinomial.sample();
	        m_k[k_new]++;
	        Z_new[c]=k_new;
	        
	    }// for each cascade
	      
		return Z_new;
	}//sampleZ

	/*
	 * Computes log Prob(z_c=k| content_c)
	 */
	private double[] computeLogProbContent(Model model, CascadeData data,
            int c, Counters counters) {
      
	    int n_features=model.n_features;
	    List<WordOccurrence> cascadeContent=data.getCascadeContent(c);
	    double Phi[][]=model.getPhi();
        
        double logProbContent[]=new double[n_features];
        int contentLength=0;
        for(WordOccurrence wo:cascadeContent){
            int w=wo.word;
            int n_w_c=wo.cnt;
            contentLength+=n_w_c;
            for(int k=0;k<n_features;k++){
                logProbContent[k]+=n_w_c*Math.log(Phi[w][k]);
            }
        }//for each word
        
        for(int k=0;k<n_features;k++){
            logProbContent[k]-=contentLength*counters.Phi_k[k];
        }
        return logProbContent;
    }//computeLogProbContent

    /*
	 * Compute log Prob(Z_c=k| events_c)
	 */
	private double[] computeLogProbEvents(Model model, CascadeData data,
            int c,SparseIntMatrix2D Y, Counters counters) {

	    int n_features=model.n_features;
       
        double A[][]=model.getA();
        double S[][]=model.getS();
        
        double firstComponent[]=new double[n_features];
        double secondComponent[]=new double[n_features];
        
        double thirdComponent_A[]=new double[n_features];
        double thirdComponent_B[]=new double[n_features];
        double thirdComponent_C[]=new double[n_features];
        double thirdComponent_D[]=new double[n_features];
        double thirdComponent_E[]=new double[n_features];

        
        
        List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);       
        
        int n_events_cascade=cascadeEvents.size();
    
        double F_c[]=counters.F_curr[c];
        
        for(int k=0;k<n_features;k++){
            firstComponent[k]=(n_events_cascade-1)*Math.log(F_c[k]);
            thirdComponent_A[k]=counters.tilde_S_c_k[c][k]*counters.A_c_k[c][k];
            thirdComponent_C[k]=counters.tilde_A_c_k[c][k]*counters.S_c_k[c][k];
            thirdComponent_E[k]=(counters.S_k[k]-counters.S_c_k[c][k])*
                                    (data.t_max*counters.A_c_k[c][k]-counters.tilde_A_c_k[c][k]);
        }
        
        // scan over events
        for(int e=0;e<n_events_cascade;e++){
            CascadeEvent ce=cascadeEvents.get(e);
            int u=ce.node;
            double t_u=ce.timestamp;

            for(int k=0;k<n_features;k++){
               if(e>0){
                   //TODO  check
                   // id influencer
                   int v = (int) (Y.get(c, u));
                   if(A[v][k]*S[u][k]>0)
                       secondComponent[k]+=Math.log(A[v][k]*S[u][k]);
                   else throw new RuntimeException();
               }
               thirdComponent_B[k]+=A[u][k]*counters.tilde_S_c_u_k[c].get(u,k); 
               thirdComponent_D[k]+=A[u][k]*t_u*counters.S_c_u_k[c].get(u,k); 
                
            }//for each k
          
        }// for each event
        
        double thirdComponent[]=new double[n_features];
        for(int k=0;k<n_features;k++){
            thirdComponent[k]=-F_c[k]*(thirdComponent_A[k]
                                        -thirdComponent_B[k]
                                        -thirdComponent_C[k]
                                        +thirdComponent_D[k]
                                        +thirdComponent_E[k]);
        }
        
      double logProbEvents[]=new double[n_features];
 
      for(int k=0;k<n_features;k++){
          logProbEvents[k]=firstComponent[k]+secondComponent[k]+thirdComponent[k];
      }
        
      return logProbEvents;
    }//computeLogProbEvents

	/*
	 * Y(c,u) contains the id of the influencer node.
	 */
    private SparseIntMatrix2D sampleY(Model model, CascadeData data,
            SparseIntMatrix2D Y, int[] z, int[] m_v) {

        int n_nodes=data.n_nodes;
        int n_cascades=data.n_cascades;
        
        double A[][]=model.getA();
        
        
        SparseIntMatrix2D Y_new =new SparseIntMatrix2D(n_cascades,n_nodes);
        
        Multinomial multinomial;
        for(int c=0;c<n_cascades;c++){
            List<CascadeEvent> cascadeEvents=data.getCascadeEvents(c);
            int n_events_cascade=cascadeEvents.size();
            
            // get the active component for cascade c
            int k_c=z[c];
           
            // loop on events
            for(int e=0;e<n_events_cascade;e++){
                CascadeEvent ce=cascadeEvents.get(e);
                int u=ce.node;
                double t_u=ce.timestamp;
 
                if(e>0){ // skip first activation
                   
                    /*
                     * Remove the old influencer from counters
                     */
                    int old_influencer=Y.get(c,u);
                    m_v[old_influencer]=Math.max(m_v[old_influencer]-1,0);
                    
                    //the number of previous events is (e-1) and they go from 0 to e
                    double probInfluencers[]= new double[e];


					// loop on previous events
					for (int e1 = 0; e1 < e; e1++) {
						CascadeEvent e_prime = cascadeEvents.get(e1);
						int v = e_prime.node;
						double t_v = e_prime.timestamp;
						if (t_u - t_v == 0) {
							probInfluencers[e1] = 0;
						} else if (m_v[v] + Beta[v] > 0)
							probInfluencers[e1] = A[v][k_c]
									* (m_v[v] + Beta[v]);
						else
							throw new RuntimeException();
					}// for each previous event
                    
                    multinomial=new Multinomial(probInfluencers);
                   
                    //sample index of the event that triggered the activation
                    int e1=multinomial.sample();
                    
                    //get the corresponding node
                    int v=cascadeEvents.get(e1).node;
                    
                    //update counter of influence for v
                    m_v[v]++;
                    //set influencer
                    Y_new.set(c,u,v);
                }// (e>0)
                
              
            
            }//for each event
            
            
        }// for each cascade

        
        return Y_new;
	}//sampleY

   
}//GibbsSampler
