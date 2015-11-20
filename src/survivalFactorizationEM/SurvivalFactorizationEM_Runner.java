package survivalFactorizationEM;


import java.io.FileInputStream;

import java.util.Properties;

import data.CascadeData;

public class SurvivalFactorizationEM_Runner {

    public static void main(String[] args) throws Exception {
     /*  MemeTracker
        args=new String[]{"-e","resources/datasets/memeTracker/cleaned_debug_activations",
                           "-c","resources/datasets/memeTracker/cleaned_debug_content",
                           "-k","3",
                           "-o","resources/model"
                            };
                            */
        
 /*       args=new String[]{"-e","resources/datasets/twitter/activations",
                "-k","16",
                "-o","resources/model16Twitter"
                 };*/
        
/*        args=new String[]{"-e","resources/datasets/weibo/weibo_dpu_activations.txt",
                "-k","16",
                "-o","resources/model16Weibo"
                 };*/
        
        System.out.println("*** Survival Factorization EM ***");

        String file_events = null;
        String file_content=null;
        int nFactors=SurvivalFactorizationEM_Configuration.DEFAULT_N_FACTORS;
        int nMaxIterations = SurvivalFactorizationEM_Configuration.DEFAULT_N_ITERATIONS;
        String outputFile =  SurvivalFactorizationEM_Configuration.DEFAULT_OUTPUT;

        if (args.length == 0) {
            printUsage();
            return;
        }
        
        /*

        for (int i = 0; i < args.length; i++) {
            if (args[i].equalsIgnoreCase("--help")) {
                printUsage();
                return;
            }

            if (args[i].equals("-e")) {
                file_events = args[i + 1];
                i++;
            }
            
            if (args[i].equals("-c")) {
                file_content = args[i + 1];
                i++;
            }

            if (args[i].equals("-k")) {
                nFactors = Integer.parseInt(args[i + 1]);
                i++;
            }

            if (args[i].equals("-o")) {
                outputFile = args[i + 1];
                i++;
            }

            if (args[i].equals("-maxIt")) {
                nMaxIterations = Integer.parseInt(args[i + 1]);
                i++;
            }
   

        } // for each arg
   */
        
        
    	System.out.println("Reading parameters...");

    	final String conf = args[0];

    	Properties prop = new Properties();
    	prop.load(new FileInputStream(conf));
    	
    	if (!prop.containsKey("event_file") && !prop.containsKey("content_file")){
    		System.out.println("Cascade Data must be specified. Please specify either events or contents");
            return;
    	}
    	file_events = prop.getProperty("event_file");
	
    	file_content= prop.getProperty("content_file");

    	
		if (prop.containsKey("n_factors"))
			nFactors = Integer.parseInt(prop.getProperty("n_factors"));

		if (prop.containsKey("output") )
			outputFile = prop.getProperty("output"); 
		else
			outputFile = SurvivalFactorizationEM_Configuration.DEFAULT_OUTPUT;

		if (prop.containsKey("max_iterations") )
			Integer.parseInt(prop.getProperty("max_iterations"));

        
        
        CascadeData cascadeData=new CascadeData(file_events, file_content);
        cascadeData.getInfo();

        SurvivalFactorizationEM_Learner inf = new SurvivalFactorizationEM_Learner();

        SurvivalFactorizationEM_Model model=inf.build(cascadeData, nFactors, nMaxIterations);
       
        System.out.println("Saving the model");
        model.store(outputFile);
        System.out.println("DONE");
    }


    private static void printUsage() {
        System.out.println("-e <eventLog> -c <cascadeContent> -k <nFactors> -o <output> -maxIt <maxIt> ");
    }
}

