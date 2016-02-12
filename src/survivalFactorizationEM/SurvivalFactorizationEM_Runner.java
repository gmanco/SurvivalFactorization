package survivalFactorizationEM;


import java.io.FileInputStream;

import java.util.Properties;

import data.CascadeData;

public class SurvivalFactorizationEM_Runner {

    public static void main(String[] args) throws Exception {
    	
 //   	args =new String[]{"resources/twitter_run.properties"};
        
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
            
    	System.out.print("Reading parameters...");

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


		if (prop.containsKey("max_iterations") )
			nMaxIterations = Integer.parseInt(prop.getProperty("max_iterations"));

		if (prop.containsKey("output") )
			outputFile = prop.getProperty("output");
		else
			outputFile = file_events + "_" + file_content + "_" + nFactors + ".model";
        
    	System.out.println("Done.");

        
        CascadeData cascadeData=new CascadeData(file_events, file_content);
        cascadeData.getInfo();

//        SurvivalFactorizationEM_Learner inf = new SurvivalFactorizationEM_Learner();
        SurvivalFactorizationEM_LearnerOPT inf = new SurvivalFactorizationEM_LearnerOPT();

        SurvivalFactorizationEM_Model model=inf.build(cascadeData, nFactors, nMaxIterations);
       
        System.out.println("Saving the model");
        model.store(outputFile);
        System.out.println("DONE");
    }


    private static void printUsage() {
        System.out.println("-e <eventLog> -c <cascadeContent> -k <nFactors> -o <output> -maxIt <maxIt> ");
    }
}

