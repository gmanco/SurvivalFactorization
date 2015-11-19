package survivalFactorizationEM;


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
        
        args=new String[]{"-e","resources/datasets/twitter/activations",
                "-k","3",
                "-o","resources/model3Twitter"
                 };
        
        System.out.println("*** Survival Factorization EM ***");

        String file_events = null;
        String file_content=null;
        int nFactors=-1;
        int nMaxIterations = SurvivalFactorizationEM_Configuration.DEFAULT_N_ITERATIONS;
        String outputFile = null;

        if (args.length == 0) {
            printUsage();
            return;
        }

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
        
        if(outputFile==null){
            System.out.println("Please specify output file");
            return;
        }
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

