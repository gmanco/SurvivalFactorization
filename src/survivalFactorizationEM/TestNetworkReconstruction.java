package survivalFactorizationEM;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.PrintWriter;

public class TestNetworkReconstruction {

    
    public static void main(String[] args) throws Exception{
        String modelPath=args[0];
        String testFile=args[1];
        String outputFile=args[2];
        
        System.out.print("Loading model...");
        SurvivalFactorizationEM_Model model=SurvivalFactorizationEM_Model.readFromFile(modelPath);
        System.out.println("Done");
        
        
        PrintWriter pw=new PrintWriter(new FileWriter(outputFile));
        BufferedReader br=new BufferedReader(new FileReader(testFile));
        String line=br.readLine();
        
        pw.println("Prediction\tClass");
        
        String tokens[];
        while(line!=null){
            tokens=line.split("\t");
            int source=Integer.parseInt(tokens[0])-1;
            int destination=Integer.parseInt(tokens[1])-1;
            String label=tokens[2];
            double score=computeLinkPrediction(model,source,destination);
            
            pw.println(""+score+"\t"+label);
            
            line=br.readLine();
        }
        
        pw.flush();
        pw.close();
        br.close();
        System.out.println("DONE");
    }

    
    private static double computeLinkPrediction(
            SurvivalFactorizationEM_Model model, int source, int destination) {
        double score=0.0;
        for(int k=0;k<model.nFactors;k++){
            score+=model.pi[k]*model.S[source][k]*model.A[destination][k];
        }
        return score;
    }
}
