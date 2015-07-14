package survivalFactorization;

import data.CascadeData;

public class Main {

    /**
     * @param args
     */
    public static void main(String[] args) {
        String file_events= "/Users/barbieri/Dropbox/shared ICAR/SurvivalFactorization/exp/meme_tracker/debug/cleaned_debug_activations";
        String file_content="/Users/barbieri/Dropbox/shared ICAR/SurvivalFactorization/exp/meme_tracker/debug/cleaned_debug_content";
    
        CascadeData data=new CascadeData(file_events, file_content);
        data.getInfo();
        
        GibbsSampler sampler=new GibbsSampler(SamplerSettings.getDefaultSettings());
        sampler.runInference(data, 5);
    }//main

}//Main
