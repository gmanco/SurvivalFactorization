package survivalFactorizationEM;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import data.CascadeData;
import data.CascadeEvent;
import data.WordOccurrence;
import utils.Weka_Utils;

public class SurvivalFactorizationEM_Learner {

    private int save_step = SurvivalFactorizationEM_Configuration.DEFAULT_SAVE_STEP;

    public SurvivalFactorizationEM_Model build(CascadeData cascadeData,
            int nFactors, int nMaxIterations) throws Exception {

        System.out.println("Survival Factorization EM  with " + nFactors
                + " latent factors");

        int nWords = cascadeData.getNWords();
        int nVertices = cascadeData.getNNodes();

        // init model
        SurvivalFactorizationEM_Model model = new SurvivalFactorizationEM_Model(
                nVertices, nWords, nFactors);

        model = iterateEM(cascadeData, model, nMaxIterations);
        return model;
    }// build

    private SurvivalFactorizationEM_Model iterateEM(CascadeData cascadeData,
            SurvivalFactorizationEM_Model model, int nMaxIterations) {

        System.out.println("Learning phase:\tstarting");
        System.out.println("#Iteration\tChanges\tTime");

        int saveIteration = 0;
        long initTime = System.currentTimeMillis();

        int nCascades = cascadeData.getNCascades();
        int nFactors = model.getNFactors();
        double[][] gamma = new double[nCascades][nFactors];

        for (int c = 0; c < nCascades; c++) {
            for (int k = 0; k < nFactors; k++) {
                gamma[c][k] = Math.max(
                        SurvivalFactorizationEM_Configuration.eps,
                        SurvivalFactorizationEM_Configuration.randomGen
                        .nextDouble());
            }
            Weka_Utils.normalize(gamma[c]);
        }

        double logLikelihood = computeLogLikelihood(cascadeData, model, gamma);
        System.out.println("Init loglikelihood\t" + logLikelihood);

        for (int iterationsDone = 1; iterationsDone <= nMaxIterations; iterationsDone++) {

            gamma = E_Step(cascadeData, model);

            model = M_Step(cascadeData, model, gamma);

            System.out.println(iterationsDone + "\t"
                    + (System.currentTimeMillis() - initTime));

            saveIteration++;
            if (saveIteration == save_step) {
                saveIteration = 0;
                logLikelihood = computeLogLikelihood(cascadeData, model, gamma);
                System.out.println("Iteration:" + iterationsDone
                        + "\t Likelihood \t" + logLikelihood);
                model.store("tmp_" + iterationsDone + ".model");
            }

        } // for each iteration

        long learningTime = System.currentTimeMillis() - initTime;
        System.out.println("Learning Phase: DONE  ("
                + ((double) learningTime / 1000) + " secs)\n");

        return model;

    }// iterate

    
    private double computeLogLikelihood(CascadeData cascadeData,
            SurvivalFactorizationEM_Model model, double gamma[][]) {
        double logLikelihood = 0;

        double first_term = 0.0;
        double second_term = 0.0;
        double third_term = 0.0;
        double fourth_term = 0.0;
        double fifth_term = 0.0;
        double sixth_term = 0.0;

        int nFactors = model.nFactors;

        double log_pi[] = new double[nFactors];

        for (int k = 0; k < nFactors; k++)
            if (model.pi[k] != 0.0)
                log_pi[k] = Math.log(model.pi[k]);
            else
                throw new RuntimeException("Prior is zero");

        List<CascadeEvent> eventsCurrCascade;
        List<CascadeEvent> prevEventsCurrCascade = new ArrayList<CascadeEvent>();
        Set<Integer> inactiveVertices = new HashSet<Integer>();

        List<WordOccurrence> contentCurrCascade;

        for (int c = 0; c < cascadeData.getNCascades(); c++) {

            // compute first term
            for (int k = 0; k < nFactors; k++) {
                first_term += gamma[c][k] * log_pi[k];
            }

            eventsCurrCascade = cascadeData.getCascadeEvents(c);
            prevEventsCurrCascade.clear();
            inactiveVertices.clear();
            inactiveVertices.addAll(cascadeData.getNodeIds());
            double cumulativeAprev[] = new double[nFactors];
            for (CascadeEvent currentEvent : eventsCurrCascade) {

                for (CascadeEvent prevEvent : prevEventsCurrCascade) {

                    double delta_c_uv = currentEvent.timestamp
                            - prevEvent.timestamp;

                    for (int k = 0; k < nFactors; k++) {

                        // all checks
                        if (cumulativeAprev[k] <= 0)
                            throw new RuntimeException();
                        if (model.A[prevEvent.node][k] <= 0)
                            throw new RuntimeException();
                        if (model.S[currentEvent.node][k] <= 0)
                            throw new RuntimeException();

                        // update second term
                        double etaCurr_k = model.A[prevEvent.node][k]
                                / cumulativeAprev[k];
                        double logSA_k = Math.log(model.A[prevEvent.node][k])
                                + Math.log(model.S[currentEvent.node][k]);
                        second_term += gamma[c][k] * etaCurr_k * logSA_k;

                        // update third term
                        third_term += gamma[c][k] * delta_c_uv
                                * model.S[currentEvent.node][k]
                                        * model.A[prevEvent.node][k];

                    } // for each k

                } // for each previous event

                // add current Event to previous
                prevEventsCurrCascade.add(currentEvent);
                inactiveVertices.remove(currentEvent.node);
                // update cumulative for computing eta
                for (int k = 0; k < nFactors; k++) {
                    cumulativeAprev[k] += model.A[currentEvent.node][k];
                }
            } // for each cascade event

            // compute fourth term (inactive nodes)
            for (int inactiveNode : inactiveVertices) {

                for (CascadeEvent currentEvent : eventsCurrCascade) {
                    for (int k = 0; k < nFactors; k++) {
                        fourth_term += gamma[c][k]
                                * (cascadeData.t_max - currentEvent.timestamp)
                                * model.S[inactiveNode][k]
                                        * model.A[currentEvent.node][k];
                    }
                }

            }

            // compute likelihood content
            contentCurrCascade = cascadeData.getCascadeContent(c);
            int n_c = cascadeData.getLenghtOfCascadeContent(c);
            for (WordOccurrence wo : contentCurrCascade) {
                for (int k = 0; k < nFactors; k++) {
                    if (model.Phi[wo.word][k] <= 0)
                        throw new RuntimeException();

                    fifth_term += gamma[c][k] * wo.cnt
                            * Math.log(model.Phi[wo.word][k])
                            - n_c * model.Phi[wo.word][k];

                }

            }

        } // for each cascade

        // compute fifth term
        int N = cascadeData.getNNodes();
        for (int u = 0; u < cascadeData.getNNodes(); u++) {
            for (int k = 0; k < nFactors; k++) {
                sixth_term += N * (model.S[u][k] + model.A[u][k]);
            }
        }

        double K_squared_plus_one = model.nFactors * model.nFactors + 1.0;
        for (int w = 0; w < cascadeData.getNWords(); w++) {
            for (int k = 0; k < nFactors; k++) {
                fifth_term += K_squared_plus_one * (model.Phi[w][k]);
            }
        }

        logLikelihood = first_term + second_term - third_term - fourth_term
                + fifth_term - sixth_term;

        return logLikelihood;
    }// computeLogLikelihood

    /*
     * E-Step: updates gamma_{c,k}
     */
    private double[][] E_Step(CascadeData cascadeData,
            SurvivalFactorizationEM_Model model) {

        double gammaNew[][] = new double[cascadeData.n_cascades][model.nFactors];
        double[] llkEvents;
        double[] llkContent;
        for (int c = 0; c < cascadeData.getNCascades(); c++) {
            llkEvents = computeLogLikelihoodEvents(cascadeData, c, model);
            llkContent = computeLogLikelihoodContent(cascadeData, c, model);
            for (int k = 0; k < model.nFactors; k++)
                gammaNew[c][k] = llkEvents[k] + llkContent[k];
            gammaNew[c] = Weka_Utils.logs2probs(gammaNew[c]);
        }
        return gammaNew;
    }

    private double[] computeLogLikelihoodEvents(CascadeData cascadeData,
            int cascadeIndex, SurvivalFactorizationEM_Model model) {
        double logLikelihoodEvents[] = new double[model.nFactors];

        List<CascadeEvent> prevEvents = new ArrayList<CascadeEvent>();

        double sumA[] = new double[model.nFactors];
        HashSet<Integer> inactiveNodes = new HashSet<Integer>(
                cascadeData.getNodeIds());
        for (CascadeEvent currentEvent : cascadeData
                .getCascadeEvents(cascadeIndex)) {

            for (CascadeEvent prevEvent : prevEvents) {
                double delta = currentEvent.timestamp - prevEvent.timestamp;
                for (int k = 0; k < model.nFactors; k++) {
                    logLikelihoodEvents[k] -= 
                            delta * model.S[currentEvent.node][k]
                                    * model.A[prevEvent.node][k];
                }
            }

            for (int k = 0; k < model.nFactors; k++) {
            		logLikelihoodEvents[k] += Math.log(model.S[currentEvent.node][k] * sumA[k]);
            }
            
            prevEvents.add(currentEvent);
            inactiveNodes.remove(currentEvent.node);
            for (int k = 0; k < model.nFactors; k++) {
                sumA[k] += model.A[currentEvent.node][k];
            }
        }

        for (int inactiveNode : inactiveNodes) {
            for (CascadeEvent currentEvent : cascadeData
                    .getCascadeEvents(cascadeIndex)) {
                for (int k = 0; k < model.nFactors; k++) {
                    logLikelihoodEvents[k] -= (cascadeData.t_max
                            - currentEvent.timestamp)
                            * model.S[currentEvent.node][k]
                                    * model.A[inactiveNode][k];
                }
            }
        }
        return logLikelihoodEvents;
    }// computeLogLikelihoodEvents

    private double[] computeLogLikelihoodContent(CascadeData cascadeData,
            int cascadeIndex, SurvivalFactorizationEM_Model model) {

        double logLikelihoodContent[] = new double[model.nFactors];

        int n_c = cascadeData.getLenghtOfCascadeContent(cascadeIndex);
        for (WordOccurrence wo : cascadeData.getCascadeContent(cascadeIndex)) {
            for (int k = 0; k < model.nFactors; k++) {
                logLikelihoodContent[k] += wo.cnt
                        * Math.log(model.Phi[wo.word][k])
                        - n_c * model.Phi[wo.word][k];
            }
        }

        return logLikelihoodContent;
    }// computeLogLikelihoodContent

    /*
     * M-Step updates the model
     */
    private SurvivalFactorizationEM_Model M_Step(CascadeData cascadeData,
            SurvivalFactorizationEM_Model model, double[][] gamma) {
       
        double A_new_num[][]=new double[cascadeData.n_nodes][model.nFactors];
        double A_new_den[][]=new double[cascadeData.n_nodes][model.nFactors];
        double A_new[][]=new double[cascadeData.n_nodes][model.nFactors];
        
        double S_new_num[][]=new double[cascadeData.n_nodes][model.nFactors];
        double S_new_den[][]=new double[cascadeData.n_nodes][model.nFactors];
        double S_new[][]=new double[cascadeData.n_nodes][model.nFactors];
        
        double Phi_new_num[][]=new double[cascadeData.n_words][model.nFactors];
        double Phi_new[][]=new double[cascadeData.n_words][model.nFactors];
        
        double pi_new[]=new double[model.nFactors]; 
        
        List<CascadeEvent> eventsCurrCascade;
        List<CascadeEvent> prevEventsCurrCascade = new ArrayList<CascadeEvent>();
        Set<Integer> inactiveVertices = new HashSet<Integer>();

        List<WordOccurrence> contentCurrCascade;
        int length_all_traces=0;
        for (int c = 0; c < cascadeData.getNCascades(); c++) {

            for(int k=0;k<model.nFactors;k++)
                pi_new[k]+=gamma[c][k];
            
            eventsCurrCascade = cascadeData.getCascadeEvents(c);
            prevEventsCurrCascade.clear();
            inactiveVertices.clear();
            inactiveVertices.addAll(cascadeData.getNodeIds());
            double cumulativeAprev[] = new double[model.nFactors];
            for (CascadeEvent currentEvent : eventsCurrCascade) {

                for (CascadeEvent prevEvent : prevEventsCurrCascade) {

                    double delta_c_uv = currentEvent.timestamp
                            - prevEvent.timestamp;

                    for (int k = 0; k < model.nFactors; k++) {

                        double etaCurr_k = model.A[prevEvent.node][k]/ cumulativeAprev[k];
                       
                        //update S_num
                        S_new_num[currentEvent.node][k]+=etaCurr_k*gamma[c][k];
                        //update S_den considering activations
                        S_new_den[currentEvent.node][k]+=gamma[c][k]*delta_c_uv* model.A[prevEvent.node][k];
                        
                        //update A
                        A_new_num[prevEvent.node][k]+=etaCurr_k*gamma[c][k];
                        A_new_den[prevEvent.node][k]+=gamma[c][k]*delta_c_uv* model.S[currentEvent.node][k];

                        
                    } // for each k

                } // for each previous event

                // add current Event to previous
                prevEventsCurrCascade.add(currentEvent);
                inactiveVertices.remove(currentEvent.node);
                // update cumulative for computing eta
                for (int k = 0; k < model.nFactors; k++) {
                    cumulativeAprev[k] += model.A[currentEvent.node][k];
                }
            } // for each cascade event

            // for each inactive node (update S_den)
            for (int inactiveNode : inactiveVertices) {
                for (CascadeEvent currentEvent : eventsCurrCascade) {
                    double delta_c_uv=cascadeData.t_max-currentEvent.timestamp;
                    for (int k = 0; k < model.nFactors; k++) {
                        //update S_den considering non-activations
                        S_new_den[inactiveNode][k]+=gamma[c][k]*delta_c_uv* model.A[currentEvent.node][k];
                        //update A_den considering non-activation
                        A_new_den[currentEvent.node][k]+=gamma[c][k]*delta_c_uv* model.S[inactiveNode][k];
                    }
                }

            }// for each inactive node
            
            // update Phi_new & normalization factor for Phi
            contentCurrCascade = cascadeData.getCascadeContent(c);
            length_all_traces+= cascadeData.getLenghtOfCascadeContent(c);
            for (WordOccurrence wo : contentCurrCascade) {
                for (int k = 0; k < model.nFactors; k++) {
                    Phi_new_num[wo.word][k]+=gamma[c][k]*wo.cnt;
                }
            }
        }// for each cascade
        
        
        // now compute the new model parameters        
        double k_squared_plus_one=model.nFactors*model.nFactors+1;
        for (int k = 0; k < model.nFactors; k++) {
            for(int u=0;u<cascadeData.n_nodes;u++){
                //update S
                S_new[u][k]=S_new_num[u][k]/(S_new_den[u][k]+cascadeData.n_nodes);
                
                if(A_new_den[u][k]==0.0)
                    throw new RuntimeException();
                //update A
                A_new[u][k]=A_new_num[u][k]/(A_new_den[u][k]+cascadeData.n_nodes);
            }
            for(int w=0;w<cascadeData.n_words;w++){
                Phi_new[w][k]=Phi_new_num[w][k]/(length_all_traces+k_squared_plus_one);
            }       
        }
        
        Weka_Utils.normalize(pi_new);
        
        
        //build the new model
        SurvivalFactorizationEM_Model updatedModel=new SurvivalFactorizationEM_Model(cascadeData.n_nodes,cascadeData.n_words,model.nFactors);
        updatedModel.setA(A_new);
        updatedModel.setS(S_new);
        updatedModel.setPhi(Phi_new);
        updatedModel.setPi(pi_new);
        return updatedModel;
        
    }//M-step


}
