package survivalFactorizationEM;

import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import data.CascadeData;
import data.CascadeEvent;
import data.WordOccurrence;
import utils.ArrayUtilities;
import utils.Weka_Utils;

public class SurvivalFactorizationEM_LearnerOPT {

	private int save_step = SurvivalFactorizationEM_Configuration.DEFAULT_SAVE_STEP;

	public SurvivalFactorizationEM_Model build(CascadeData cascadeData, int nFactors, int nMaxIterations)
			throws Exception {

		System.out.println("Survival Factorization EM  with " + nFactors + " latent factors");

		System.out.println("Max number of  Iterations: " + nMaxIterations);

		int nWords = cascadeData.getNWords();
		int nVertices = cascadeData.getNNodes();

		// init model
		SurvivalFactorizationEM_Model model = new SurvivalFactorizationEM_Model(nVertices, nWords, nFactors);

		model = iterateEM(cascadeData, model, nMaxIterations);
		return model;
	}// build

	private SurvivalFactorizationEM_Model iterateEM(CascadeData cascadeData, SurvivalFactorizationEM_Model model,
			int nMaxIterations) throws Exception {

		System.out.println("Learning phase: starting at " + new Date());

		int saveIteration = 0;
		long initTime = System.currentTimeMillis();

		int nVertices = cascadeData.getNNodes();
		int nCascades = cascadeData.getNCascades();
		int nFactors = model.getNFactors();
		double[][] gamma = new double[nCascades][nFactors];
		SurvivalFactorizationEM_Model newmodel;

		SurvivalFactorizationEM_ModelCounters counters = new SurvivalFactorizationEM_ModelCounters(nVertices, nFactors);

		for (int c = 0; c < nCascades; c++) {
			for (int k = 0; k < nFactors; k++) {
				gamma[c][k] = Math.max(SurvivalFactorizationEM_Configuration.eps,
						SurvivalFactorizationEM_Configuration.randomGen.nextDouble());
			}
			Weka_Utils.normalize(gamma[c]);
		}

		double logLikelihood = computeLogLikelihood(cascadeData, model, counters, gamma);

		System.out.format("Init loglikelihood\t %.8f\n", logLikelihood);
		double prevlogLikelihood = logLikelihood;
		double improvement = 1;
		int iterationsDone = 1;
		do {// for each iteration
			gamma = E_Step(cascadeData, model, counters);

			newmodel = M_Step(cascadeData, model, gamma, counters);

			saveIteration++;
			if (saveIteration == save_step) {
				saveIteration = 0;
				logLikelihood = computeLogLikelihood(cascadeData, newmodel, counters, gamma);
				improvement = (prevlogLikelihood - logLikelihood) / prevlogLikelihood;
				if (improvement < 0) {
					try {
						model.store("tmp_" + iterationsDone + ".model");
					} catch (Exception e) {
						e.printStackTrace();
					}
					throw new RuntimeException(String.format(
							"Likelihood increasing\n\tCurrent Likelihood:\t %.5f\n\tPrevious Likelihood:\t %.5f\n\tImprovement:\t %.5f",
							logLikelihood, prevlogLikelihood, improvement));
				}
				prevlogLikelihood = logLikelihood;
				model = newmodel;

				System.out.format("Iteration: %d\tLikelihood\t%.5f (Improvement: %.4f)\n", iterationsDone,
						logLikelihood, improvement);

				try {
					model.store("tmp_" + iterationsDone + ".model");
				} catch (Exception e) {
					e.printStackTrace();
				}
			}
			iterationsDone++;
		} while (iterationsDone <= nMaxIterations && Math.abs(improvement) > SurvivalFactorizationEM_Configuration.eps);

		long learningTime = System.currentTimeMillis() - initTime;
		System.out.format("Learning Phase: DONE  (%d iterations, %.0f secs)\n", iterationsDone,
				((double) learningTime / 1000));

		return model;

	}// iterate

	private double computeLogLikelihood(CascadeData cascadeData, SurvivalFactorizationEM_Model model,
			SurvivalFactorizationEM_ModelCounters counters, double gamma[][]) {

		// FIXME: This LogLikelihood uses the logsumexp trick. IT DOES NOT
		// IMPLEMENT THE EXPECTED LIKELIHOOD
		double logLikelihood = 0;
		double llkEvents[];
		double llkContent[];

		double logpriors = 0;

		int nFactors = model.nFactors;
		double constant_sixth_term = 2 * model.nFactors * model.nFactors + 2.0;

		double exponentialTerms[] = new double[nFactors];
		double log_pi[] = new double[nFactors];

		for (int k = 0; k < nFactors; k++) {
			if (model.pi[k] > 0.0)
				log_pi[k] = Math.log(model.pi[k]);
			else
				throw new RuntimeException("Prior is \t" + model.pi[k]);
		}

		for (int c = 0; c < cascadeData.getNCascades(); c++) {
			llkEvents = computeLogLikelihoodEvents(cascadeData, c, model, counters);
			llkContent = computeLogLikelihoodContent(cascadeData, c, model);
			for (int k = 0; k < model.nFactors; k++) {
				exponentialTerms[k] = llkEvents[k] + llkContent[k] + log_pi[k];
			}
			logLikelihood += Weka_Utils.logSumOfExponentials(exponentialTerms);

		}

		// compute sixth_term term
		int N = cascadeData.getNNodes();
		for (int u = 0; u < cascadeData.getNNodes(); u++) {
			for (int k = 0; k < nFactors; k++) {
				logpriors += N * (model.S[u][k] + model.A[u][k]);
			}
		}

		// TODO @manco:check
		for (int w = 0; w < cascadeData.getNWords(); w++) {
			for (int k = 0; k < nFactors; k++) {
				logpriors -= constant_sixth_term * (model.Phi[w][k]) + Math.log(model.Phi[w][k]);
			}
		}

		logLikelihood += logLikelihood - logpriors;

		return logLikelihood;

	}// computeLogLikelihood

	/*
	 * E-Step: updates gamma_{c,k}
	 */
	private double[][] E_Step(CascadeData cascadeData, SurvivalFactorizationEM_Model model,
			SurvivalFactorizationEM_ModelCounters counters) {

		double gammaNew[][] = new double[cascadeData.n_cascades][model.nFactors];
		double[] llkEvents;
		double[] llkContent;
		double log_pi[] = new double[model.nFactors];
		for (int k = 0; k < model.nFactors; k++) {
			if (model.pi[k] <= 0) {
				ArrayUtilities.print(model.pi);
				throw new RuntimeException();
			}
			log_pi[k] = Math.log(model.pi[k]);
		}

		counters.cumulateS(model);
		for (int c = 0; c < cascadeData.getNCascades(); c++) {
			llkEvents = computeLogLikelihoodEvents(cascadeData, c, model, counters);
			llkContent = computeLogLikelihoodContent(cascadeData, c, model);
			for (int k = 0; k < model.nFactors; k++) {
				gammaNew[c][k] = llkEvents[k] + llkContent[k] + log_pi[k];
			}
			gammaNew[c] = Weka_Utils.logs2probs(gammaNew[c]);
		}
		return gammaNew;
	}

	private double[] computeLogLikelihoodEvents(CascadeData cascadeData, int cascadeIndex,
			SurvivalFactorizationEM_Model model, SurvivalFactorizationEM_ModelCounters counters) {
		double logLikelihoodEvents[] = new double[model.nFactors];
		
		// INFO: counters are updated and used on a local (cascade) basis
		counters.updateCountersOnCascade(cascadeData, cascadeIndex, model);


		for (int k = 0; k < model.nFactors; k++) {
			logLikelihoodEvents[k] = counters.L_c_k[k] - (counters.S_k[k] - counters.S_c_k[k])
					* (cascadeData.t_max * counters.A_c_k[k] - counters.tilde_A_c_k[k]);

			CascadeEvent prevEvent = null;
			for (CascadeEvent currentEvent : cascadeData.getCascadeEvents(cascadeIndex)) {

				if (prevEvent != null) {
					if (counters.A_c_u_k[currentEvent.node][k] <= 0) {
						throw new RuntimeException("A_c_u_k contains a negative value");
					}
					logLikelihoodEvents[k] += Math.log(counters.A_c_u_k[currentEvent.node][k])
							- model.S[currentEvent.node][k]
									* (currentEvent.timestamp * counters.A_c_u_k[currentEvent.node][k]
											- counters.tilde_A_c_u_k[currentEvent.node][k]);
				}

				prevEvent = currentEvent;
			}
		}

		return logLikelihoodEvents;

	}// computeLogLikelihoodEvents

	private double[] computeLogLikelihoodContent(CascadeData cascadeData, int cascadeIndex,
			SurvivalFactorizationEM_Model model) {

		double logLikelihoodContent[] = new double[model.nFactors];

		int n_c = cascadeData.getLenghtOfCascadeContent(cascadeIndex);
		if (n_c != 0) {
			for (WordOccurrence wo : cascadeData.getCascadeContent(cascadeIndex)) {
				for (int k = 0; k < model.nFactors; k++) {
					if (model.Phi[wo.word][k] <= 0) {
						ArrayUtilities.print(model.Phi[wo.word]);
						throw new RuntimeException();
					}
					logLikelihoodContent[k] += wo.cnt * Math.log(model.Phi[wo.word][k]);
				}
			}
			// TODO @manco:check
			for (int w = 0; w < cascadeData.getNWords(); w++) {
				for (int k = 0; k < model.nFactors; k++) {
					logLikelihoodContent[k] -= n_c * model.Phi[w][k];
				}
			}
		}

		return logLikelihoodContent;
	}// computeLogLikelihoodContent

	/*
	 * M-Step updates the model
	 */
	private SurvivalFactorizationEM_Model M_Step(CascadeData cascadeData, SurvivalFactorizationEM_Model model,
			double[][] gamma, SurvivalFactorizationEM_ModelCounters counters) throws Exception {

		double A_new_num[][] = new double[cascadeData.n_nodes][model.nFactors];
		double A_new_den[][] = new double[cascadeData.n_nodes][model.nFactors];
		double A_new[][] = new double[cascadeData.n_nodes][model.nFactors];

		double S_new_num[][] = new double[cascadeData.n_nodes][model.nFactors];
		double S_new_den[][] = new double[cascadeData.n_nodes][model.nFactors];
		double S_new[][] = new double[cascadeData.n_nodes][model.nFactors];

		double Phi_new_num[][] = new double[cascadeData.n_words][model.nFactors];
		double Phi_new[][] = new double[cascadeData.n_words][model.nFactors];

		double pi_new[] = new double[model.nFactors];
		Arrays.fill(pi_new, SurvivalFactorizationEM_Configuration.eps);

		List<CascadeEvent> eventsCurrCascade;
		Set<Integer> inactiveVertices = new HashSet<Integer>();

		counters.cumulateS(model);

		List<WordOccurrence> contentCurrCascade;
		int length_all_traces = 0;
		for (int c = 0; c < cascadeData.getNCascades(); c++) {

			// FIXME: counters are updated and used on a local (cascade) basis
			counters.updateCountersOnCascade(cascadeData, c, model);

			for (int k = 0; k < model.nFactors; k++)
				pi_new[k] += gamma[c][k];

			eventsCurrCascade = cascadeData.getCascadeEvents(c);
			inactiveVertices.clear();
			inactiveVertices.addAll(cascadeData.getNodeIds());

			CascadeEvent prevEvent = null;

			for (CascadeEvent currentEvent : eventsCurrCascade) {

				for (int k = 0; k < model.nFactors; k++) {


					if (prevEvent != null) {
						// update S_num 
						S_new_num[currentEvent.node][k] += gamma[c][k];
						// update S_den considering activations
						S_new_den[currentEvent.node][k] += gamma[c][k]
								* (currentEvent.timestamp * counters.A_c_u_k[currentEvent.node][k]
										- counters.tilde_A_c_u_k[currentEvent.node][k]);

						if (S_new_den[currentEvent.node][k] <= 0) {
							throw new RuntimeException(
									"Negative value for S on node " + currentEvent.node + " of cascade " + c);
						}
					}

					// FIXME: the optimization should be alternated and the
					// counters should be updated in the meantime
					// update A_num
					A_new_num[currentEvent.node][k] += gamma[c][k] * model.A[currentEvent.node][k]
							* counters.R_c_u_k[currentEvent.node][k];
					// update A_den
					A_new_den[currentEvent.node][k] += gamma[c][k]
							* ((counters.tilde_S_c_k[k] - counters.tilde_S_c_u_k[currentEvent.node][k])
									- currentEvent.timestamp
											* (counters.S_c_k[k] - counters.S_c_u_k[currentEvent.node][k])
									+ cascadeData.t_max * (counters.S_k[k] - counters.S_c_k[k])
									- currentEvent.timestamp * (counters.S_k[k] - counters.S_c_k[k]));

					if (A_new_den[currentEvent.node][k] < 0) {
						throw new RuntimeException(
								"Negative value for A on node " + currentEvent.node + " of cascade " + c);
					}

				} // for each k
				inactiveVertices.remove(currentEvent.node);
				prevEvent = currentEvent;
			} // for each cascade event

			// for each inactive node (update S_den)
			for (int inactiveNode : inactiveVertices) {
				for (int k = 0; k < model.nFactors; k++) {
					// update S_den considering non-activations
					S_new_den[inactiveNode][k] += gamma[c][k] * (cascadeData.t_max * counters.A_c_u_k[inactiveNode][k]
							- counters.tilde_A_c_u_k[inactiveNode][k]);
				}

			} // for each inactive node

			// update Phi_new & normalization factor for Phi
			contentCurrCascade = cascadeData.getCascadeContent(c);
			length_all_traces += cascadeData.getLenghtOfCascadeContent(c);
			for (WordOccurrence wo : contentCurrCascade) {
				for (int k = 0; k < model.nFactors; k++) {
					Phi_new_num[wo.word][k] += gamma[c][k] * wo.cnt;
				}
			}

		} // for each cascade

		// now compute the new model parameters
		double two_k_squared_plus_2 = 2 * model.nFactors * model.nFactors + 2;
		for (int k = 0; k < model.nFactors; k++) {

			for (int u = 0; u < cascadeData.n_nodes; u++) {
				// update S
				S_new[u][k] = (S_new_num[u][k] + 1) / (S_new_den[u][k] + cascadeData.n_nodes);

				// update A
				A_new[u][k] = (A_new_num[u][k] + 1) / (A_new_den[u][k] + cascadeData.n_nodes);
			}
			for (int w = 0; w < cascadeData.n_words; w++) {
				Phi_new[w][k] = (Phi_new_num[w][k] + 1.0) / (length_all_traces + two_k_squared_plus_2);
			}
		}

		Weka_Utils.normalize(pi_new);

		// build the new model
		SurvivalFactorizationEM_Model updatedModel = new SurvivalFactorizationEM_Model(cascadeData.n_nodes,
				cascadeData.n_words, model.nFactors);
		updatedModel.setA(A_new);
		updatedModel.setS(S_new);
		updatedModel.setPhi(Phi_new);
		updatedModel.setPi(pi_new);

		// counters.update(cascadeData,updatedModel);

		return updatedModel;

	}// M-step

}
