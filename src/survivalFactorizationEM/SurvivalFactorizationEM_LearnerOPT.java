package survivalFactorizationEM;

import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import utils.ArrayUtilities;
import utils.ExpLogUtils;
import utils.Randoms;
import utils.Weka_Utils;
import data.CascadeData;
import data.CascadeEvent;
import data.WordOccurrence;

public class SurvivalFactorizationEM_LearnerOPT {

	private static void resetMatrix(double[][] m) {
		for (final double[] element : m)
			Arrays.fill(element, 0);
	}

	private final int save_step;
	private final double eps;
	private final Randoms randomGen;

	private double[][] A_new_num;
	private double[][] A_new_den;
	private double[][] A_new;

	private double[][] S_new_num;
	private double[][] S_new_den;
	private double[][] S_new;

	private double[][] Phi_new_num;
	private double[] Phi_new_den;
	private double[][] Phi_new;

	private double[] pi_new;

	private double[][] llk_exponentialTerms;

	private double exponentialPrior_rate;
	private double priorGamma_scale;
	private double priorGamma_rate;

	private final boolean enablePriorComponent;
	private final boolean enableContentLikelihood;
	private final boolean storeTemporaryModels;

	public SurvivalFactorizationEM_LearnerOPT() {
		eps = SurvivalFactorizationEM_Configuration.eps;
		save_step = SurvivalFactorizationEM_Configuration.DEFAULT_SAVE_STEP;
		randomGen = SurvivalFactorizationEM_Configuration.randomGen;
		enablePriorComponent = true;
		enableContentLikelihood = true;
		storeTemporaryModels = false;
	}

	public SurvivalFactorizationEM_Model build(CascadeData cascadeData,
			int nFactors, int nMaxIterations) throws Exception {

		System.out.println("Survival Factorization EM  with " + nFactors
				+ " latent factors");

		System.out.println("Max number of  Iterations: " + nMaxIterations);

		final int nWords = cascadeData.getNWords();
		final int nVertices = cascadeData.getNNodes();

		exponentialPrior_rate = nVertices;
		priorGamma_scale = 1D + 1D / nFactors;
		priorGamma_rate = 2 * nFactors * nFactors + 2;

		// init model
		SurvivalFactorizationEM_Model model = new SurvivalFactorizationEM_Model(
				nVertices, nWords, nFactors);

		buildUtilityVariables(cascadeData, model);

		model = iterateEM(cascadeData, model, nMaxIterations);

		return model;
	}

	private void buildUtilityVariables(CascadeData cascadeData,
			SurvivalFactorizationEM_Model model) {

		A_new_num = new double[cascadeData.n_nodes][model.nFactors];
		A_new_den = new double[cascadeData.n_nodes][model.nFactors];
		A_new = new double[cascadeData.n_nodes][model.nFactors];

		S_new_num = new double[cascadeData.n_nodes][model.nFactors];
		S_new_den = new double[cascadeData.n_nodes][model.nFactors];
		S_new = new double[cascadeData.n_nodes][model.nFactors];

		Phi_new_num = new double[cascadeData.n_words][model.nFactors];
		Phi_new_den = new double[model.nFactors];
		Phi_new = new double[cascadeData.n_words][model.nFactors];

		pi_new = new double[model.nFactors];

		llk_exponentialTerms = new double[cascadeData.n_cascades][model.nFactors];
	}

	// Something's gone wrong
	// private double computeExpectedLogLikelihood(CascadeData cascadeData,
	// SurvivalFactorizationEM_Model model,
	// SurvivalFactorizationEM_ModelCounters counters, double gamma[][]) {
	//
	// final int nFactors = model.nFactors;
	// final int nCascades = cascadeData.getNCascades();
	// final double log_pi[] = new double[nFactors];
	//
	// double logLikelihood = 0;
	//
	// for (int k = 0; k < nFactors; k++)
	// if (model.pi[k] >= 0.0)
	// log_pi[k] = Math.log(model.pi[k]);
	// else
	// throw new RuntimeException("Prior is \t" + model.pi[k]);
	//
	// counters.cumulateS(model);
	//
	// for (int c = 0; c < nCascades; c++) {
	// counters.updateCountersOnCascade(cascadeData, c, model);
	//
	// for (int k = 0; k < nFactors; k++) {
	// if (gamma[c][k] == 0)
	// continue;
	//
	// double sum = log_pi[k];
	//
	// CascadeEvent prevEvent = null;
	// double prevSum = 0;
	//
	// for (final CascadeEvent currentEvent : cascadeData
	// .getCascadeEvents(c)) {
	//
	// if (prevEvent != null) {
	// sum += Math.log(model.S[currentEvent.node][k]);
	// sum += prevSum;
	//
	// double acukDiff = currentEvent.timestamp
	// * counters.A_c_u_k[currentEvent.node][k]
	// - counters.tilde_A_c_u_k[currentEvent.node][k];
	//
	// if (acukDiff < 0)
	// if (acukDiff < SurvivalFactorizationEM_ModelCounters.NEGATIVE_TOLLERANCE)
	// throw new RuntimeException("Inconsistent value");
	// else
	// acukDiff = 0;
	//
	// sum -= model.S[currentEvent.node][k] * acukDiff;
	// }
	//
	// prevSum += counters.A_c_u_k[currentEvent.node][k] == 0 ? 0
	// : Math.log(model.A[currentEvent.node][k])
	// * model.A[currentEvent.node][k]
	// / counters.A_c_u_k[currentEvent.node][k];
	//
	// prevEvent = currentEvent;
	// }
	//
	// double ackDiff = cascadeData.t_max * counters.A_c_k[k]
	// - counters.tilde_A_c_k[k];
	//
	// if (ackDiff < 0)
	// if (ackDiff < SurvivalFactorizationEM_ModelCounters.NEGATIVE_TOLLERANCE)
	// throw new RuntimeException("Inconsistent value");
	// else
	// ackDiff = 0;
	//
	// double skDiff = counters.S_k[k] - counters.S_c_k[k];
	//
	// if (skDiff < 0)
	// if (skDiff < SurvivalFactorizationEM_ModelCounters.NEGATIVE_TOLLERANCE)
	// throw new RuntimeException("Inconsistent value");
	// else
	// skDiff = 0;
	//
	// sum -= skDiff * ackDiff;
	//
	// if (enableContentLikelihood) {
	// final int n_c = cascadeData.getLenghtOfCascadeContent(c);
	//
	// if (n_c != 0) {
	// for (final WordOccurrence wo : cascadeData
	// .getCascadeContent(c))
	// sum += wo.cnt * Math.log(model.Phi[wo.word][k]);
	//
	// for (int w = 0; w < cascadeData.getNWords(); w++)
	// sum -= n_c * model.Phi[w][k];
	// }
	// }
	//
	// logLikelihood += sum;
	// }
	// }
	//
	// // Prior component
	// if (enablePriorComponent) {
	// double logpriors1 = 0;
	// double logpriors2 = 0;
	//
	// for (int u = 0; u < cascadeData.getNNodes(); u++)
	// for (int k = 0; k < nFactors; k++)
	// logpriors1 += model.S[u][k] + model.A[u][k];
	//
	// logpriors1 *= cascadeData.getNNodes();
	//
	// if (enableContentLikelihood) {
	// for (int w = 0; w < cascadeData.getNWords(); w++)
	// for (int k = 0; k < nFactors; k++)
	// logpriors2 += model.Phi[w][k]
	// + Math.log(model.Phi[w][k]);
	//
	// logpriors2 *= 2 * nFactors * nFactors + 2;
	// }
	//
	// logLikelihood -= logpriors1 + logpriors2;
	// }
	//
	// return logLikelihood;
	// }

	private double computeLogLikelihood(CascadeData cascadeData,
			SurvivalFactorizationEM_Model model,
			SurvivalFactorizationEM_ModelCounters counters, double gamma[][]) {

		final int nFactors = model.nFactors;
		final int nCascades = cascadeData.getNCascades();
		final double log_pi[] = new double[nFactors];

		double logLikelihood = 0;

		for (int k = 0; k < nFactors; k++)
			if (model.pi[k] >= 0.0)
				log_pi[k] = Math.log(model.pi[k]);
			else
				throw new RuntimeException("Prior is \t" + model.pi[k]);

		counters.cumulateS(model);

		for (int c = 0; c < nCascades; c++) {
			final double[] llkEvents = computeLogLikelihoodEvents(cascadeData,
					c, model, counters);
			final double[] llkContent = computeLogLikelihoodContent(
					cascadeData, c, model);

			for (int k = 0; k < model.nFactors; k++)
				llk_exponentialTerms[c][k] = llkEvents[k] + llkContent[k]
						+ log_pi[k];

			logLikelihood += ExpLogUtils.logSumOfExps(llk_exponentialTerms[c]);
		}

		// Prior component
		if (enablePriorComponent) {
			double logpriors1 = 0;
			double logpriors2 = 0;

			for (int u = 0; u < cascadeData.getNNodes(); u++)
				for (int k = 0; k < nFactors; k++)
					logpriors1 += model.S[u][k] + model.A[u][k];

			logpriors1 *= -exponentialPrior_rate;

			if (enableContentLikelihood) {
				double logSum = 0;
				double sum = 0;

				for (int w = 0; w < cascadeData.getNWords(); w++)
					for (int k = 0; k < nFactors; k++) {
						logSum += Math.log(model.Phi[w][k]);
						sum += model.Phi[w][k];
					}

				logpriors2 = (priorGamma_scale - 1) * logSum - priorGamma_rate
						* sum;
			}

			logLikelihood += logpriors1 + logpriors2;
		}

		return logLikelihood;
	}

	private double[] computeLogLikelihoodContent(CascadeData cascadeData,
			int cascadeIndex, SurvivalFactorizationEM_Model model) {

		if (!enableContentLikelihood)
			return new double[model.nFactors];

		final double[] logLikelihoodContent = new double[model.nFactors];
		final int n_c = cascadeData.getLenghtOfCascadeContent(cascadeIndex);

		if (n_c != 0) {
			for (final WordOccurrence wo : cascadeData
					.getCascadeContent(cascadeIndex))
				for (int k = 0; k < model.nFactors; k++) {
					if (model.Phi[wo.word][k] <= 0) {
						ArrayUtilities.print(model.Phi[wo.word]);
						throw new RuntimeException();
					}

					logLikelihoodContent[k] += wo.cnt
							* Math.log(model.Phi[wo.word][k]);
				}

			for (int w = 0; w < cascadeData.getNWords(); w++)
				for (int k = 0; k < model.nFactors; k++)
					logLikelihoodContent[k] -= n_c * model.Phi[w][k];
		}

		return logLikelihoodContent;
	}

	private double[] computeLogLikelihoodEvents(CascadeData cascadeData,
			int cascadeIndex, SurvivalFactorizationEM_Model model,
			SurvivalFactorizationEM_ModelCounters counters) {

		final double logLikelihoodEvents[] = new double[model.nFactors];

		counters.updateCountersOnCascade(cascadeData, cascadeIndex, model);

		loop: for (int k = 0; k < model.nFactors; k++) {
			if (Double.isInfinite(counters.L_c_k[k])) {
				logLikelihoodEvents[k] = counters.L_c_k[k];
				continue;
			}

			double ackDiff = cascadeData.t_max * counters.A_c_k[k]
					- counters.tilde_A_c_k[k];

			if (ackDiff < 0)
				if (ackDiff < SurvivalFactorizationEM_ModelCounters.NEGATIVE_TOLERANCE)
					throw new RuntimeException("Inconsistent value");
				else
					ackDiff = 0;

			logLikelihoodEvents[k] = counters.L_c_k[k]
					- (counters.S_k[k] - counters.S_c_k[k]) * ackDiff;

			CascadeEvent prevEvent = null;

			for (final CascadeEvent currentEvent : cascadeData
					.getCascadeEvents(cascadeIndex)) {
				if (prevEvent != null) {
					if (counters.A_c_u_k[currentEvent.node][k] < 0)
						throw new RuntimeException(
								"A_c_u_k contains a negative value");

					if (counters.A_c_u_k[currentEvent.node][k] == 0) {
						logLikelihoodEvents[k] = Double.NEGATIVE_INFINITY;
						continue loop;
					}

					double acukDiff = currentEvent.timestamp
							* counters.A_c_u_k[currentEvent.node][k]
							- counters.tilde_A_c_u_k[currentEvent.node][k];

					if (acukDiff < 0)
						if (acukDiff < SurvivalFactorizationEM_ModelCounters.NEGATIVE_TOLERANCE)
							throw new RuntimeException("Inconsistent value");
						else
							acukDiff = 0;

					logLikelihoodEvents[k] += Math
							.log(counters.A_c_u_k[currentEvent.node][k])
							- model.S[currentEvent.node][k] * acukDiff;
				}

				prevEvent = currentEvent;
			}
		}

		return logLikelihoodEvents;
	}

	/*
	 * E-Step: updates gamma_{c,k}
	 */
	private double[][] E_Step(CascadeData cascadeData,
			SurvivalFactorizationEM_Model model) {

		final double gammaNew[][] = new double[cascadeData.n_cascades][model.nFactors];

		for (int c = 0; c < cascadeData.n_cascades; c++) {
			System.arraycopy(llk_exponentialTerms[c], 0, gammaNew[c], 0,
					model.nFactors);

			ExpLogUtils.expNormalization(gammaNew[c]);

			for (int k = 0; k < model.nFactors; k++) {
				if (gammaNew[c][k] <= 1E-100)
					gammaNew[c][k] = 0;

				if (gammaNew[c][k] >= 1 - 1E-100)
					gammaNew[c][k] = 1;
			}
		}

		return gammaNew;
	}

	private SurvivalFactorizationEM_Model iterateEM(CascadeData cascadeData,
			SurvivalFactorizationEM_Model model, int nMaxIterations)
					throws Exception {

		System.out.println("Learning phase: starting at " + new Date());

		int saveIteration = 0;
		final long initTime = System.currentTimeMillis();

		final int nVertices = cascadeData.getNNodes();
		final int nCascades = cascadeData.getNCascades();
		final int nFactors = model.getNFactors();

		double[][] gamma = new double[nCascades][nFactors];

		final SurvivalFactorizationEM_ModelCounters counters = new SurvivalFactorizationEM_ModelCounters(
				nVertices, nFactors);

		for (int c = 0; c < nCascades; c++) {
			for (int k = 0; k < nFactors; k++)
				gamma[c][k] = randomGen.nextDouble();

			Weka_Utils.normalize(gamma[c]);
		}

		double improvement = 1;
		int iterationsDone = 1;

		double prevlogLikelihood = computeLogLikelihood(cascadeData, model,
				counters, gamma);

		System.out.format("Init loglikelihood\t\t%.8f\r\n", prevlogLikelihood);

		do {// for each iteration
			gamma = E_Step(cascadeData, model);
			model = M_Step(cascadeData, model, gamma, counters);
			final double logLikelihood = computeLogLikelihood(cascadeData,
					model, counters, gamma);

			saveIteration++;

			if (saveIteration == save_step) {
				saveIteration = 0;

				if (prevlogLikelihood > logLikelihood) {
					improvement = (logLikelihood - prevlogLikelihood)
							/ prevlogLikelihood;

					System.out
					.format("######\tIteration: %d\tLogLikelihood\t%.5f (Loss: %.6f)\r\n",
							iterationsDone, logLikelihood, improvement);
					System.out.println(String.format("\tCurrent Likelihood:"
							+ "\t%.5f\n\tPrevious Likelihood:\t%.5f",
							logLikelihood, prevlogLikelihood));
				} else {
					improvement = (prevlogLikelihood - logLikelihood)
							/ prevlogLikelihood;

					System.out
					.format("Iteration: %d\tLogLikelihood\t%.5f (Improvement: %.6f)\r\n",
							iterationsDone, logLikelihood, improvement);
				}

				prevlogLikelihood = logLikelihood;

				if (storeTemporaryModels)
					try {
						model.store("tmp_" + iterationsDone + ".model");
					} catch (final Exception e) {
						e.printStackTrace();
					}
			}

			++iterationsDone;
		} while (iterationsDone <= nMaxIterations
				&& Math.abs(improvement) > eps);

		final long learningTime = System.currentTimeMillis() - initTime;
		System.out.format("Learning Phase: DONE  (%d iterations, %.0f secs)\n",
				iterationsDone, (double) learningTime / 1000);

		return model;
	}

	/*
	 * M-Step updates the model
	 */
	private SurvivalFactorizationEM_Model M_Step(CascadeData cascadeData,
			SurvivalFactorizationEM_Model model, double[][] gamma,
			SurvivalFactorizationEM_ModelCounters counters) throws Exception {

		resetUtilityVariables();
		counters.cumulateS(model);

		final Set<Integer> inactiveVertices = new HashSet<Integer>();

		for (int c = 0; c < cascadeData.getNCascades(); c++) {
			counters.updateCountersOnCascade(cascadeData, c, model);

			for (int k = 0; k < model.nFactors; k++)
				pi_new[k] += gamma[c][k];

			final List<CascadeEvent> eventsCurrCascade = cascadeData
					.getCascadeEvents(c);
			inactiveVertices.clear();
			inactiveVertices.addAll(cascadeData.getNodeIds());

			CascadeEvent prevEvent = null;

			for (final CascadeEvent currentEvent : eventsCurrCascade) {
				for (int k = 0; k < model.nFactors; k++) {
					if (prevEvent != null) {
						S_new_num[currentEvent.node][k] += gamma[c][k];

						double acukDiff = currentEvent.timestamp
								* counters.A_c_u_k[currentEvent.node][k]
										- counters.tilde_A_c_u_k[currentEvent.node][k];

						if (acukDiff < 0)
							if (acukDiff < SurvivalFactorizationEM_ModelCounters.NEGATIVE_TOLERANCE)
								throw new RuntimeException("Inconsistent value");
							else
								acukDiff = 0;

						S_new_den[currentEvent.node][k] += gamma[c][k]
								* acukDiff;
					}

					double prod = gamma[c][k]
							* counters.R_c_u_k[currentEvent.node][k];

					if (Double.isNaN(prod) || Double.isInfinite(prod))
						prod = 0;

					A_new_num[currentEvent.node][k] += model.A[currentEvent.node][k]
							* prod;

					A_new_den[currentEvent.node][k] += gamma[c][k]
							* (counters.tilde_S_c_k[k]
									- counters.tilde_S_c_u_k[currentEvent.node][k]
											- currentEvent.timestamp
											* (counters.S_c_k[k] - counters.S_c_u_k[currentEvent.node][k])
											+ cascadeData.t_max
											* (counters.S_k[k] - counters.S_c_k[k]) - currentEvent.timestamp
											* (counters.S_k[k] - counters.S_c_k[k]));
				}

				inactiveVertices.remove(currentEvent.node);
				prevEvent = currentEvent;
			}

			// for each inactive node (update S_den)
			for (final int inactiveNode : inactiveVertices)
				for (int k = 0; k < model.nFactors; k++) {
					double acukDiff = cascadeData.t_max
							* counters.A_c_u_k[inactiveNode][k]
									- counters.tilde_A_c_u_k[inactiveNode][k];

					if (acukDiff < 0)
						if (acukDiff < SurvivalFactorizationEM_ModelCounters.NEGATIVE_TOLERANCE)
							throw new RuntimeException("Inconsistent value");
						else
							acukDiff = 0;

					S_new_den[inactiveNode][k] += gamma[c][k] * acukDiff;
				}

			// update Phi_new
			if (enableContentLikelihood) {
				final List<WordOccurrence> contentCurrCascade = cascadeData
						.getCascadeContent(c);

				for (final WordOccurrence wo : contentCurrCascade)
					for (int k = 0; k < model.nFactors; ++k)
						Phi_new_num[wo.word][k] += gamma[c][k] * wo.cnt;

				for (int k = 0; k < model.nFactors; ++k)
					Phi_new_den[k] += gamma[c][k]
							* cascadeData.getLenghtOfCascadeContent(c);
			}
		}

		// now compute the new model parameters

		for (int k = 0; k < model.nFactors; k++) {
			for (int u = 0; u < cascadeData.n_nodes; u++) {
				S_new[u][k] = S_new_num[u][k]
						/ (exponentialPrior_rate + S_new_den[u][k]);

				if (S_new[u][k] == 0)
					S_new[u][k] = 1e-100;

				A_new[u][k] = A_new_num[u][k]
						/ (exponentialPrior_rate + A_new_den[u][k]);

				if (A_new[u][k] == 0)
					A_new[u][k] = 1e-100;
			}

			if (enableContentLikelihood)
				for (int w = 0; w < cascadeData.n_words; w++)
					Phi_new[w][k] = (priorGamma_scale - 1 + Phi_new_num[w][k])
							/ (priorGamma_rate + Phi_new_den[k]);
		}

		Weka_Utils.normalize(pi_new);

		// build the new model
		final SurvivalFactorizationEM_Model updatedModel = new SurvivalFactorizationEM_Model(
				cascadeData.n_nodes, cascadeData.n_words, model.nFactors);

		updatedModel.setA(A_new);
		updatedModel.setS(S_new);
		updatedModel.setPhi(Phi_new);
		updatedModel.setPi(pi_new);

		return updatedModel;
	}

	private void resetUtilityVariables() {
		resetMatrix(A_new_num);
		resetMatrix(A_new_den);

		resetMatrix(S_new_num);
		resetMatrix(S_new_den);

		resetMatrix(Phi_new_num);
		Arrays.fill(Phi_new_den, 0);

		Arrays.fill(pi_new, 0);
	}
}