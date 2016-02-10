package survivalFactorizationEM;

import java.util.HashSet;
import java.util.List;
import java.util.ListIterator;
import java.util.Set;

import data.CascadeData;
import data.CascadeEvent;

public class SurvivalFactorizationEM_ModelCounters {
	double[][] A_c_u_k;
	double[][] tilde_A_c_u_k;
	double[] A_c_k;
	double[] tilde_A_c_k;

	double[][] R_c_u_k;

	double[] S_k;
	double[] S_c_k;
	double[] tilde_S_c_k;
	double[][] S_c_u_k;
	double[][] tilde_S_c_u_k;

	double[] L_c_k;

	public int nVertices;
	public int nFactors;

	public SurvivalFactorizationEM_ModelCounters(int n, int k) {
		this.nFactors = k;
		this.nVertices = n;

		resetCounters();
		S_k = new double[nFactors];
	}

	/*
	 * public void update(CascadeData cascadeData, SurvivalFactorizationEM_Model
	 * model) { resetCounters(); List<CascadeEvent> eventsCurrCascade;
	 * Set<Integer> inactiveVertices = new HashSet<Integer>();
	 * 
	 * for (int k = 0; k < nFactors; k++) { for (int n = 0; n < nVertices; n++)
	 * S_k[k] = model.S[n][k]; }
	 * 
	 * for (int c = 0; c < cascadeData.getNCascades(); c++) {
	 * 
	 * eventsCurrCascade = cascadeData.getCascadeEvents(c);
	 * inactiveVertices.clear(); //FIXME: mi sembra estremamente inefficiente.
	 * inactiveVertices.addAll(cascadeData.getNodeIds()); CascadeEvent prevEvent
	 * = null; for (CascadeEvent currentEvent : eventsCurrCascade) { for (int k
	 * = 0; k < nFactors; k++) { int n = currentEvent.node; double time =
	 * currentEvent.timestamp;
	 * 
	 * if (prevEvent != null) { A_c_u_k[c][n][k] += model.A[prevEvent.node][k] +
	 * A_c_u_k[c][prevEvent.node][k]; tilde_A_c_u_k[c][n][k] += time *
	 * model.A[prevEvent.node][k] + tilde_A_c_u_k[c][prevEvent.node][k]; }
	 * 
	 * tilde_A_c_k[c][k] += time * model.A[n][k]; A_c_k[c][k] += model.A[n][k];
	 * 
	 * //FIXME: skip first active node S_c_k[c][k] += model.S[n][k];
	 * tilde_S_c_k[c][k] += time * model.S[n][k]; L_c_k[c][k] +=
	 * Math.log(model.S[n][k]);
	 * 
	 * inactiveVertices.remove(currentEvent.node); prevEvent = currentEvent; } }
	 * // for each inactive node (update S_den) for (int inactiveNode :
	 * inactiveVertices) { for (int k = 0; k < model.nFactors; k++) {
	 * tilde_S_c_k[c][k] += cascadeData.t_max * model.S[inactiveNode][k]; } } }
	 * 
	 * for (int c = 0; c < cascadeData.getNCascades(); c++) {
	 * 
	 * eventsCurrCascade = cascadeData.getCascadeEvents(c);
	 * ListIterator<CascadeEvent> li =
	 * eventsCurrCascade.listIterator(eventsCurrCascade.size()); CascadeEvent
	 * prevEvent = null;
	 * 
	 * while (li.hasPrevious()) { CascadeEvent currentEvent = li.previous(); if
	 * (prevEvent != null) for (int k = 0; k < nFactors; k++) {
	 * if(A_c_u_k[c][prevEvent.node][k]==0.0) throw new RuntimeException();
	 * R_c_u_k[c][currentEvent.node][k] += 1 / (A_c_u_k[c][prevEvent.node][k]) +
	 * R_c_u_k[c][prevEvent.node][k]; } prevEvent = currentEvent; } } }
	 */

	public void updateCountersOnCascade(CascadeData cascadeData, int cascadeIndex,
			SurvivalFactorizationEM_Model model) {
		resetCounters();
		List<CascadeEvent> eventsCurrCascade;

		eventsCurrCascade = cascadeData.getCascadeEvents(cascadeIndex);
		CascadeEvent prevEvent = null;
		for (CascadeEvent currentEvent : eventsCurrCascade) {
			int n = currentEvent.node;
			double time = currentEvent.timestamp;

			for (int k = 0; k < nFactors; k++) {
				if (prevEvent != null) {

					A_c_u_k[n][k] = model.A[prevEvent.node][k] + A_c_u_k[prevEvent.node][k];
					tilde_A_c_u_k[n][k] = time * model.A[prevEvent.node][k] + tilde_A_c_u_k[prevEvent.node][k];

					S_c_u_k[n][k] = S_c_u_k[prevEvent.node][k];
					tilde_S_c_u_k[n][k] = tilde_S_c_u_k[prevEvent.node][k];

					if (model.S[n][k] <= 0) {
						throw new RuntimeException("S contains a negative value");
					}
					L_c_k[k] += Math.log(model.S[n][k]);
				}

				S_c_u_k[n][k] += model.S[n][k];
				tilde_S_c_u_k[n][k] += time * model.S[n][k];

				tilde_A_c_k[k] += time * model.A[n][k];
				A_c_k[k] += model.A[n][k];

				S_c_k[k] += model.S[n][k];
				tilde_S_c_k[k] += time * model.S[n][k];

			}
			prevEvent = currentEvent;
		}

		ListIterator<CascadeEvent> li = eventsCurrCascade.listIterator(eventsCurrCascade.size());
		prevEvent = null;

		while (li.hasPrevious()) {
			CascadeEvent currentEvent = li.previous();
			if (prevEvent != null)
				for (int k = 0; k < nFactors; k++) {
					if (A_c_u_k[prevEvent.node][k] == 0.0)
						throw new RuntimeException();
					R_c_u_k[currentEvent.node][k] += 1 / (A_c_u_k[prevEvent.node][k]) + R_c_u_k[prevEvent.node][k];
				}
			prevEvent = currentEvent;
		}
	}


	private void resetCounters() {
		A_c_u_k = new double[nVertices][nFactors];
		tilde_A_c_u_k = new double[nVertices][nFactors];
		A_c_k = new double[nFactors];
		tilde_A_c_k = new double[nFactors];

		R_c_u_k = new double[nVertices][nFactors];
		S_c_u_k = new double[nVertices][nFactors];
		tilde_S_c_u_k = new double[nVertices][nFactors];
		S_c_k = new double[nFactors];
		tilde_S_c_k = new double[nFactors];
		L_c_k = new double[nFactors];

	}

	public void cumulateS(SurvivalFactorizationEM_Model model) {
		S_k = new double[nFactors];

		for (int k = 0; k < nFactors; k++) {
			for (int n = 0; n < nVertices; n++)
				S_k[k] += model.S[n][k];
		}

	}
}
