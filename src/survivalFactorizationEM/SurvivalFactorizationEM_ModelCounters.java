package survivalFactorizationEM;

import java.util.Arrays;
import java.util.List;
import java.util.ListIterator;

import data.CascadeData;
import data.CascadeEvent;

public class SurvivalFactorizationEM_ModelCounters {
	public static final double NEGATIVE_TOLERANCE = -1e-5;

	private static void resetMatrix(double[][] m) {
		for (final double[] element : m)
			Arrays.fill(element, 0);
	}

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
		nVertices = n;
		nFactors = k;

		buildCounters();

		S_k = new double[nFactors];
	}

	private void buildCounters() {
		A_c_u_k = new double[nVertices][nFactors];
		tilde_A_c_u_k = new double[nVertices][nFactors];
		A_c_k = new double[nFactors];
		tilde_A_c_k = new double[nFactors];

		R_c_u_k = new double[nVertices][nFactors];

		S_k = new double[nFactors];

		S_c_k = new double[nFactors];
		tilde_S_c_k = new double[nFactors];

		S_c_u_k = new double[nVertices][nFactors];
		tilde_S_c_u_k = new double[nVertices][nFactors];

		L_c_k = new double[nFactors];
	}

	public void cumulateS(SurvivalFactorizationEM_Model model) {
		Arrays.fill(S_k, 0);

		for (int k = 0; k < nFactors; k++)
			for (int n = 0; n < nVertices; n++)
				S_k[k] += model.S[n][k];

	}

	private void resetCounters() {
		resetMatrix(A_c_u_k);
		resetMatrix(tilde_A_c_u_k);
		Arrays.fill(A_c_k, 0);
		Arrays.fill(tilde_A_c_k, 0);

		resetMatrix(R_c_u_k);

		Arrays.fill(S_c_k, 0);
		Arrays.fill(tilde_S_c_k, 0);

		resetMatrix(S_c_u_k);
		resetMatrix(tilde_S_c_u_k);

		Arrays.fill(L_c_k, 0);
	}

	public void updateCountersOnCascade(CascadeData cascadeData,
			int cascadeIndex, SurvivalFactorizationEM_Model model) {

		resetCounters();
		final List<CascadeEvent> eventsCurrCascade = cascadeData
				.getCascadeEvents(cascadeIndex);
		CascadeEvent prevEvent = null;

		for (final CascadeEvent currentEvent : eventsCurrCascade) {
			final int n = currentEvent.node;
			final double time = currentEvent.timestamp;

			for (int k = 0; k < nFactors; k++) {
				if (prevEvent != null) {
					if (model.S[n][k] <= 0)
						L_c_k[k] = Double.NEGATIVE_INFINITY;
					else
						L_c_k[k] += Math.log(model.S[n][k]);

					A_c_u_k[n][k] = model.A[prevEvent.node][k]
							+ A_c_u_k[prevEvent.node][k];
					tilde_A_c_u_k[n][k] = prevEvent.timestamp
							* model.A[prevEvent.node][k]
							+ tilde_A_c_u_k[prevEvent.node][k];

					S_c_u_k[n][k] = S_c_u_k[prevEvent.node][k];
					tilde_S_c_u_k[n][k] = tilde_S_c_u_k[prevEvent.node][k];
				}

				S_c_u_k[n][k] += model.S[n][k];
				tilde_S_c_u_k[n][k] += time * model.S[n][k];

				tilde_A_c_k[k] += time * model.A[n][k];
				A_c_k[k] += model.A[n][k];

				S_c_k[k] += model.S[n][k];
				tilde_S_c_k[k] += time * model.S[n][k];

				if (cascadeData.t_max * A_c_u_k[n][k] < tilde_A_c_u_k[n][k])
					throw new RuntimeException("Inconsistent value");

				if (time * A_c_u_k[n][k] - tilde_A_c_u_k[n][k] < NEGATIVE_TOLERANCE)
					throw new RuntimeException("Inconsistent value: "
							+ (time * A_c_u_k[n][k] - tilde_A_c_u_k[n][k]));

				if (cascadeData.t_max * S_c_u_k[n][k] - tilde_S_c_u_k[n][k] < NEGATIVE_TOLERANCE)
					throw new RuntimeException(
							"Inconsistent value: "
									+ (cascadeData.t_max * S_c_u_k[n][k] - tilde_S_c_u_k[n][k]));
			}

			prevEvent = currentEvent;
		}

		final ListIterator<CascadeEvent> li = eventsCurrCascade
				.listIterator(eventsCurrCascade.size());
		CascadeEvent nextEvent = null;

		while (li.hasPrevious()) {
			final CascadeEvent currentEvent = li.previous();

			if (nextEvent != null)
				for (int k = 0; k < nFactors; k++)
					if (A_c_u_k[nextEvent.node][k] == 0.0) {
						if (A_c_u_k[currentEvent.node][k] != 0.0)
							throw new RuntimeException("Inconsistent value");
						else
							R_c_u_k[currentEvent.node][k] = 0;
					} else
						R_c_u_k[currentEvent.node][k] += 1
								/ A_c_u_k[nextEvent.node][k]
								+ R_c_u_k[nextEvent.node][k];

			nextEvent = currentEvent;
		}
	}
}