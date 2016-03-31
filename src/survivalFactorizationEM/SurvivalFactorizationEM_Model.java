package survivalFactorizationEM;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import utils.Weka_Utils;

public class SurvivalFactorizationEM_Model implements Serializable {

	private static final long serialVersionUID = -19283648276165L;

	private static void copyMatrix(double[][] src, double[][] dest) {
		final int n = src.length;
		final int m = n > 0 ? src[0].length : -1;

		if (n != dest.length || n > 0 && m != dest[0].length)
			throw new RuntimeException("Unable to copy");

		for (int i = 0; i < n; ++i)
			System.arraycopy(src[i], 0, dest[i], 0, m);
	}

	public static SurvivalFactorizationEM_Model readFromFile(String filename)
			throws Exception {

		ObjectInputStream in = null;
		SurvivalFactorizationEM_Model o = null;

		try {
			in = new ObjectInputStream(new GZIPInputStream(new FileInputStream(
					new File(filename))));
			o = (SurvivalFactorizationEM_Model) in.readObject();
			in.close();
			return o;
		} catch (final Exception ex) {
			try {
				ex.printStackTrace();
				in.close();
			} catch (final Exception ex2) {
				ex2.printStackTrace();
			}

			throw ex;
		}
	}

	public int nVertices;
	public int nWords;
	public int nFactors;

	public double A[][];
	public double S[][];

	public double Phi[][];

	public double pi[];

	public SurvivalFactorizationEM_Model(int nVertices, int nWords, int nFactors) {
		this.nVertices = nVertices;
		this.nWords = nWords;
		this.nFactors = nFactors;

		init();
	}

	public double[][] getA() {
		return A;
	}

	public int getNFactors() {
		return nFactors;
	}

	public int getNVertices() {
		return nVertices;
	}

	public int getNWords() {
		return nWords;
	}

	public double[][] getPhi() {
		return Phi;
	}

	public double[] getPi() {
		return pi;
	}

	public double[][] getS() {
		return S;
	}

	private void init() {
		A = new double[nVertices][nFactors];
		S = new double[nVertices][nFactors];

		Phi = new double[nWords][nFactors];
		pi = new double[nFactors];

		final double meanExpVertices = 1d / nVertices;
		final double meanExpWords = 1d / (nFactors * nFactors + 1d);

		for (int k = 0; k < nFactors; ++k) {
			for (int u = 0; u < nVertices; ++u) {
				A[u][k] = SurvivalFactorizationEM_Configuration.randomGen
						.nextExp(meanExpVertices);
				S[u][k] = SurvivalFactorizationEM_Configuration.randomGen
						.nextExp(meanExpVertices);
			}

			for (int w = 0; w < nWords; ++w)
				Phi[w][k] = SurvivalFactorizationEM_Configuration.randomGen
				.nextExp(meanExpWords);

			pi[k] = SurvivalFactorizationEM_Configuration.randomGen
					.nextDouble();
		}

		Weka_Utils.normalize(pi);
	}

	public void setA(double[][] a) {
		copyMatrix(a, A);
	}

	public void setNFactors(int nFactors) {
		this.nFactors = nFactors;
	}

	public void setNVertices(int nVertices) {
		this.nVertices = nVertices;
	}

	public void setNWords(int nWords) {
		this.nWords = nWords;
	}

	public void setPhi(double[][] phi) {
		copyMatrix(phi, Phi);
	}

	public void setPi(double[] pi) {
		System.arraycopy(pi, 0, this.pi, 0, pi.length);
	}

	public void setS(double[][] s) {
		copyMatrix(s, S);
	}

	public void store(String fileName) throws Exception {
		ObjectOutputStream out = null;

		try {
			final GZIPOutputStream gz = new GZIPOutputStream(
					new FileOutputStream(new File(fileName)));

			out = new ObjectOutputStream(gz);
			out.writeObject(this);

			out.flush();
			gz.flush();
			gz.finish();
			out.close();
		} catch (final Exception ex) {
			try {
				ex.printStackTrace();
				out.close();
			} catch (final Exception ex2) {
				ex2.printStackTrace();
			}

			throw ex;
		}
	}
}