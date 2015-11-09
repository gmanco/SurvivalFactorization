package survivalFactorizationEM;

import java.io.Serializable;
import java.util.Arrays;

import utils.Weka_Utils;

public class SurvivalFactorizationEM_Model implements Serializable {

    /**
     * 
     */
    private static final long serialVersionUID = 1L;

    public int nVertices;
    public int nWords;
    public int nFactors;

    public double A[][];
    public double S[][];
    public double Phi[][];
    public double pi[];

    public SurvivalFactorizationEM_Model(int nVertices, int nWords,
            int nFactors) {
        this.nVertices = nVertices;
        this.nWords = nWords;
        this.nFactors = nFactors;
        init();
    }

    public int getNVertices() {
        return nVertices;
    }

    public void setNVertices(int nVertices) {
        this.nVertices = nVertices;
    }

    public int getNWords() {
        return nWords;
    }

    public void setNWords(int nWords) {
        this.nWords = nWords;
    }

    public int getNFactors() {
        return nFactors;
    }

    public void setNFactors(int nFactors) {
        this.nFactors = nFactors;
    }

    public double[][] getA() {
        return A;
    }

    public void setA(double[][] a) {
        A = a;
    }

    public double[][] getS() {
        return S;
    }

    public void setS(double[][] s) {
        S = s;
    }

    public double[][] getPhi() {
        return Phi;
    }

    public void setPhi(double[][] phi) {
        Phi = phi;
    }

    public double[] getPi() {
        return pi;
    }

    public void setPi(double[] pi) {
        this.pi = pi;
    }

    private void init() {
        this.A = new double[nVertices][nFactors];
        this.S = new double[nVertices][nFactors];
        this.Phi = new double[nWords][nFactors];
        this.pi = new double[nFactors];

        double meanExpVertices = 1d / nVertices;
        double meanExpWords = 1d / (nFactors * nFactors + 1);
        for (int k = 0; k < nFactors; k++) {
            for (int u = 0; u < nVertices; u++) {
                A[u][k] = SurvivalFactorizationEM_Configuration.randomGen
                        .nextExp(meanExpVertices);
                S[u][k] = SurvivalFactorizationEM_Configuration.randomGen
                        .nextExp(meanExpVertices);
            }
            for (int w = 0; w < nWords; w++) {
                Phi[w][k] = SurvivalFactorizationEM_Configuration.randomGen
                        .nextExp(meanExpWords);

            }
            pi[k] = SurvivalFactorizationEM_Configuration.randomGen.nextDouble();
        }
        Weka_Utils.normalize(pi);
    }// init

    
    public boolean store(String fileName){
        throw new UnsupportedOperationException();
    }
    
    
}
