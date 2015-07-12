package survivalFactorization;

import java.io.Serializable;
import java.util.Arrays;

public class HyperParameters implements Serializable {
	
    /**
     * 
     */
    private static final long serialVersionUID = 1L;
    double a; 
	double b; 
	double[] C;
	double [] D;

	public HyperParameters(double a2, double b2, double[] c2, double[] d2) {
		this.a = a2;
		this.b = b2;
		this.C = c2; 
		this.D = d2;
	}

	public HyperParameters(HyperParameters h) {
        this.a=h.a;
        this.b=h.b;
        this.C=Arrays.copyOf(h.C, h.C.length);
        this.D=Arrays.copyOf(h.D, h.D.length);
    }

    public double getA() {
		return a;
	}

	public double getB() {
		return b;
	}

	public double[] getC() {
		return C;
	}

	public double[] getD() {
		return D;
	}

}//HyperParameters
