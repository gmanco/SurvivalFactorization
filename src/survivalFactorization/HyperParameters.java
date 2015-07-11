package survivalFactorization;

import java.io.Serializable;

public class HyperParameters implements Serializable {
	
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
