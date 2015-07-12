package utils;


public class MatrixUtilities {

	
	public static void fill(double[][]a, double val){
		for(int i=0;i<a.length;i++)
			for(int j=0;j<a[0].length;j++)
				a[i][j]=val;
	}
	
	
	public static void copy(double[][] src, double[][] dest) {
		for (int i = 0; i < src.length; i++)
			System.arraycopy(src[i], 0, dest[i], 0, dest[i].length);
	}

	public static void copy(int[][] src, int[][] dest) {
		for (int i = 0; i < src.length; i++)
			System.arraycopy(src[i], 0, dest[i], 0, dest[i].length);
	}

	public static void print(double[][] m) {
		System.out.println();
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[i].length; j++) {
				System.out.printf("%1$f", m[i][j]);
				System.out.print(" ");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	public static void print(int[][] m) {
		System.out.println();
		for (int i = 0; i < m.length; i++) {
			for (int j = 0; j < m[i].length; j++) {
				System.out.print (m[i][j]);
				System.out.print(" ");
			}
			System.out.println();
		}
		System.out.println();
	}
	
	

	public static double[][] product(double[][] a, double[][] b) {
		double prod[][] = new double[a.length][b[0].length];
		for (int i = 0; i < prod.length; i++)
			for (int j = 0; j < prod[0].length; j++)
				for (int k = 0; k < a[0].length; k++)
					prod[i][j] += a[i][k] * b[k][j];
		return prod;
	}

	public static double[][] transpose(double[][] a) {
		double[][] tmp = new double[a[0].length][a.length];

		for (int j = 0; j < a[0].length; j++)
			for (int i = 0; i < a.length; i++)
				tmp[j][i] = a[i][j];

		return tmp;
	}

	public static int[][] transpose(int[][] a) {
		int[][] tmp = new int[a[0].length][a.length];

		for (int j = 0; j < a[0].length; j++)
			for (int i = 0; i < a.length; i++)
				tmp[j][i] = a[i][j];

		return tmp;
	}
}
