package utils;

public class TestDistributions {
	
	public static void main(String[] args){
		int seed = 12345;
		Randoms r = new Randoms(seed);
		System.out.println(r.nextGamma(2, 2));
		
		
	}
}
