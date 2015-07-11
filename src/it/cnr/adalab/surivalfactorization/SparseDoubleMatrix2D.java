package it.cnr.adalab.surivalfactorization;

import java.util.HashMap;

public class SparseDoubleMatrix2D {
	
	HashMap<Integer,Double>[] byRows;
	HashMap<Integer,Double>[] byCols;
	

	public SparseDoubleMatrix2D(int n_rows, int n_cols) {
		byRows = new HashMap[n_rows]; 
		byCols = new HashMap[n_cols];
	}

	public HashMap<Integer, Double> getColumn(int c) {
		return byCols[c];
	}
	
	public HashMap<Integer, Double> getRow(int r) {
		return byRows[r];
	}

	public void set(int row, int col, double val){
		if (byRows[row] == null)
			byRows[row] = new HashMap<Integer,Double>();
		byRows[row].put(col, val);
		if (byCols[col] == null)
			byCols[col] = new HashMap<Integer,Double>();
		byCols[col].put(row, val);
	}
	
	public double get(int row, int col){
		if (byRows[row] != null && byRows[row].containsKey(col))
			return byRows[row].get(col);
		else
			return 0;		
	}

}
