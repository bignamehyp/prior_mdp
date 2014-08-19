package MDPCore;

import Jama.Matrix;

public class SparseMultiplierThread extends Thread{
	public int startIndex;
	public int endIndex;
 	public int [][] ns;
	public double [][] t;  
	public Matrix v;
	public Matrix v2;
	public int sparseness;
	
	public void run(){
		for(int i = startIndex; i < endIndex; i++){
			double prod = 0;
			for(int j = 0; j < sparseness; j++)
				prod+= t[i][j] * v.get(ns[i][j], 0);
			v2.set(i, 0, prod);				
		}
	}
}
