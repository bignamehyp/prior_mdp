package randomDots;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;

import Jama.Matrix;

public class Utilities {
	static int numSamples;
	static void setNumSamples(int n){
		numSamples = n;
	}
	
	static int [] getIdx2D(int i){
		int [] idx2D = new int[2];
		int l = numSamples;
		idx2D[0] = 0;
		while( i  - l >= 0){
			i -= l;
			idx2D[0]++;
			l--;
		}
		idx2D[1] = i;
		return idx2D;
	}
	
	static int getIdx1D(int nb, int na){
		return nb * numSamples + na - (nb -1) * nb /2;
	}
	
	
	static void printPolicy(int [] Policy, String filename){
		PrintStream fstream;
			
			try {
				fstream = new PrintStream("Matlab/policy2D_" + filename + ".txt");
				for(int n = 1; n < numSamples - 1; n++){
				
					for(int j = 0; j <= n; j++){
						int s = getIdx1D(n- j, j);
						fstream.print(String.valueOf(Policy[s]) + " ");
					}
					for(int j = 0; j < (numSamples - n -2)/2; j++){
						fstream.print("NaN ");
					}
					for(int j = 0; j < (numSamples - n -1)/2; j++){
						fstream.print("NaN ");
					}
					fstream.println();
				}
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		

		
		static void printValue(Matrix v, String filename){
			PrintStream fstream;
			
			try {
				fstream = new PrintStream("Matlab/value_" + filename + ".txt");
			
			for(int n = 1; n < numSamples - 1; n++){
			
				for(int j = 0; j <= n; j++){
					int s = getIdx1D(n- j, j);
					fstream.printf("%.4g ",v.get(s,0));
				}
				for(int j = 0; j < (numSamples - n - 2 )/2; j++){
					fstream.print("NaN ");
				}
				for(int j = 0; j < (numSamples - n -1 )/2; j++){
					fstream.print("NaN ");
				}
				fstream.println();
			}
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		static void computeBoundary(int [] Policy, String filename){
			PrintStream fstream;
			try {
				fstream = new PrintStream("Matlab/policy_" + filename + ".txt");
			
			for(int n = 1;  n < numSamples - 1; n++){
				int nS = 0;
				int nR = 0;
				int nL = 0;
				for(int j = 0; j <= n; j++){
					int s = getIdx1D(n - j, j);
					switch(Policy[s]) {
					case 0:
						nS++;
						break;
					case 1:
						nR++;
						break;
					case 2:
						nL++;
						break;
				
					}					
				}	
				//Output the range and the percentage of sampling area with respect to alpha, with respect to n. 
				fstream.printf("%d %d %d\n", nL, n - nR, n + 1 - nS - nR - nL);
			}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
}
