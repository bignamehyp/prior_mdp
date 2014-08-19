package test;

import java.util.ArrayList;

import Jama.Matrix;
import MDPCore.MDP;

public class GridWorld {
	static int h = 4;
	static int w = 4;

	static int idx2D(int i, int j){
		if(i < 0)
			i = 0;
		if(j < 0)
			j = 0;
		if(i >= h)
			i = h -1;
		if(j >= w )
			j = w -1;
		return j + i * w;
	}
	
	
	static MDP buildMDP(){
		int numStates = h * w;
		int numActions = 4;

		double [][] L = new double[numStates][numActions];
		ArrayList<double[][]> T = new ArrayList<double[][]>();
		boolean [] terminate = new boolean[numStates];
		terminate[idx2D(0,0)] = true;
		terminate[idx2D(h-1,w-1)] = true;

		for(int u = 0; u < numActions; u++){
			double [][] P = new double[numStates][numStates];
			for(int i = 0; i < h; i++){	
				for(int j = 0; j < w; j++){
					if (terminate[idx2D(i,j)] == false){
						L[idx2D(i,j)][u] = -1;
						switch(u){	//To use policy iteration, each policy needs to be proper
						case 0://Moving Right
							P[idx2D(i,j)][idx2D(i,j+1)] += 0.7;	
							P[idx2D(i,j)][idx2D(i,j-1)] += 0.1;
							P[idx2D(i,j)][idx2D(i+1,j)] += 0.1;
							P[idx2D(i,j)][idx2D(i-1,j)] += 0.1;
							break;
						case 1://Moving Left
							P[idx2D(i,j)][idx2D(i,j-1)] += 0.7;
							P[idx2D(i,j)][idx2D(i,j+1)] += 0.1;
							P[idx2D(i,j)][idx2D(i+1,j)] += 0.1;
							P[idx2D(i,j)][idx2D(i-1,j)] += 0.1;
							break;
						case 2://Moving Down
							P[idx2D(i,j)][idx2D(i+1,j)] += 0.55;
							P[idx2D(i,j)][idx2D(i,j-1)] += 0.15;
							P[idx2D(i,j)][idx2D(i,j+1)] += 0.15;
							P[idx2D(i,j)][idx2D(i-1,j)] += 0.15;
							break;
						case 3://Moving Up
							P[idx2D(i,j)][idx2D(i-1,j)] += 0.55;
							P[idx2D(i,j)][idx2D(i+1,j)] += 0.15;
							P[idx2D(i,j)][idx2D(i,j+1)] += 0.15;
							P[idx2D(i,j)][idx2D(i,j-1)] += 0.15;
							break;
						}	
					}
					else{
						P[idx2D(i,j)][idx2D(i,j)] = 1;
						L[idx2D(i,j)][u] = 0;
					}
				}
			}
			T.add(P);
		}		
		MDP m = new MDP(T, L, terminate, numStates, numActions, 1);
		return m;
	}
	
	static void printValue2D(Matrix v){
		for(int i = 0; i < h; i++){
			for(int j = 0; j < w; j++){
				System.out.printf("%.2f\t",v.get(idx2D(i,j), 0));
			}
			System.out.println();
		}
	}
	

	static void printPolicy2D(int [] v){
		for(int i = 0; i < h; i++){
			for(int j = 0; j < w; j++){
				switch(v[idx2D(i,j)]){
				case 0:
					System.out.print("Right\t");
					break;
				case 1:
					System.out.print("Left\t");
					break;
				case 2:
					System.out.print("Down\t");
					break;	
				case 3:
					System.out.print("Up\t");
					break;
				}
			}
			System.out.println();
		}
	}
	
	static public void main(String[] args){
		double lambda = 0.5;
		double epislon = 0.1;
		int maxIters = 10000;
		if(args.length >= 2){
			w = Integer.valueOf(args[0]);
			h = Integer.valueOf(args[1]);
		}
		if(args.length >= 4){
			lambda = Double.valueOf(args[2]);
			maxIters = Integer.valueOf(args[3]);
		}
		if(args.length >= 5)
			epislon = Double.valueOf(args[4]);
		
		long start = System.currentTimeMillis();
		MDP m = buildMDP();
		m.valueIteration();
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("Value iteration takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		
		if(w <= 10){
			Matrix value = m.getValue();
			printValue2D(value);
			int [] policy = m.getPolicy();
			printPolicy2D(policy);
		}
		
		m = buildMDP();
		start = System.currentTimeMillis();
		m.policyIteration();
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("Policy iteration takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		if(w <= 10){
			Matrix value = m.getValue();
			printValue2D(value);
			int [] policy = m.getPolicy();
			printPolicy2D(policy);
		}

		m = buildMDP();
		start = System.currentTimeMillis();		
		m.SARSA(lambda, maxIters,epislon);
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("SARSA takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		if(w <= 10){
			Matrix value = m.getValue();
			printValue2D(value);
			int [] policy = m.getPolicy();
			printPolicy2D(policy);
		}
		
		m = buildMDP();
		start = System.currentTimeMillis();		
		m.QLearning(lambda,maxIters,epislon);
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("Q Learning takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		if(w <= 10){
			Matrix value = m.getValue();
			printValue2D(value);
			int [] policy = m.getPolicy();
			printPolicy2D(policy);
		}
	}
}
