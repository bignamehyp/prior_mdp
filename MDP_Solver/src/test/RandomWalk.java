package test;

import java.util.ArrayList;

import Jama.Matrix;
import MDPCore.MDP;

public class RandomWalk {
	
	static int numStates = 5;
	static int numActions = 1;
	static MDP buildMDP(){
		
		double [][] P = new double[numStates][numStates];
		double [][] L = new double[numStates][numActions];
		P[0][0] = 0.5;
		P[0][1] = 0.5;
		P[numStates - 1][numStates -1] = 1;
		for(int i = 1; i < numStates - 1; i++){
			P[i][i+1] = P[i][i-1] = 0.5;
			L[i][0] = -1;
		}
		L[0][0] = -1;
		L[numStates - 1][0] = 0;
		
		boolean [] terminate = new boolean[numStates];
		
		terminate[numStates - 1] =true;
		
		ArrayList<double [][]> T;
		T = new ArrayList<double [][]>();
		T.add(P);
		MDP m = new MDP(T,L,terminate,numStates,numActions,1);
		return m;
	}
	
	static void printValue(Matrix m, int nS){
		for(int i = 0; i < nS; i++)
			System.out.printf("%5.2f ", m.get(i,0));
		System.out.println();
	}
	
	/**
	 * Random walk on a directed Graph
	 * Let h_i be the hitting time between non-termination state i and termination state 0
	 * we have h_1 = 1 + (h_2 + 0)/2
	 * h_2 = 1 + (h_3 + h_1)/2
	 * h_N = 1 + (h_N-1 + h_N)/2
	 * The solution to the above system of equations is 
	 * h_i = i*(2N + 1 - i)
	 * where N + 1 is the total number of states in the graph
	 * @param args
	 * args[0] number of states
	 * args[1] lambda for TD(lambda)
	 * args[2] max number of teration for SARSA and QLearning
	 * args[3] epsilon for epsilon-greedy algorithm in Q learning 
	 */
	static public void main(String [] args){
		double lambda = 1;
		double epislon = 0.1;
		int maxIters = 1000;
		if(args.length >= 1)
			numStates = Integer.valueOf(args[0]);
		if(args.length >= 3){
			lambda = Double.valueOf(args[1]);
			maxIters = Integer.valueOf(args[2]);
		}
		if(args.length >= 4)
			epislon = Double.valueOf(args[3]);
			
		MDP m = buildMDP();
		long start = System.currentTimeMillis();
		m.valueIteration();
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("Value iteration takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		Matrix value =  m.getValue();
		printValue(value, numStates);
		
		m = buildMDP();
		start = System.currentTimeMillis();
		m.policyIteration();
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("Policy iteration takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		Matrix value2 =  m.getValue();
		printValue(value2, numStates);
		
		m = buildMDP();
		start = System.currentTimeMillis();
		m.SARSA(lambda, maxIters,epislon);
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("SARSA takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		Matrix value3 = m.getValue();
		printValue(value3, numStates);

		m = buildMDP();
		start = System.currentTimeMillis();
		m.QLearning(lambda,maxIters,epislon);
		System.out.println("---------------------------------------------------------------------------------------");
		System.out.println("Q Learning takes " + ( System.currentTimeMillis() - start ) * 0.001 + "s" );
		Matrix value4 = m.getValue();
		printValue(value4, numStates);		
	}

}
