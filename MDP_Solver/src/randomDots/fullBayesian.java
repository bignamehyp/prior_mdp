package randomDots;

import java.util.ArrayList;

import umontreal.iro.lecuyer.probdist.BetaDist;
import MDPCore.MDP;
import MDPCore.SparseMDP;

public class fullBayesian {
	// [R_S R_P R_N] = [-0.1 100 -20] fits data the best
	// [R_S R_P R_N] = [-1 50 -100]  gives the longer RT for error choices
	static double R_S = -0.2; 
	static double R_P = 100; //450  for Ssp_data
	static double R_N = 0;
	//static double prior = 0.8;
	static double alpha_init = 1;
	static double beta_init = 1;
	//static double w_prior = 0.2;
	static int sparseness = 2;
	static int numSamples = 500;
	static int numActions = 3;
	static int numStates = numSamples * (numSamples + 1) / 2;
	
	static double [][] buildReward(){
		double [][] L = new double[numStates][numActions];
		//Reward for sampling
		for(int i = 0; i < numStates ; i++){
			L[i][0] = R_S;
		}
		//Reward for picking right or left

		for(int i = 0; i < numStates ; i++){
			int [] idx2D = Utilities.getIdx2D(i);
			BetaDist b = new BetaDist(idx2D[1] + alpha_init, 
					idx2D[0] +  beta_init);
			double prob_P  = 1.0 - b.cdf(0.5);
			L[i][1] = prob_P * R_P + (1 - prob_P) * R_N;
			L[i][2] = prob_P * R_N + (1 - prob_P) * R_P;
			}				
		for(int ia = 0; ia < numActions; ia++)
			for(int i = 0; i < numSamples; i++){
				int s = Utilities.getIdx1D(numSamples - i - 1, i);
					L[s][ia] = 0;
			}
		return L;
	}
	
	static MDP buildMDP(){
		double [][] P = new double[numStates][numStates];
		//When action is keep sampling
		for(int i = 0; i < numSamples - 1; i++){ 
			for(int j = 0; j < numSamples - 1 - i ; j++){
				int s = Utilities.getIdx1D(i,j);
	        //Probability that the next belief is b2R(next observation is O_R)?
				double P_b2R = ( j + alpha_init ) / (i + j  + alpha_init + beta_init);
			//Probability that the next belief is b2L(next observation is O_L)?
				double P_b2L = 1 - P_b2R;    
				int b2R = Utilities.getIdx1D(i , j + 1);
				int b2L  = Utilities.getIdx1D(i + 1, j);
				P[s][b2R] = P_b2R;
				P[s][b2L] = P_b2L;
			}
		}
		
		boolean [] terminate = new boolean[numStates];
		
		for(int i = 0; i < numStates; i++){	
			terminate[i] = false;
		}
		
		for(int i = 0; i < numSamples; i++){
			int s = Utilities.getIdx1D(numSamples - i - 1, i);
			P[s][s] = 1;
			terminate[s] = true;			
		}
			
		ArrayList<double [][]> T;
		T = new ArrayList<double [][]>();
		T.add(P);
		
		double [][] P2 = new double[numStates][numStates];
		for(int i = 0; i < numStates; i++){	
			for(int j = 0; j < numStates - 1; j ++)
				P2[i][j] = 0;
			P2[i][numStates-1] = 1;
		}
		T.add(P2); // when the action it Right
		T.add(P2); //when the action is Left;
		
		
		double [][] L = buildReward();

		MDP m = new MDP(T,L,terminate,numStates,numActions,1);
		return m;
	}
	
	static SparseMDP buildSparseMDP(){
		int [][] nextState = new int [numStates][sparseness];
		double [][] nextP = new double[numStates][sparseness];
		//When action is keep sampling
		for(int i = 0; i < numSamples - 1; i++){ 
			for(int j = 0; j < numSamples - 1 - i ; j++){
				int s = Utilities.getIdx1D(i,j);
	        //Probability that the next belief is b2R(next observation is O_R)?
				nextP[s][0] = ( j + alpha_init) / (i + j  + alpha_init + beta_init);
			//Probability that the next belief is b2L(next observation is O_L)?
				nextP[s][1]  = 1 - nextP[s][0];    
				nextState[s][0] = Utilities.getIdx1D(i, j + 1);
				nextState[s][1] = Utilities.getIdx1D(i + 1, j);
			}
		}
		boolean [] terminate = new boolean[numStates];
		for(int i = 0; i < numSamples; i++){
				int s = Utilities.getIdx1D(numSamples - i - 1, i);
				nextP[s][0] = 1;
				nextState[s][0] = s;
				terminate[s] = true;
		}
		
		ArrayList<int [][]> intT;
		intT = new ArrayList<int [][]>();
		intT.add(nextState);		
		ArrayList<double [][]> T;
		T = new ArrayList<double [][]>();
		T.add(nextP);
		
		int [][] nextState2 = new int [numStates][sparseness];
		double [][] nextP2 = new double[numStates][sparseness];
		for(int i = 0; i < numStates; i++){	
			nextState2[i][0] = numStates -1;
			nextP2[i][0] = 1;
			nextP2[i][1] = 0;
		}
		intT.add(nextState2);
		intT.add(nextState2);
		T.add(nextP2); // when the action it Right
		T.add(nextP2); //when the action is Left;
		
		
		double [][] L = buildReward();

		SparseMDP m = new SparseMDP(sparseness, intT, T,L,terminate,numStates,numActions,1);
		return m;
	}
	
	
	
	static public void main(String [] args){
		if(args.length >= 3){
			R_S = Double.valueOf(args[0]);
			R_P = Double.valueOf(args[1]);
			R_N = Double.valueOf(args[2]);
		}
		if(args.length >= 4){
			numSamples = Integer.valueOf(args[3]);
			numStates = numSamples * (numSamples + 1)/2;
		}
		Utilities.setNumSamples(numSamples);
		SparseMDP m = buildSparseMDP();
		m.valueIteration();
		int [] Policy = m.getPolicy();
		String rewardLable = String.format("%.1f_%.1f_%.1f", R_S,R_P,R_N);
		Utilities.computeBoundary(Policy,rewardLable);
		//Utilities.printPolicy(Policy, rewardLable);
		//Utilities.printValue(m.getValue(), rewardLable);
		System.out.println("Done");
	}
	
}
