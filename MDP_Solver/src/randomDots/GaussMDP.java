package randomDots;

import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;

import umontreal.iro.lecuyer.probdist.NormalDist;
import MDPCore.SparseMDP;

public class GaussMDP {
	// [R_S R_P R_N] = [-0.1 100 x0] fits data the best
	// [R_S R_P R_N] = [-1 50 -100] gives the longer RT for error choices
	static double R_S = -0.1;
	static double R_P = 1000; // 450 for Ssp_data
	static double R_N = 0;

	static double maxO = 5;
	static double minO = -5;// o/sigma ~ N(\mu,1) where mu\in[-1,1]
	static double sigma = 1;
	static int sparseness = 101;// number of possible values of observation o_t
	static int numSamples = 301;
	static int numActions = 3;
	static int numStates = (sparseness - 1)* numSamples * (numSamples + 1) / 2 + numSamples;

	static double[][] buildReward() {
		double[][] L = new double[numStates][numActions];
		// Reward for sampling
		for (int i = 0; i < numStates; i++) {
			L[i][0] = R_S;
		}
		// Reward for picking right or left

		for (int i = 0; i < numStates; i++) {
			int[] muTime = getIdx2D(i);
			int t = muTime[1];
			double mu =computeMu(muTime[0],t);

			NormalDist n = new NormalDist(mu, sigma / Math.sqrt(t));

			double prob_P = 1.0 - n.cdf(0);
			L[i][1] = prob_P * R_P + (1 - prob_P) * R_N;
			L[i][2] = prob_P * R_N + (1 - prob_P) * R_P;
		}
		// for terminating states
		for (int ia = 0; ia < numActions; ia++) {
			for (int i = 0; i < numStatesAtT(numSamples); i++) {
				int s = getIdx1D(i, numSamples);
				L[s][ia] = 0;
			}
		}
		return L;
	}
	
	static int numStatesAtT(int t){
		return t * (sparseness - 1) + 1;
	}
	
	
	static double computeMu(int idx1, int idx2){
		return (maxO - minO) * idx1 / (numStatesAtT(idx2) -1) + minO;
	}
	
	static int[] getIdx2D(int i) {
		int[] muTime = new int[2];

		int t = 1;

		while (i - numStatesAtT(t) >= 0) {
			i -= numStatesAtT(t);
			t++;
		}
		muTime[0] = i;
		muTime[1] = t;
		return muTime;
	}

	static int getIdx1D(int i, int t) {
		if (i >= numStatesAtT(t)) {
			return -1;
		}
		return i + (t - 1) * t / 2 * (sparseness - 1) + t - 1;
	}

	static void sparseTran(double[][] nextP, int[][] nextState, int s) {
		int[] muTime = getIdx2D(s);
		int t = muTime[1];
		double mu =computeMu(muTime[0],t);
		//O be observation and x be hidden variable
		//P[0|x] = Normal(x, sigma^2)
		//Current distribution over Pr[x] = Normal(mu, sigma^2/t)
		//E[P[O|x]] = P[O|mu, t] = Normal(mu, sigma^2/t + sigma^2)
		NormalDist n = new NormalDist(mu, sigma * Math.sqrt(1.0 + 1.0/t));
		double sumP = 0;
		for (int i = 0; i < sparseness; i++) {
			double o = computeMu(i,1);
			nextP[s][i] = n.density(o);
			sumP += nextP[s][i];
			nextState[s][i] = getIdx1D(muTime[0] + i, t + 1);
		}
		for (int i = 0; i < sparseness; i++) {
			nextP[s][i] /= sumP;
		}

	}

	static SparseMDP buildSparseMDP() {
		int[][] nextState = new int[numStates][sparseness];
		double[][] nextP = new double[numStates][sparseness];
		// When action is keep sampling
		for (int s = 0; s < numStates; s++) {
			sparseTran(nextP, nextState, s);
		}

		boolean[] terminate = new boolean[numStates];
		for (int i = 0; i < numStatesAtT(numSamples); i++) {
			int s = getIdx1D(i, numSamples);
			for(int j = 0; j < sparseness; j++){
				nextP[s][j] = 0;
				nextState[s][j] = s;
			}
			nextP[s][0] = 1;
			terminate[s] = true;
		}

		ArrayList<int[][]> intT;
		intT = new ArrayList<int[][]>();
		intT.add(nextState);
		ArrayList<double[][]> T;
		T = new ArrayList<double[][]>();
		T.add(nextP);

		int[][] nextState2 = new int[numStates][sparseness];
		double[][] nextP2 = new double[numStates][sparseness];
		for (int i = 0; i < numStates; i++) {
			nextState2[i][0] = numStates - 1;
			nextP2[i][0] = 1;
		}
		for (int a = 1; a < numActions; a++) {
			intT.add(nextState2);
			T.add(nextP2);
		}
		double[][] L = buildReward();

		SparseMDP m = new SparseMDP(sparseness, intT, T, L, terminate,
				numStates, numActions, 1);
		return m;
	}

	static void computeBoundary(int[] Policy, String filename) {
		PrintStream fstream;
		try {
			fstream = new PrintStream("Matlab/GaussPolicy_" + filename + ".txt");

			for (int t = 1; t < numSamples; t++) {
				int nR = 0;
				int nL = 0;
				for (int j = 0; j < numStatesAtT(t); j++) {
					int s = getIdx1D(j, t);
					if (Policy[s] == 0) {
						// System.out.print("S ");
					} else {
						if (Policy[s] == 1) {
							nR++;
							// System.out.print("R ");
						} else {
							// System.out.print("L ");
							nL++;
						}
					}
				}
				fstream.printf("%g %g %d\n", computeMu(nL,t), computeMu(numStatesAtT(t) - 1 - nR ,t),t);
			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	static public void main(String[] args) {
		if (args.length >= 3) {
			R_S = Double.valueOf(args[0]);
			R_P = Double.valueOf(args[1]);
			R_N = Double.valueOf(args[2]);
		}
		if (args.length >= 4) {
			numSamples = Integer.valueOf(args[3]);
			numStates = numSamples * (numSamples + 1) / 2;
		}
		SparseMDP m = buildSparseMDP();
		m.valueIteration();
		int[] Policy = m.getPolicy();
		String rewardLable = String.format("%.1f_%.1f_%.1f", R_S, R_P, R_N);
		computeBoundary(Policy, rewardLable);
		System.out.println("Done");
	}

}
