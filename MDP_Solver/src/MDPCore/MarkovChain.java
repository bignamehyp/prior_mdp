package MDPCore;

import java.util.Random;
import java.util.Vector;


public class MarkovChain {
	private int numStates;
	private double accumP[][];
	private Random generator;	
	/**
	 * Constructor for MarkovChain
	 * @param p: the transition matrix of the Markov Chain
	 * @param nS: number of states
	 */
	public MarkovChain(double [][] p, int nS){
		numStates = nS;
		accumP = new double[nS][nS];
		for(int i = 0; i < nS; i++){
			accumP[i][0] = p[i][0];
			for(int j = 1; j < nS; j++){
				accumP[i][j] = accumP[i][j-1] + p[i][j];
			}
		}
		generator = new Random(); //Two Random objects created within the same millisecond will have the same sequence of random numbers.
	}
	
	/**
	 * Simulate a Markov chain trajectory 
	 * @param initS: initial state
	 * @param length: the length of the output trajectory
	 * @return an array of state trajectory with initial state = initS
	 */
	public int[] trajectory(int initS, int length){
		int [] path = new int[length];
		path[0] = initS;
		for(int t = 1; t < length; t++){
			path[t] = nextState(path[t-1]);
		}
		return path;
	}
	
	/**
	 * Simulate a Markov Chain trajectory until it hits a termination state
	 * @param initS: initial state
	 * @param isTerminate: a binary vector indicating which state is terminal
	 * @return: a vector of state trajectory starting from initS ending at one of the termination state
	 */
	public Vector<Integer> trajectory(int initS, boolean [] isTerminate){
		Vector<Integer> path = new Vector<Integer>();
		path.add(initS);
		do{
			path.add(nextState(path.lastElement()));
			
		}while(isTerminate[path.lastElement()]==false);
		return path;
	}
	
	
	
	/**
	 * Generate the next state determined on the transition probability matrix
	 * @param initS: current state
	 * @return nextState 
	 */
	public int nextState(int initS){
		int next = numStates/2;
		int l = 0;
		int r = numStates - 1;
		//Binary Search
		double p = generator.nextDouble();
		while(l <= r){//we want accumP[initS][l - 1] < p <= accumP[initS][r + 1]
			next = (l + r)/2;
			if(p < accumP[initS][next] )
				r = next - 1;
			else if( p > accumP[initS][next])
				l = next + 1;
			else if(l == next)
				return next;
			else
				r = next;
		}
		return l;
	}	
	
	/**
	 * Generate the next state determined on the transition probability matrix
	 * @param initS: current state
	 * @param nS, number of states
	 * @param t, transition matrix with size nS by nS 
	 * @return nextState 
	 */
	public int nextState(int initS, int nS, double[][] t){
		double [] AP = new double[nS];
		AP[0] = t[initS][0];
		for(int j = 1; j < nS; j++){
			AP[j] = AP[j-1] + t[initS][j];
		}
		int next = nS/2;
		int l = 0;
		int r = nS - 1;
		//Binary Search
		double p = generator.nextDouble();
		while(l <= r){//we want accumP[initS][l - 1] < p <= accumP[initS][r + 1]
			next = (l + r)/2;
			if(p < AP[next] )
				r = next - 1;
			else if( p > AP[next])
				l = next + 1;
			else if(l == next)
				return next;
			else
				r = next;
		}
		return l;
	}
}
