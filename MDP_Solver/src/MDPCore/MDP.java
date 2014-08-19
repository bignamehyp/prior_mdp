package MDPCore;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Vector;

import Jama.Matrix;

/**
 * Discrete time MDP solver 
 * @author Yanping Huang Email huangyp@cs.washington.edu
 * Notations and conventions in this code follows those in Neuro-dynamic programming
 *
 */
public class MDP {
	private int numStates;
	private int numActions;
	private ArrayList<Matrix> transition; //Transition matrix with size numStates by numStates by numActions T_u(s,s')
	private boolean[] isTerminate;
	private Matrix reward; //matrix with size numStates by numActions R^u(s)	
	private Matrix value; //vector with size numStates V^*(s)
	private int[] optPolicy; //One by numStates vector. Optimal policy function of states
	private double gamma; //Discount factor
	private double tol = 1.0e-12;
	private double alphaRate = 1.0; //the learning rate alpha(t) = (1/t)^alphaRate
	private Random generator;
	/**
	 * Constructor of MDP
	 * @param p: ArrayList of 2D transition matrix, size numStates by numStates by numActions
	 * @param r: 2D reward matrix with size numStates by numActions R^u(s)	
	 * @param t: binary array indicating which state is terminal.
	 * @param nS: number of states
	 * @param nA: number of available actions
	 * @param df: discount factor 
	 */
	public MDP(ArrayList<double[][]> p, double [][] r, boolean [] t, int nS, int nA, double df){
		numStates = nS;
		numActions = nA;
		transition = new ArrayList<Matrix>();
		for(int i = 0; i < nA; i++){
			transition.add(new Matrix(p.get(i)));	
		}
		
		isTerminate = new boolean[nS];
		for(int i = 0; i < nS; i++)
			isTerminate[i] = t[i];		
		reward = new Matrix(r);
		value = new Matrix(numStates,1,0);
		optPolicy = new int[numStates];
		gamma = df;
		generator = new Random(); //Two Random objects created within the same millisecond will have the same sequence of random numbers.
	}

	public double getTolerance(){
		return tol;
	}

	public void setTolerance(double t){
		tol = t;
	}

	private void setValue(int s, double v){
		value.set(s, 0, v);
	}
	
	public Matrix getValue(){
		return value;
	}
	
	public double getValue(int s){
		if(s < 0 || s >= numStates){
			System.out.println("Illegal query to state " + s + "  out of boundary. Reset state to 0");
			s = 0;
		}
		return value.get(s, 0);
	}

	public void printValue(){
		for(int i = 0; i < numStates; i++){
			System.out.printf("%.2f ",value.get(i, 0));
		}
		System.out.println();
	}

	public int getPolicy(int s){
		if(s < 0 || s >= numStates){
			System.out.println("Illegal query to state " + s + "  out of boundary. Reset state to 0");
			s = 0;
		}
		return optPolicy[s];
	}
	
	public int [] getPolicy(){
		return optPolicy;
	}
	
	private void setPolicy(int s, int u){
		optPolicy[s] = u;
	}
	
	public void printPolicy(){
		for(int i = 0; i < numStates; i++){
			System.out.print(optPolicy[i] + " ");
		}
		System.out.println();
	}

	/**
	 * L reward
	 * alpha discount factor
	 * P transition matrix
	 * v value function
	 * Compute H = L + alpha * P * v;
	 * for ia = 1:nu  H(:,ia) = L(:,ia) + alpha*P(:,:,ia)*v; end
	 * Update v and compute policy based on hamiltonian
	 * [value, opPolicy] = max(H,[],2);
	 * @return v, policy
	 */
	private void computeHamiltonian(){
		Matrix hamiltonian  = reward.copy(); 
		for (int ia = 0; ia < numActions; ia++){
			hamiltonian.setMatrix(0, numStates-1, ia, ia, reward.getMatrix(0, numStates - 1, ia, ia).plus(transition.get(ia).times(value.times(gamma))));
		}
		for (int ia = 0; ia < numActions; ia++)
			for (int is = 0; is < numStates; is++)	
				if(isTerminate[is]){
					hamiltonian.set(is, ia, reward.get(is, ia));
				}
		policyImprovement(hamiltonian);

	}
	
	private int getRandomState(){
		return generator.nextInt(numStates);
	}
	
	private int getRandomAction(){
		return generator.nextInt(numActions);
	}
	
	private int greedyAction(int s){
		return optPolicy[s];
	}
	
	/**
	 * epsilon-greedy action selector
	 * @param s: current state s
	 * @param epislon: choose the optimal action with probability 1 - epislon, chose random action w.p. epilson
	 * @return random action a
	 */
	private int eGreedyAction(int s, double epislon){
		int a = getPolicy(s);
		if(epislon > 0){
			double p = generator.nextDouble();								
			if(p <= epislon)
				a = getRandomAction();
		}
		return a;
	}

	private boolean policyImprovement(Matrix m){
		for(int is = 0; is < numStates; is++){
			double max = Double.NEGATIVE_INFINITY;
			int argmax = -1;
			for(int ia = 0; ia < numActions; ia++){
				if(m.get(is, ia) > max){
					max = m.get(is, ia); 
					argmax = ia;
				}				
			}
			if(argmax == -1){
				System.out.println("Unapte to perform policy improvement for state s = " +  is +  "!");
				return false;
			}
			setValue(is, max);
			setPolicy(is, argmax);
		}
		return true;
	}


	private double learningRate(int iter, double exp){
		return Math.pow(1.0/iter, exp);
	}

	private Matrix TransitionGivenPolicy(){
		Matrix PPi = new Matrix(numStates, numStates, 0);
		for (int i = 0; i < numStates; i++){
			for(int j = 0; j < numStates; j++){
				PPi.set(i, j, transition.get(getPolicy(i)).get(i, j));
			}
		}
		return PPi;
	}
	/**
	 * Value Iteration
	 * Update V_{t+1} = B V_t
	 * where B is the Bellman operator
	 * Until V_{t+1} = V_t
	 * Time complexity O(|S|^2|A|) for each iterations and unknown number of iterations until convergence
	 */
	public void valueIteration(){
		while(true){
			Matrix valueOld = value.copy();
			computeHamiltonian();						
			if(valueOld.minusEquals(value).normF() < tol){
				break;
			}
		}	
	}


	/**
	 * Policy Iteration
	 * Caution! that every policy needs to be proper.
	 * Under a proper policy each state has positive probabilities to move to terminal state 
	 * Starting from an arbitrary policy pi_0
	 * Policy evaluation: solve the solution to v_1 = B_0 v_1
	 * where B_0 the Bellman operator with fixed policy pi_0
	 * Policy improvement: find the policy pi_1 that satisfies B v_1 = B_1 v_1
	 * Iterate until v_{t+1}  = v_t
	 * Time complexity: O(|S|^3 + |S|^2|A|) for each iteration and unknown number of iterations until convergence
	 */	
	public void policyIteration(){
		//Find nonterminal states
		int numNonT = 0;
		for(int i = 0; i < numStates; i++)
			if(isTerminate[i]==false)
				numNonT++;
		int [] TStates = new int [numStates - numNonT];
		int [] nonT = new int [numNonT];
		int tempN = 0;
		int tempT = 0;
		for(int i = 0; i < numStates; i++){
			if(isTerminate[i]==false)
				nonT[tempN++] = i;
			else
				TStates[tempT++] = i;
		}
		while(true){
			Matrix valueOld = value.copy();
			//Policy evaluation

			Matrix PPi = TransitionGivenPolicy();
			Matrix LPi = new Matrix(numStates, 1, 0);
			for (int i = 0; i < numStates; i++){				
				LPi.set(i, 0, reward.get(i, getPolicy(i)));					
			}

			if (gamma < 1){//discounted problem
				value = Matrix.identity(numStates, numStates).minus(PPi.times(gamma)).solve(LPi);				
			}
			else{//non-discounted problem with terminal states.
				//Solve the equation V = (I - P_NN)^-1 (L_N + P_NT * LT)
				Matrix IPNN =  Matrix.identity(numNonT, numNonT).minus(PPi.getMatrix(nonT, nonT)); // I - P_NN
				if(IPNN.rank() < numNonT){
					System.out.println("I - P_NN is singular, exist a nonproper policy, fail to perform policy iteration");
					break;
				}
				Matrix LPNT = PPi.getMatrix(nonT, TStates).times(LPi.getMatrix(TStates, 0,0)).plus(LPi.getMatrix(nonT,0,0));//L_N + P_NT * L_T 
				value.setMatrix(nonT, 0,0,  IPNN.solve(LPNT));			
				for(int i = 0; i < numStates; i++){
					if(isTerminate[i])
						setValue(i, reward.get(i,0));
						
				}
			}			
			//Policy Improvement
			computeHamiltonian();
			if(valueOld.minusEquals(value).normF() < tol){
				break;
			}
		}
	}


	/**
	 * State-Actions-Reward-State-Action is an algorithm for learning a MDP policy, 
	 * which is an on-policy temporal difference control algorithm
	 * In general TD(lambda), V(s_k) = \E[\sum_{t=k}^\infty lambda^{t-k} d_t] + v(s_k)
	 * Reference: NDP Sec 5.3 and 5.4
	 * Reference: Intro. RL Sec.6.4 Sarsa: On-Policy TD Control
	 * @param lambda: TD(lambda)
	 * @param maxIter: max number of iteration until it converges
	 * @param epislon: epislon greedy policy
	 */
	public void SARSA(double lambda, int maxIter, double epislon){
		Matrix Q = new Matrix(numStates, numActions, 0);
		int iter = 0;
		int [][] numUpdates = new int[numStates][numActions];//the number of times state action pair has been visited
		while(iter < maxIter * numActions * numStates){
			iter++;
			int [][] lastUpdates = new int[numStates][numActions];
			Matrix PPi = TransitionGivenPolicy();
			MarkovChain mc = new MarkovChain(transition.get(0).getArray(),  numStates);
 			int s = getRandomState();
			int a = eGreedyAction(s,epislon);
			int t = 0;
			//Vector<Integer> path = mc.trajectory(s0, isTerminate);			
			double [][] e = new double[numStates][numActions];  // eligible trace
			//for(int t = 1; t < path.size(); t++){
			do{//Repeat for each episode	
				t++;
				int s2 = mc.nextState(s,numStates, transition.get(a).getArray());
				int a2 =eGreedyAction(s2,epislon);
				double r = reward.get(s, a);
				double delta = r  + gamma * Q.get(s2,a2) - Q.get(s, a);
				numUpdates[s][a]++;

				e[s][a] *= Math.pow(gamma * lambda, t - lastUpdates[s][a]);
				e[s][a] += 1.0;		
				Q.set(s, a,  Q.get(s, a) +  learningRate(numUpdates[s][a],alphaRate) * delta * e[s][a]);
				lastUpdates[s][a] = t;
			 
				s = s2;
				a = a2;
			}while(isTerminate[s] == false);//end of for t

			if(policyImprovement(Q) == false)
				break;			
		}//end of while
	}

	/**
	 * Compute the model-free (transition probabilities unknown) value function
	 * using q-learning, which is a off-policy temporal difference control algorithm
	 * Reference: NDP, Sec 5.6 Q-Learning
	 * Reference: Intro. RL, Sec 6.5 Q-Learning: Off-Policy TD Control
	 * @param lambda: TD(lambda)
	 * @param maxIter: max number of iteration until it converges
	 * @param epislon: epislon greedy policy
	 */
	public void QLearning(double lambda, int maxIter, double epislon){//Model Free Off Policy learning
		Matrix Q = new Matrix(numStates, numActions, 0);
		int [][] numUpdates = new int [numStates][numActions];
		int iter = 0;
		MarkovChain mc = new MarkovChain(transition.get(0).getArray(), numStates);
		while( iter < maxIter * numActions * numStates){//Repeat for each episode	
			iter++;
			int [][] lastUpdates = new int[numStates][numActions];
			double [][] e = new double[numStates][numActions];
			int s = getRandomState();
			int a = eGreedyAction(s,epislon);
			int t = 0;
			do{//For each step of episode
				t++;
				int s2 = mc.nextState(s,numStates, transition.get(a).getArray());
				int a2 = eGreedyAction(s2, epislon);
				numUpdates[s][a]++;
				double r = reward.get(s, a);
				
				double delta = r  + gamma * getValue(s2) - Q.get(s, a);
				e[s][a] *= Math.pow(gamma * lambda, t - lastUpdates[s][a]);
				e[s][a] += 1.0;		
				Q.set(s, a,  Q.get(s, a) + learningRate(numUpdates[s][a],alphaRate) * delta * e[s][a] );
				
				if(a2 != getPolicy(s2)){
					e = new double[numStates][numActions];
					for(int [] row : lastUpdates)
						Arrays.fill(row,t);
				} 
				lastUpdates[s][a] = t;

				s = s2;
				a = a2;
			}while(isTerminate[s] == false);
			
			if(policyImprovement(Q) == false)
				break;
		}
	}
	
	 
 
}
