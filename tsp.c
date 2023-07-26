#include "tsp.h"

//returns the current amount of seconds passed since the start of the program
double second()
{
	return ((double)clock() / (double)CLK_TCK);
}//second

//Returns a random number bewteen 0.0 and 1.0
double random01() { return ((double)rand() / RAND_MAX); }

//Returns the cost between two nodes by using Euclidian distance formula
double dist(int i, int j, instance* inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j];
	if (!inst->integer_costs) return sqrt(dx * dx + dy * dy);
	int dis = sqrt(dx * dx + dy * dy) + 0.499999999; 					// approximation to the nearest integer 
	return dis + 0.0;
}//dist

// sorts an array of integers in increasing order 
void bubblesort(int* v, int n) {

	for (int i = 0; i < n - 1; i++) {
		for (int k = 0; k < n - 1 - i; k++) {
			if (v[k] > v[k + 1]) {
				int temp = v[k];
				v[k] = v[k + 1];
				v[k + 1] = temp;
			}//if
		}//for
	}//for

}//bubblesort

//sorts an array of solution structures in increasing order of cost
void bubblesort_solution(solution* v, int n) {

	for (int i = 0; i < n - 1; i++) {
		for (int k = 0; k < n - 1 - i; k++) {
			if (v[k].cost > v[k + 1].cost) {
				double temp = v[k].cost;
				v[k].cost = v[k + 1].cost;
				v[k + 1].cost = temp;
			}
		}
	}

}

//Check if our solution contains a single tour
int check_tour(int* succ, int nnodes) {
	int* uncov = (int*)malloc(nnodes * sizeof(int));

	//inizialization uncov
	for (int i = 0; i < nnodes; i++) {
		uncov[i] = 0;
	}//for
	int i = 0;
	while (i < nnodes) {
		if (uncov[succ[i]] == 0) { //first time I visit such node
			uncov[succ[i]] = 1;
			i++;
		}
		else {  //node visit more than once -> not a tour
			free(uncov);
			return 0;
		}//if-else
	}//while
	free(uncov);
	return 1;	// succ contains a tour
}//check_tour

// Calls GNUplot commands to print the solution saved as the best_sol in the instance, or just the points
void plot_solution(instance* inst, int plot_edges)
{
	FILE* fin = fopen("data.dat", "w");
	char coordinates[180];
	for (int i = 0; i < inst->nnodes - 1; i++)
	{
		if (plot_edges) {
			int k = inst->best_sol[i];
			int h = inst->best_sol[i + 1];
			sprintf(coordinates, "%lf %lf\n%lf %lf\n\n", inst->xcoord[k], inst->ycoord[k],
				inst->xcoord[h], inst->ycoord[h]);
		}
		else {
			sprintf(coordinates, "%f %f\n\n", inst->xcoord[i], inst->ycoord[i]);
		}//if-else
		fputs(coordinates, fin);
	}//for
	if (plot_edges) {
		//plot the edge that connects the first node with the last one
		sprintf(coordinates, "%lf %lf\n%lf %lf\n\n", inst->xcoord[inst->best_sol[inst->nnodes - 1]],
			inst->ycoord[inst->best_sol[inst->nnodes - 1]],
			inst->xcoord[inst->best_sol[0]], inst->ycoord[inst->best_sol[0]]);
		fputs(coordinates, fin);
	}//if
	fclose(fin);
	system("gnuplot commands.txt");
}//plot_solution

//precomputes all the distances between nodes and allocates memory and default values for an instance and solution
void preprocessing(instance* inst, solution* sol)
{
	//preprocessiong instance
	inst->cst = (double*)malloc(inst->nnodes * inst->nnodes * sizeof(double));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = 0; j < inst->nnodes; j++)
		{
			inst->cst[i * inst->nnodes + j] = dist(i, j, inst); // precomputes all the distances between nodes
		}//for
	}//for
	inst->best_sol = (int*)malloc(inst->nnodes * sizeof(int));
	inst->best_value = INF;

	// preprocessing solution 
	sol->n = inst->nnodes;
	sol->succ = (int*)malloc(inst->nnodes * sizeof(int)); // best array with successors
	sol->cost = INF;

	// start the time
	inst->tstart = second();
}//preprocessing


// update the incumbent, best solution found so far
void update_incumbent(instance* inst, double cost, int* succ) {
	//update cost
	inst->best_value = cost;

	//update the best solution found so far
	int n = 0;
	inst->best_sol[0] = 0;
	for (int i = 1; i < inst->nnodes; i++)
	{
		inst->best_sol[i] = succ[n];
		n = succ[n];
	}//for
}//update_incumbent

// update the incumbent and plot the solution
void postprocessing(instance* inst, solution* sol) {
	if ((VERBOSE > 60) && check_tour(sol->succ, inst->nnodes)) printf("Solution is a tour\n");
	//update incumbent
	update_incumbent(inst, sol->cost, sol->succ);
	if (VERBOSE >= 10) printf("$STAT %10.6lf\n", inst->best_value);

	//plot the solution
	plot_solution(inst, 1);
}//postprocessing
//converts successor array into optimal node sequence array
void succtopath(int* succ, int* path, int n)
{
	int k = 0;
	path[0] = 0;
	for (int i = 1; i < n; i++)
	{
		path[i] = succ[k];
		k = succ[k];
	}//for
}//succtopath

//converts optimal sequence array into successor array
void pathtosucc(int* path, int* succ, int n)
{
	for (int i = 0; i < n - 1; i++)
	{
		succ[path[i]] = path[i + 1];
	}
	succ[path[n - 1]] = path[0];
}

//returns cost of path from node i to node j
double cost(int i, int j, instance* inst)
{
	return inst->cst[i * inst->nnodes + j];
}//cost

//given the successor vector, compute the cost of the associated solution
double compute_solution_cost(int* succ, instance* inst)
{
	double c = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		c += cost(i, succ[i], inst);
	}//for
	return c;
}//compute_solution_cost

//Computes the GRASP nearest neighbor heuristic on all starting nodes, given a second choice probability
double multistart_grasp_nearest_neighbor(instance* inst, int* best_succ, const double p, int do2opt)
{
	double best_cost = INF;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int)); // local array with successors

	for (int k = 0; (k < inst->nnodes) && (second() - inst->tstart < inst->time_limit); k++) { //for every starting node k

		//compute solution for current node
		double current_cost = grasp_nearest_neighbor(inst, succ, k, p, 0);

		//update the incumbent if improved
		if (current_cost < best_cost) {
			//save best solution so far
			best_cost = current_cost;
			for (int i = 0; i < inst->nnodes; i++) best_succ[i] = succ[i];
		}//if

	}//for
	free(succ);

	//do 2 optimality moves on solution if specified
	if ((VERBOSE > 50) && (do2opt)) printf("****Starting 2opt moves on the best solution found****\n");
	while ((second() - inst->tstart < inst->time_limit) && (do2opt) && (two_optimality_move(inst, best_succ)));
	best_cost = compute_solution_cost(best_succ, inst);

	if (VERBOSE > 40) printf("Best solution found: %10.6lf, time spent = %6.2lf seconds\n", best_cost, second() - inst->tstart);
	return best_cost;
}//multistart_grasp_nearest_neighbor

//Computes the GRASP nearest neighbor heuristic, given a starting node and a second choice probability
double grasp_nearest_neighbor(instance* inst, int* succ, const int startnode, const double p, int do2opt) {

	if (VERBOSE > 50) printf("Starting GRASP Nearest Neighbor heuristic with starting node %d and probability %1.3lf\n", startnode + 1, p);

	int n_uncov = inst->nnodes; //number of uncovered nodes
	int* uncov = (int*)malloc(n_uncov * sizeof(int)); // array with uncovere nodes
	for (int i = 0; i < n_uncov; i++) { //initialize uncovered nodes
		uncov[i] = i;
	}//for

	int current_node = startnode; //the analyzed node

	uncov[startnode] = uncov[--n_uncov]; // starting node is marked as covered

	for (int i = 1; i < inst->nnodes; i++) { //loop through all nodes in order to find the successor of each one

		double min_dist = INF; // the minimum distance is initialized to infinity
		int index_best = -1; //index of the best successor node
		double min_dist2 = INF; // the second smallest distance is initialized to infinity
		int index_best2 = -1;//index of the second best successor node

		for (int j = 0; j < n_uncov; j++) { //loop through all candidates to become the successor of node i
			if (cost(current_node, uncov[j], inst) < min_dist && cost(current_node, uncov[j], inst) < min_dist2) { // update both  distances
				min_dist2 = min_dist;
				index_best2 = index_best;
				min_dist = cost(current_node, uncov[j], inst);
				index_best = j;
			}
			else if (cost(current_node, uncov[j], inst) > min_dist && cost(current_node, uncov[j], inst) < min_dist2) { //only update second best distances
				min_dist2 = cost(current_node, uncov[j], inst);
				index_best2 = j;
			}//if-else
		}//for

		if ((random01() < p) && (index_best2 != -1)) { // update with the second best choice (if available)
			succ[current_node] = uncov[index_best2];
			uncov[index_best2] = uncov[--n_uncov];
		}
		else { // update with the best choice
			succ[current_node] = uncov[index_best];
			uncov[index_best] = uncov[--n_uncov];
		}//if-else
		current_node = succ[current_node]; // update to next node

	}//for
	succ[current_node] = startnode; // closing the cycle
	free(uncov);

	//do 2 optimality moves on solution if specified
	if ((VERBOSE > 50) && (do2opt)) printf("****Starting 2opt moves on found solution****\n");
	while ((second() - inst->tstart < inst->time_limit) && (do2opt) && (two_optimality_move(inst, succ)));
	double cost = compute_solution_cost(succ, inst);

	if (VERBOSE > 50) printf("Solution cost for starting node %d: %10.6lf\n", startnode + 1, cost);
	return cost;
}//grasp_nearest_neighbor

//computes the GRASP Extra Mileage starting from the farthest node couple in the graph
double farthest_grasp_extra_mileage(instance* inst, int* best_succ, const double p, int do2opt)
{
	int* uncov = (int*)malloc(inst->nnodes * sizeof(int)); // array with uncovered nodes
	int n_uncov = inst->nnodes; //number of uncovered nodes
	double maxdist = 0.0; //initialize highest distance found so far
	int n1 = -1; int n2 = -1; //initialize the farthest node couple

	//inizialization of array succ, succ[i] != -1 means that the node i is covered
	for (int i = 0; i < inst->nnodes; i++) {
		best_succ[i] = -1;
	}//for

	//inizialization of uncovered
	for (int i = 0; i < n_uncov; i++) { //initialize uncovered nodes
		uncov[i] = i;
	}//for

	//search for fathest node couple
	for (int i = 0; i < inst->nnodes - 1; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			if (cost(i, j, inst) > maxdist)
			{
				n1 = i;
				n2 = j;
				maxdist = cost(i, j, inst);
			}//if
		}//for
	}//for
	//initialize the starting subtour
	best_succ[n1] = n2;
	best_succ[n2] = n1;
	uncov[n1] = uncov[--n_uncov];
	uncov[n2] = uncov[--n_uncov];

	if (VERBOSE > 50) printf("Starting GRASP Extra Mileage heuristic with starting subtour %d<->%d and probability %1.6lf\n", n1 + 1, n2 + 1, p);
	double total_cost = grasp_extra_mileage(inst, best_succ, p, uncov, n_uncov, 0);
	if (VERBOSE > 50) printf("Solution cost for this subtour: %10.6lf\n", total_cost);

	//do 2 optimality moves on solution if specified
	if ((VERBOSE > 50) && (do2opt)) printf("****Starting 2opt moves on the best solution found****\n");
	while ((second() - inst->tstart < inst->time_limit) && (do2opt) && (two_optimality_move(inst, best_succ)));
	total_cost = compute_solution_cost(best_succ, inst);

	if (VERBOSE > 40) printf("Best solution found: %10.6lf, time spent = %6.2lf seconds\n", total_cost, second() - inst->tstart);
	return total_cost;
}//farthest_grasp_extra_mileage

//computes the GRASP Extra Mileage solution for every possible starting node and the farthest one from it
double multistart_grasp_extra_mileage(instance* inst, int* best_succ, const double p, int do2opt)
{
	double best_cost = INF;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int)); // array with successors
	int* uncov = (int*)malloc(inst->nnodes * sizeof(int)); // array with uncovered nodes

	for (int i = 0; (i < inst->nnodes) && (second() - inst->tstart < inst->time_limit - 5); i++) //leave 5 seconds for refinement at the end
	{
		int n_uncov = inst->nnodes; //number of uncovered nodes
		double maxdist = 0.0;
		int h = -1;  //here the farthest node from i

		//inizialization of array succ, succ[i] != -1 means that the node i is covered
		for (int k = 0; k < inst->nnodes; k++) succ[k] = -1;

		//inizialization of uncovered
		for (int k = 0; k < n_uncov; k++) uncov[k] = k;

		//look for farthest node from the current node i in exam
		for (int j = 0; j < inst->nnodes; j++)
		{
			if (cost(i, j, inst) > maxdist)
			{
				h = j;
				maxdist = cost(i, j, inst);
			}//if
		}//for
		//initialize the starting subtour
		succ[i] = h;
		succ[h] = i;
		uncov[i] = uncov[--n_uncov];
		uncov[h] = uncov[--n_uncov];

		//compute solution for current starting subtour
		if (VERBOSE >= 50) printf("Starting GRASP Extra Mileage heuristic with starting subtour %d<->%d and probability %1.6lf\n", i + 1, h + 1, p);
		double current_cost = grasp_extra_mileage(inst, succ, p, uncov, n_uncov, 0);
		if (VERBOSE >= 50) printf("Solution cost for this subtour: %10.6lf\n", current_cost);

		//update incumbent if improved
		if (current_cost < best_cost)
		{
			//save best solution so far
			best_cost = current_cost;
			for (int k = 0; k < inst->nnodes; k++) best_succ[k] = succ[k];
		}//if

	}//for
	free(succ);
	free(uncov);

	//do 2 optimality moves on solution if specified
	if ((VERBOSE > 50) && (do2opt)) printf("****Starting 2opt moves on the best solution found****\n");
	while ((second() - inst->tstart < inst->time_limit) && (do2opt) && (two_optimality_move(inst, best_succ)));
	best_cost = compute_solution_cost(best_succ, inst);

	if (VERBOSE > 40) printf("Best solution found: %10.6lf, time spent = %6.2lf seconds\n", best_cost, second() - inst->tstart);
	return best_cost;
}//multistart_grasp_extra_mileage

//computes the GRASP extra mileage heuristic, given a partial tour and the array of uncovered nodes, with second choice probability
double grasp_extra_mileage(instance* inst, int* succ, double p, int* uncov, int n_uncov, int do2opt)
{
	while (n_uncov > 0) { // loop for each node left to cover

		double delta = INF; // cost change if we add the best node h to the solution 
		double delta2 = INF; // cost change if we add the second best node h to the solution 
		int h_best = -1, i_best = -1, h_best2 = -1, i_best2 = -1; //inizialization of the indices to -1
		for (int h = 0; h < n_uncov; h++) { // for each uncovered node
			for (int i = 0; i < inst->nnodes; i++) { // for each node in the partial solution
				if (succ[i] == -1) {
					continue;
				}
				else {
					if (cost(i, uncov[h], inst) + cost(uncov[h], succ[i], inst) - cost(i, succ[i], inst) < delta && cost(i, uncov[h], inst) + cost(uncov[h], succ[i], inst) - cost(i, succ[i], inst) < delta2) { //update both 
						//update the second best solution
						delta2 = delta;
						h_best2 = h_best;
						i_best2 = i_best;

						//update the best solution
						delta = cost(i, uncov[h], inst) + cost(uncov[h], succ[i], inst) - cost(i, succ[i], inst);
						h_best = h;
						i_best = i;
					}
					else if (cost(i, uncov[h], inst) + cost(uncov[h], succ[i], inst) - cost(i, succ[i], inst) > delta && cost(i, uncov[h], inst) + cost(uncov[h], succ[i], inst) - cost(i, succ[i], inst) < delta2) {
						//update only the second best solution
						delta2 = cost(i, uncov[h], inst) + cost(uncov[h], succ[i], inst) - cost(i, succ[i], inst);
						h_best2 = h;
						i_best2 = i;
					}//if-else

				}//if-else
			}//for
		}//for
		if ((random01() < p) && (h_best2 != -1)) { // update with the second best choice (if available)
			int temp = succ[i_best2];
			succ[i_best2] = uncov[h_best2];
			succ[uncov[h_best2]] = temp;
			uncov[h_best2] = uncov[--n_uncov];
		}
		else { //update with the best choice
			int temp = succ[i_best];
			succ[i_best] = uncov[h_best];
			succ[uncov[h_best]] = temp;
			uncov[h_best] = uncov[--n_uncov];
		}//if-else

	}//while

	//do 2 optimality moves on solution if specified
	if ((VERBOSE > 50) && (do2opt)) printf("****Starting 2opt moves on the best solution found****\n");
	while ((second() - inst->tstart < inst->time_limit) && (do2opt) && (two_optimality_move(inst, succ)));

	//compute cost of the current solution
	double current_cost = compute_solution_cost(succ, inst);

	return current_cost;
}//grasp_extra_mileage

//inverts the path in the successor array from node start to node end
void invert_succ(int* succ, int start, int end) {
	int current = start;
	int next = succ[start];
	while (current != end)
	{
		int temp = succ[next];
		succ[next] = current;
		current = next;
		next = temp;
	}
}// invert_succ

//refines the solution in the successor array with a 2opt move, if possible
int two_optimality_move(instance* inst, int* succ) {

	double delta = INF;  // the difference between the new cost anche the old one 

	//origins of the best crossing pair we can remove
	int best_a = -1;
	int best_b = -1;

	//double loop for each pair (a,b)
	for (int a = 0; a < inst->nnodes - 1; a++) {
		for (int b = a + 1; b < inst->nnodes; b++) {
			double temp_cost = cost(a, b, inst) + cost(succ[a], succ[b], inst)
				- cost(a, succ[a], inst) - cost(b, succ[b], inst);  //delta cost between the current solution and the old one 
			//update the solution if we have a lower cost
			if ((temp_cost < delta) && ((int)temp_cost < 0)) {  //the second condition is due to the approximation of very small negative numbers
				delta = temp_cost;
				best_a = a;
				best_b = b;
			}//if
		}//for
	}//for

	int a_prime = succ[best_a];
	int b_prime = succ[best_b];

	if (VERBOSE > 70) printf("best a: %d, best b: %d, a prime: %d, b prime: % d\n", best_a, best_b, a_prime, b_prime);
	if ((delta > INF / 2)) {
		return 0; // don't do the optimality move because there are no crossing pairs to remove 
	}

	// step 1: invert the path from a_prime to b
	if (VERBOSE > 70) printf("inverting path from node %d to node %d\n", a_prime, best_b);
	invert_succ(succ, a_prime, best_b);

	// step 2: change the successor of a and a_prime 
	succ[best_a] = best_b;
	succ[a_prime] = b_prime;

	return 1; // two optimality move executed successfully
}//two_optimality_move

//Variant of 2opt for Tabu Search metaheuristic which takes tabu list into account
void two_optimality_move_tabu(instance* inst, int* succ, int tenure, int* tabu_list, int curr_iter) {
	double delta = INF;  // the difference between the new cost anche the old one 

	//origins of the best crossing pair we can remove
	int best_a = -1;
	int best_b = -1;

	//double loop for each pair (a,b)
	for (int a = 0; a < inst->nnodes - 1; a++) {
		if (curr_iter - tabu_list[a] <= tenure) continue;
		for (int b = a + 1; b < inst->nnodes; b++) {
			if (curr_iter - tabu_list[b] <= tenure) continue;
			double temp_cost = cost(a, b, inst) + cost(succ[a], succ[b], inst)
				- cost(a, succ[a], inst) - cost(b, succ[b], inst);  //delta cost between the current solution and the old one 
			//update the solution if we have a lower cost
			if ((temp_cost < delta)) {
				delta = temp_cost;
				best_a = a;
				best_b = b;
			}//if
		}//for
	}//for

	int a_prime = succ[best_a];
	int b_prime = succ[best_b];

	tabu_list[best_a] = curr_iter;
	tabu_list[best_b] = curr_iter;

	if (VERBOSE > 70) printf("best a: %d, best b: %d, a prime: %d, b prime: % d\n", best_a, best_b, a_prime, b_prime);

	// step 1: invert the path from a_prime to b
	if (VERBOSE > 70) printf("inverting path from node %d to node %d\n", a_prime, best_b);
	invert_succ(succ, a_prime, best_b);

	// step 2: change the successor of a and a_prime 
	succ[best_a] = best_b;
	succ[a_prime] = b_prime;


}//two_optimality_move_tabu

//applies Tabu Search metaheuristic with alternating tenure for intensification and diversification
double tabu_search(instance* inst, int* best_succ, int tenure_shift_frequency) {

	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	int min_tenure = 10;
	int max_tenure;

	//inizialization of tenure 
	if (100 < inst->nnodes / 10)max_tenure = 100;
	else max_tenure = inst->nnodes / 10;
	int* tabu_list = (int*)malloc(inst->nnodes * sizeof(int));

	//inizialization of tabu list
	for (int i = 0; i < inst->nnodes; i++)
	{
		tabu_list[i] = -INF;
	}
	int curr_iter = 1;

	//obtain an initial solution 
	double best_cost = multistart_grasp_nearest_neighbor(inst, succ, 0.1, 0);

	if (VERBOSE > 40) printf("Starting tabu search procedure\n");
	// alternation of diversification and intensification phase starting with intensificatio phase
	int tenure = max_tenure;
	while (second() - inst->tstart < inst->time_limit) {
		if (curr_iter % tenure_shift_frequency == 0) { // we change the phase every tot iterations
			if (tenure == max_tenure) tenure = min_tenure; // change into intensification phase
			else tenure = max_tenure; // change into diversification phase
		}//if

		//apply tabu variant of 2opt move
		two_optimality_move_tabu(inst, succ, tenure, tabu_list, curr_iter);

		//update best solution if improved
		double current_cost = compute_solution_cost(succ, inst);
		if (current_cost < best_cost) {
			best_cost = current_cost;
			for (int i = 0; i < inst->nnodes; i++) best_succ[i] = succ[i];
		}//if

		if (VERBOSE > 70) printf("Solution cost of iteration %d: %10.6lf", curr_iter, current_cost);
		curr_iter++;
	}//while
	free(succ);
	if (VERBOSE > 40) printf("tabu search executed in %6.3lf seconds with %d iterations\n ", second() - inst->tstart, curr_iter);
	return compute_solution_cost(best_succ, inst);

}//tabu_search

//executes a five opt kick on the solution passed as parameter
void five_opt_kick(int* succ, int nnodes) {
	if (VERBOSE > 60) printf("Executing 5-opt kick\n");
	int node[5];  // the nodes chosen at random in order to eliminate 5 edges
	int* path = (int*)malloc(nnodes * sizeof(int));  //the path from start to end

	succtopath(succ, path, nnodes);
	int i = 0;

	// chose at random the index of five nodes
	while (i < 5) {
		int next = random01() * (nnodes - 1);
		bool flag = 1;
		for (int j = 0; j < i; j++) { //check if a node pick at random is equal to a node already chosen
			if (next == node[j]) {
				flag = 0;
			}//if
		}//for
		if (flag) {
			if (VERBOSE > 60) printf("found node %d\n", path[next]);
			node[i] = next;
			i++;
		}//if
	}//while

	//order the index -> guarantee to obtain still a valid solution
	bubblesort(node, 5);

	//substitute the index with the node
	for (int j = 0; j < 5; j++) node[j] = path[node[j]];

	//reconnect with 5 new edges (no need to invert tour due to the sequence of chosen nodes)
	int temp = succ[node[0]];
	for (int j = 0; j < 4; j++) succ[node[j]] = succ[node[j + 1]];
	succ[node[4]] = temp;

	if ((VERBOSE > 60) && check_tour(succ, nnodes)) printf("Tour verified\n");
	free(path);
}//five_opt_kick

//applies Variable Neighborhood Search metaheuristic with 5opt-kick for diversification
double vns(instance* inst, int* best_succ) {

	//starting solution
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	double current_cost = multistart_grasp_nearest_neighbor(inst, succ, 0.1, 0);
	int curr_iter = 1;

	if (VERBOSE > 40) printf("Starting variable neighborhood search procedure\n");
	//Apply VNS to the solution
	while (second() - inst->tstart < inst->time_limit) {

		//phase one is intensification-> 2-opt
		while (two_optimality_move(inst, succ));
		//update incumbent
		if (compute_solution_cost(succ, inst) < current_cost) {
			current_cost = compute_solution_cost(succ, inst);
			for (int i = 0; i < inst->nnodes; i++) best_succ[i] = succ[i];
		}//if
		//phase two is diversification-> 5-opt kick
		if (VERBOSE > 70) printf("Starting 5-opt kick procedure\n");
		five_opt_kick(succ, inst->nnodes);
		if (VERBOSE > 70) printf("Phase 2 over\n");
		curr_iter++;
	}//while
	if (VERBOSE > 40) printf("VNS executed in %d iterations in %lf seconds\n", curr_iter, second() - inst->tstart);
	free(succ);
	return current_cost;

}//vns

double simulated_annealing(instance* inst, int* best_succ, double Tmax) {
	//set variables
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));
	double p, best_cost = INF; //p is the probability of choosing the 2opt move even if delta cost is > 0
	double scaler = 25.0; //average edge cost of a TSP heuristic solution with 1000 nodes
	double T = Tmax; //temperature starts off maximum and decreases as we approache the time limit

	//start with a good heuristic solution
	double current_cost = multistart_grasp_nearest_neighbor(inst, succ, 0.1, 0);

	//apply the simulated annealing modified 2opt routine while the time limit holds
	while (second() - inst->tstart < inst->time_limit) {
		if (VERBOSE > 50) printf("=================Starting a Simulated Annealing iteration=====================\n");
		//execute 2opt between 2 random nodes.
		T = Tmax * (inst->time_limit - (second() - inst->tstart)); //update temperature
		int a = random01() * (inst->nnodes - 1);
		int b = random01() * (inst->nnodes - 1);
		while (a == b) b = random01() * (inst->nnodes - 1); //select two random nodes
		double delta_cost = cost(a, b, inst) + cost(succ[a], succ[b], inst)
			- cost(a, succ[a], inst) - cost(b, succ[b], inst);
		int a_prime = succ[a];
		int b_prime = succ[b];
		if (delta_cost < 0) { //normal 2opt move
			invert_succ(succ, a_prime, b);
			succ[a] = b;
			succ[a_prime] = b_prime;
		}//if
		else {
			p = exp(-delta_cost / (scaler * T));
			if (random01() < p) { //diversification 2opt move
				invert_succ(succ, a_prime, b);
				succ[a] = b;
				succ[a_prime] = b_prime;
			}//if
		}//else
		double current_cost = compute_solution_cost(succ, inst);
		if (current_cost < best_cost) {
			best_cost = current_cost;
			for (int i = 0; i < inst->nnodes; i++) best_succ[i] = succ[i];
		}//if
	}//while

	free(succ);
	return best_cost;
}//simulated_annealing

//create a new child via crossing of two parents
void create_child(instance* inst, solution* sol, int parent_1, int parent_2, solution* child) {
	// step 1: chose at random the position in which we cut the "chromosome"
	int pos = inst->nnodes * (0.25 + random01() * 0.5);
	int* path1 = (int*)malloc(inst->nnodes * sizeof(int));
	int* path2 = (int*)malloc(inst->nnodes * sizeof(int));
	int* path_child = (int*)malloc(inst->nnodes * sizeof(int));

	//convert representation to node sequence
	succtopath(sol[parent_1].succ, path1, inst->nnodes);
	succtopath(sol[parent_2].succ, path2, inst->nnodes);


	//merge the two parents at the chosen random point
	for (int i = 0; i < pos; i++) path_child[i] = path1[i];
	for (int i = pos; i < inst->nnodes; i++) path_child[i] = path2[i];

	if (VERBOSE > 100) {
		for (int i = 0; i < inst->nnodes; i++) printf("basic path_child[%d] = %d\n", i, path_child[i]);
	}//if

	free(path1);
	free(path2);

	// step 2: repair the solution
	// repair the solution a) shortcut-> delete the twice visited nodes

	//flag all double nodes as -1
	for (int i = pos; i < inst->nnodes; i++) {
		for (int j = 0; j < pos; j++) {
			if (path_child[j] == path_child[i]) {
				path_child[i] = -1;
			}//if
		}//for
	}//for

	if (VERBOSE > 100) {
		for (int i = 0; i < inst->nnodes; i++) printf("flagged path_child[%d] = %d\n", i, path_child[i]);
	}//if

	//shift all nodes once -1 is found
	int child_nodes = inst->nnodes;

	int i = pos;
	while (i < child_nodes) {
		if (path_child[i] == -1) {
			for (int j = i; j < child_nodes - 1; j++) {
				path_child[j] = path_child[j + 1];
			}//for
			child_nodes--;
		}//if
		else i++;
	}//while

	if (VERBOSE > 100) {
		for (int i = 0; i < inst->nnodes; i++) printf("reduced path_child[%d] = %d\n", i, path_child[i]);
		printf("number of nodes in partial sol: %d\n", child_nodes);
	}//if
	//initialize sucessor vector of the child

	for (int i = 0; i < inst->nnodes; i++) {
		child->succ[i] = -1;
	}//for

	//convert representation to successor vector to use extra_mileage
	pathtosucc(path_child, child->succ, child_nodes);
	free(path_child);

	//initialize uncovered vector to use extra_mileage
	int* uncov = (int*)malloc(inst->nnodes * sizeof(int));
	int n_uncov = 0;
	for (int i = 0; i < inst->nnodes; i++) {
		if (child->succ[i] == -1) uncov[n_uncov++] = i;
	}//for

	// repair the solution b) insert the missing nodes with extra- mileage and c) do full 2-opt
	child->cost = grasp_extra_mileage(inst, child->succ, 0.0, uncov, n_uncov, 1);
	free(uncov);
}//create_child

//Create a new child via mutation
void mutate(instance* inst, solution* sol, int parent, solution* child, int n) {
	//initialize children
	int* path_child = (int*)malloc(inst->nnodes * sizeof(int));

	succtopath(sol[parent].succ, path_child, inst->nnodes);

	for (int i = 0; i < n; i++) {
		int h = random01() * (inst->nnodes - 1);
		int k = random01() * (inst->nnodes - 1);
		int temp = path_child[k];
		path_child[k] = path_child[h];
		path_child[h] = temp;
	}//for
	pathtosucc(path_child, child->succ, inst->nnodes);

	//optimize child
	while (two_optimality_move(inst, child->succ));
	child->cost = compute_solution_cost(child->succ, inst);
	free(path_child);
}//mutate

// genetic algorithm :)
double genetic(instance* inst, int* best_succ, int num_sol) {

	// allocate population and children
	solution* sol = (solution*)malloc(num_sol * sizeof(solution));
	for (int i = 0; i < num_sol; i++) {
		sol[i].succ = (int*)malloc(inst->nnodes * sizeof(int));
	}
	solution* children = (solution*)malloc(num_sol / 2 * sizeof(solution));
	for (int i = 0; i < num_sol / 2; i++) {
		children[i].succ = (int*)malloc(inst->nnodes * sizeof(int));
	}

	// create the initial population of 1000 solutions
	double cost_champion = INF;
	int index_champion = -1;
	for (int i = 0; i < num_sol; i++) {
		sol[i].cost = grasp_nearest_neighbor(inst, sol[i].succ, random01() * (inst->nnodes - 1), 0.2, 0); //create a new solution
		if (sol[i].cost < cost_champion) { //search for the champion 
			cost_champion = sol[i].cost;
			index_champion = i;
		}
	}

	if (VERBOSE > 40) printf("Created initial population\n");

	// at each iteration we update the population with 500 new solutions -> 1500 solution in total
	// so we need to eliminate 500 solution to return to 1000 solution for each iteration 
	int curr_iter = 0;

	while (second() - inst->tstart < inst->time_limit) {
		//CROSSOVER -> 40%  (in our case 400 new solutions)
		// step 1: create the child by choosing at random the parents 
		// and overwrite the solution to bias the deletion high cost (low fitness)
		// chose at random two parents with index i and j 
		int num_child = 0;

		while (num_child < num_sol * 0.4) {
			int i = random01() * (num_sol - 1);
			int j = random01() * (num_sol - 1);
			if (i != j) {
				if (VERBOSE > 70) printf("Creating child %d\n", num_child);
				create_child(inst, sol, i, j, &children[num_child]);
				if ((VERBOSE > 70) && check_tour(children[num_child].succ, inst->nnodes)) printf("Tour verified for child %d\n", num_child);
				num_child++;
			}
			else {
				continue;
			}//if-else
		}//while

		//MUTATION -> 10% (in our case 100 new solutions)

		while (num_child < num_sol * 0.5) {
			//select at random a parent to mutate
			int i = random01() * (num_sol - 1);
			//select at random the number of nodes to mutate
			int n = random01() * (inst->nnodes - 1);
			mutate(inst, sol, i, &children[num_child], n);
			num_child++;
		}//while

		//sort pupolation in order of cost
		bubblesort_solution(sol, num_sol);
		int counter = 0; int bad_sol = 0;
		while (counter < num_child) {
			//skip a bad solution randomly
			if (random01() < 0.05) {
				bad_sol++;
				continue;
			}//if

			//substitute the bad solution with a child
			sol[bad_sol % num_sol].cost = children[counter].cost;
			for (int i = 0; i < inst->nnodes; i++) sol[bad_sol % num_sol].succ[i] = children[counter].succ[i];
			counter++;
			bad_sol++;
		}//for

		//update the champion if necessary
		for (int i = 0; i < num_sol; i++) {
			if (sol[i].cost < cost_champion) { //search for the champion 
				cost_champion = sol[i].cost;
				index_champion = i;
			}//if
		}//for

		curr_iter++; //update the number of iterations
		if (VERBOSE > 40) printf("Advanced to epoch %d, current champion has cost: %lf\n", curr_iter, cost_champion);
	}//while

	//save best solution found
	for (int i = 0; i < inst->nnodes; i++) best_succ[i] = sol[index_champion].succ[i];
	if (VERBOSE > 40) printf("Ended genetic search after %d epochs and %lf seconds\n", curr_iter, second() - inst->tstart);

	//free memory
	for (int i = 0; i < num_sol; i++) {
		free(sol[i].succ);
	}//for
	free(sol);
	for (int i = 0; i < num_sol / 2; i++) {
		free(children[i].succ);
	}//for
	free(children);

	return cost_champion;
}//genetic








/*CPLEX FUNCTIONS**********************************************************************************************************/
double TSPopt_Benders(instance* inst, int* best_succ)
/**************************************************************************************************************************/
{
	// open CPLEX model and create the data structure env, lp
	int error; //contain possible errors
	CPXENVptr env = CPXopenCPLEX(&error); //data structure env: enviroment 
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); //data structure lp: store the model
	if (error) print_error("CPXcreateprob() error");

	//create the model by populating env and lp
	build_model(inst, env, lp);

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 123456);
	CPXsetdblparam(env, CPX_PARAM_TILIM, 3600.0);

	// set lower bound and upper bound
	double LB = 0;
	double UB = INF; // incumbent -> the best integer solution found so far

	while (LB < 0.9999 * UB) {

		//control time limit and update the cplex time limit
		if (second() - inst->tstart > inst->time_limit) break;
		CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit - (second() - inst->tstart));

		if (VERBOSE > 50) printf("Starting CPLEX optimizer\n");

		//set parameters
		CPXsetdblparam(env, lp, CPX_PARAM_CUTUP, UB); //cutoff limit: say to cplex to not provide solution with cost greater than UB
		CPXsetintparam(env, lp, CPX_PARAM_NODELIM, 100); //node limit: CPXmipopt will run for maximum 100 nodes

		//solve the problem
		error = CPXmipopt(env, lp);
		if (error)
		{
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error");
		}//if

		//get the best objective value and update lower bound accordingly
		double obj;
		error = CPXgetbestobjval(env, lp, &obj);
		if (error) {
			printf("CPX error code %d\n", error);
			print_error("CPXgetbestobjval() error");
		}//if
		if (obj > LB) LB = obj;


		// use the optimal solution found by CPLEX
		inst->ncols = CPXgetnumcols(env, lp); // ncols corresponds to the number of variables in the model
		int ncols = inst->ncols;
		double* xstar = (double*)malloc(ncols * sizeof(double));  // here is my solution
		if (CPXgetx(env, lp, xstar, 0, ncols - 1)) continue;	  // populate xstar with the optimal solution found by cplex; last two arguments -> start and end

		if (VERBOSE >= 50) {// print our solution				
			for (int i = 0; i < inst->nnodes; i++)
			{
				for (int j = i + 1; j < inst->nnodes; j++)
				{

					if (xstar[xpos(i, j, inst)] > 0.5) printf("  ... x(%3d,%3d) = 1\n", i + 1, j + 1); // x is in the solution if equal to 1

				}//for
			}//for
		}//if
		//compute the number of components and populate succ and comp 
		int* succ = (int*)malloc(inst->nnodes * sizeof(int)); // array with successors
		int* comp = (int*)malloc(inst->nnodes * sizeof(int)); // array with component numbers for each node
		int ncomp;											  //number of separate subtours
		build_sol(xstar, inst, succ, comp, &ncomp);

		// Bender's loop
		if (ncomp >= 2) {
			//Add SECs for all connected components found
			if (VERBOSE > 40) printf("Found %d separate subtours\n", ncomp);
			for (int k = 1; k <= ncomp; k++) {
				add_SEC(inst, env, lp, comp, succ, k, NULL); //add sec for component comp
			}//for

			//Patch up the current collection of subtours into a feasible TSP solution and update incumbent
			patchingHeuristic(inst, env, lp, succ, comp, ncomp, NULL);
		}//if
		double current_cost = compute_solution_cost(succ, inst);
		if (current_cost < UB) {
			UB = current_cost;
			for (int k = 0; k < inst->nnodes; k++) best_succ[k] = succ[k];
		}//if

		//free memory
		free(xstar);
		free(succ);
		free(comp);
		if (VERBOSE > 50) printf("Current LB = %lf, current UB = %lf\n", LB, UB);
	}//while

	if ((VERBOSE > 10) && (fabs(UB - LB) < 0.0001)) printf("Found OPTIMAL solution of cost %lf\n", UB);
	if ((VERBOSE > 10) && (fabs(UB - LB) > 0.0001)) printf("Found NOT OPTIMAL solution of cost %lf\n", UB);


	// free and close cplex model   
	CPXfreeprob(env, &lp); //free lp data structure
	CPXcloseCPLEX(&env); //free env data structure
	return UB;

}//TSPopt_Benders

double TSPopt_branch_and_cut(instance* inst, int* best_succ)
/**************************************************************************************************************************/
{
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error); // data structure envi
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);  // you need to create your own build_model in order to populate env and lp

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 0);
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit - (second() - inst->tstart));

	inst->ncols = CPXgetnumcols(env, lp);
	int ncols = inst->ncols;

	//inizialization
	double* xstar = (double*)malloc(ncols * sizeof(double));	// here is my solution
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));		// array with successors
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));		// array with component numbers for each node
	int ncomp;													// number of separate subtours

	// install a "lazyconstraint" callback to cut infeasible integer sol.s (found e.g. by heuristics) 
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE; // only want to call the callback when cplex found a candidate solution
	if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) print_error("CPXcallbacksetfunc() error");

	//initialize CPLEX with a "good" heuristic solution
	double best_cost = multistart_grasp_nearest_neighbor(inst, succ, 0.1, 1);

	//Converting a "successor" TSP solution succ[0..nnodes-1] to a CPLEX solution xheu[0..ncols-1]
	double* xheu = (double*)malloc(ncols * sizeof(double));
	for (int i = 0; i < ncols; i++) xheu[i] = 0.0;								//initialize xheu to all zeros
	for (int i = 0; i < inst->nnodes; i++) xheu[xpos(i, succ[i], inst)] = 1.0;	//set solution edges to one
	int* ind = (int*)malloc(ncols * sizeof(int));
	for (int i = 0; i < ncols; i++) ind[i] = i; //indices
	int effortlevel = CPX_MIPSTART_NOCHECK; //cplex don't check if the incumbent solution is feasible
	int beg = 0;
	//provide to cplex an incumbent solution
	if (CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");
	free(ind);
	free(xheu);

	//CPXsetintparam(env, CPX_PARAM_THREADS, 1); 	// just for debugging, only one thread at time
	CPXmipopt(env, lp); 							// cplex will eventually use the callback installed

	if (CPXgetx(env, lp, xstar, 0, ncols - 1)) { // optimal TSP sol. (or just feasible if time limit...), if any
		if (VERBOSE > 40) printf("we don't have any better solution\n");
		for (int k = 0; k < inst->nnodes; k++) best_succ[k] = succ[k];
	}
	else {
		build_sol(xstar, inst, succ, comp, &ncomp);
		best_cost = compute_solution_cost(succ, inst);
		for (int k = 0; k < inst->nnodes; k++) best_succ[k] = succ[k];
	}//if-else

	// free and close cplex model   
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);

	//free memory
	free(xstar);
	free(succ);
	free(comp);

	return best_cost;
}//TSPopt_branch_and_cut

/********************************************************************************************************/
static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle)
/********************************************************************************************************/
{
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");

	instance* inst = (instance*)userhandle; //cast
	int ncols = inst->ncols;
	int nnodes = inst->nnodes;
	double* xstar = (double*)malloc(ncols * sizeof(double));
	double objval = CPX_INFBOUND;
	//get integer candidate solution from cplex
	if (CPXcallbackgetcandidatepoint(context, xstar, 0, ncols - 1, &objval)) print_error("CPXcallbackgetcandidatepoint error"); // is used to obtain the solution xstar

	int* succ = (int*)malloc(inst->nnodes * sizeof(int)); // array with successors
	int* comp = (int*)malloc(inst->nnodes * sizeof(int)); // array with component numbers for each node
	int ncomp;											  // number of separate subtours

	build_sol(xstar, inst, succ, comp, &ncomp);
	// if the solution xstar contains more than one component we need to add SEC's 
	if (ncomp >= 2) {

		if (VERBOSE > 40) printf("Found %d separate subtours\n", ncomp);
		for (int k = 1; k <= ncomp; k++) {
			add_SEC(inst, env, lp, comp, succ, k, context);
		}//for
	}//if

	//Repair the solution and post it to CPLEX if better than incumbent, eventually this solution will be used by cplex
	if (VERBOSE > 40)printf("Repairing the temporary heuristic solution\n");
	patchingHeuristic(inst, env, lp, succ, comp, ncomp, context);
	double objheu = compute_solution_cost(succ, inst);
	double incumbent = CPX_INFBOUND;
	CPXcallbackgetinfodbl(context, CPXCALLBACKINFO_BEST_SOL, &incumbent);
	if (objheu < 0.9999 * incumbent) {
		if (VERBOSE > 40) printf("Posting the heuristic solution to CPLEX\n");
		double* xheu = (double*)malloc(ncols * sizeof(double));
		//Converting a "successor" TSP solution succ[0..nnodes-1] to a CPLEX solution xheu[0..ncols-1]
		for (int i = 0; i < ncols; i++) xheu[i] = 0.0;							//initialize xheu to all zeros
		for (int i = 0; i < nnodes; i++) xheu[xpos(i, succ[i], inst)] = 1.0;	//set solution edges to one
		int* ind = (int*)malloc(ncols * sizeof(int));
		for (int j = 0; j < ncols; j++) ind[j] = j;
		if (CPXcallbackpostheursoln(context, ncols, ind, xheu, objheu, CPXCALLBACKSOLUTION_NOCHECK)) print_error("CPXcallbackpostheursoln() error");

		//free memory
		free(ind);
		free(xheu);
	}//if
	free(xstar);
	return 0;
}//my_callback

void add_SEC(instance* inst, CPXENVptr env, CPXLPptr lp, int* comp, int* succ, int mycomp, CPXCALLBACKCONTEXTptr context) {
	//initialize parameters of the row to add
	int ncols = inst->ncols;
	int* index = (int*)malloc(ncols * sizeof(int)); //index of the column of the variable in the new constraint
	double* coeff = (double*)malloc(ncols * sizeof(double)); //coefficient of the variable in the constraint
	int nnz = 0; //number of nonzero entries in the new constraint
	double rhs = -1.0; //right hand side of new constraint

	for (int i = 0; i < inst->nnodes; i++) {
		if (comp[i] != mycomp) continue;
		rhs = rhs + 1.0; //cardinality of the set
		for (int j = i + 1; j < inst->nnodes; j++) {
			if (comp[j] != mycomp) continue;
			index[nnz] = xpos(i, j, inst);
			coeff[nnz] = 1.0;
			nnz++;
		}//for
	}//for

	//constraint definition
	char sense = 'L';
	int izero = 0;
	char** rname = (char**)malloc(1 * sizeof(char*));		// (char **) required by cplex... row name
	rname[0] = (char*)malloc(100 * sizeof(char));
	sprintf(rname[0], "component(%d)", mycomp);

	//add the SEC for current component k
	if (context == NULL) {
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, coeff, NULL, &rname[0])) print_error("CPXaddrows(): error 1");
	}
	else {
		if (CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, coeff)) print_error("CPXcallbackrejectcandidate() error"); // reject the solution and adds one cut 
	}//if-else
	free(rname[0]);
	free(rname);
	free(index);
	free(coeff);
	if (VERBOSE > 60) printf("Added row for component %d\n", mycomp);
}//add_SEC

//In cplex there is a matrix where each column corresponds to x(i, j), this function outputs the column of a given x(i, j)
/***************************************************************************************************************************/
int xpos(int i, int j, instance* inst)
/***************************************************************************************************************************/
{
	if (i == j) print_error(" i == j in xpos");
	if (i > j) return xpos(j, i, inst); // invert the index -> undirected graph
	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
	return pos;
}//xpos


/***************************************************************************************************************************/
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp)
/**************************************************************************************************************************/
{

	double zero = 0.0;
	char binary = 'B'; //because our variables can only be equal to 0 or 1

	char** cname = (char**)malloc(1 * sizeof(char*));		// (char **) required by cplex...
	cname[0] = (char*)malloc(100 * sizeof(char)); //name of the variables, normally matrix


	// add columns corresponding to binary var.s x(i,j) for i < j   -> remember that graph is undirected

	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			sprintf(cname[0], "x(%d,%d)", i + 1, j + 1);  		// print a value in a string; start from x(1,2) and not x(0,1)
			double obj = dist(i, j, inst); // cost == distance   
			double lb = 0.0; //remember variable is binary 
			double ub = 1.0;
			//create a new column corrisponding to a variable
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(" wrong CPXnewcols on x var.s"); // the third argument corresponds to the number of colums to add
			// debug
			if (CPXgetnumcols(env, lp) - 1 != xpos(i, j, inst)) print_error(" wrong position for x var.s");// in this way we are sure that xpos function is correctly implemented 

		}//for
	}//for


	// add the degree constraints 

	int* index = (int*)malloc(inst->nnodes * sizeof(int)); //index of the variables
	double* value = (double*)malloc(inst->nnodes * sizeof(double)); //value of the coefficient

	for (int h = 0; h < inst->nnodes; h++)  		// add the degree constraint on node h
	{
		double rhs = 2.0; //righthand side is 2 (remember our model) 
		char sense = 'E'; // 'E' for equality constraint 
		sprintf(cname[0], "degree(%d)", h + 1);
		int nnz = 0; //number of non zero
		//I only have to specify the constraints different from zero
		//other coefficients are automatically equal to zero
		for (int i = 0; i < inst->nnodes; i++)
		{
			if (i == h) continue;
			index[nnz] = xpos(i, h, inst); //index of the variable of nnx coeffiecient
			value[nnz] = 1.0; //value of the coefficient
			nnz++;
		}
		int izero = 0;
		//add constraints
		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])) print_error("CPXaddrows(): error 1"); // 0 stays for #cols and 1 stays for #rows
	}//for

	free(value);
	free(index);

	free(cname[0]);
	free(cname);

	if (VERBOSE >= 100) CPXwriteprob(env, lp, "model.lp", NULL);   //write the model in a file, in this case in model.lp

}//build_model

/*********************************************************************************************************************************/
void build_sol(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp) // build succ() and comp() wrt xstar()...
/*********************************************************************************************************************************/
{
	// 1) check that the degree is 2 (one entering edge and one leaving edge) 
	// 2) verify that no fractional components are present

#ifdef DEBUG
	int* degree = (int*)calloc(inst->nnodes, sizeof(int));
	for (int i = 0; i < inst->nnodes; i++)
	{
		for (int j = i + 1; j < inst->nnodes; j++)
		{
			int k = xpos(i, j, inst);
			if (fabs(xstar[k]) > EPS && fabs(xstar[k] - 1.0) > EPS) print_error(" wrong xstar in build_sol()"); //check if there is fractional component
			// the edge is in the solution -> increase degree
			if (xstar[k] > 0.5)
			{
				++degree[i];
				++degree[j];
			}//if
		}//for
	}//for
	// control if the degree is 2
	for (int i = 0; i < inst->nnodes; i++)
	{
		if (degree[i] != 2) print_error("wrong degree in build_sol()");
	}//for
	free(degree);
#endif
	// 3) compute comp, succ and ncomp

	// 3a) inizialization
	* ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++)
	{
		//all nodes unvisited
		succ[i] = -1;
		comp[i] = -1;
	}//for

	// 3b) populate comp,succ and calculate ncomp
	for (int start = 0; start < inst->nnodes; start++)
	{
		if (comp[start] >= 0) continue;  // node "start" was already visited, just skip it

		// a new component is found-> we explore this component and add all nodes to it 
		(*ncomp)++;
		int i = start;
		int done = 0;
		while (!done)  // go and visit the current component
		{
			comp[i] = *ncomp;
			done = 1;
			for (int j = 0; j < inst->nnodes; j++)
			{
				if (i != j && xstar[xpos(i, j, inst)] > 0.5 && comp[j] == -1) // the edge [i,j] is selected in xstar and j was not visited before 
				{
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}//if
			}//for
		}//while
		succ[i] = start;  // last arc to close the cycle

		// go to the next component...
	}//for
}

void patchingHeuristic(instance* inst, CPXENVptr env, CPXLPptr lp, int* succ, int* comp, const int ncomp, CPXCALLBACKCONTEXTptr context) {
	int localcomp = ncomp;
	if (VERBOSE > 50) printf("Patching up %d subtours\n", localcomp);
	while (localcomp >= 2) {

		double delta = INF;
		int best_a = -1;
		int best_b = -1;
		// look for the minimum delta cost couple
		for (int a = 0; a < inst->nnodes; a++) {
			for (int b = 0; b < inst->nnodes; b++) {

				if (comp[a] < comp[b]) {
					double temp_cost = cost(a, succ[b], inst) + cost(succ[a], b, inst)
						- cost(a, succ[a], inst) - cost(b, succ[b], inst);  //delta cost between the current solution and the old one
					if (temp_cost < delta) {
						delta = temp_cost;
						best_a = a;
						best_b = b;
					}//if
				}//if
			}//for
		}//for
		// update our solution
		int a_prime = succ[best_a];
		int b_prime = succ[best_b];
		int COMPA = comp[best_a];
		int h = b_prime;
		// fusion of two components
		while (h != best_b) {
			comp[h] = COMPA;
			h = succ[h];
		}
		comp[best_b] = COMPA;
		succ[best_a] = b_prime;
		succ[best_b] = a_prime;

		localcomp--; // we eliminate one component

		if (localcomp != 1) add_SEC(inst, env, lp, comp, succ, COMPA, context); //add SEC corresponding to new subtour	
	}//while
	while (two_optimality_move(inst, succ)); //Refine the solution with 2opt
	if (VERBOSE > 50 && check_tour(succ, inst->nnodes)) printf("Tour verified\n");
}//patchingHeuristic

//*******************************************************************
//MATHEURISTICS
//*******************************************************************

double hard_fixing(instance* inst, int* best_succ, double pfix, double cplex_tl) {

	double best_cost = INF;
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);  // you need to create your own build_model in order to populate env and lp

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 0);

	CPXsetdblparam(env, CPX_PARAM_TILIM, cplex_tl);

	inst->ncols = CPXgetnumcols(env, lp);
	int ncols = inst->ncols;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));		// array with successors
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));		// array with component numbers for each node
	int ncomp;													// number of separate subtours

	//find a good starting heuristic solution
	double current_cost = multistart_grasp_nearest_neighbor(inst, succ, 0.1, 1);
	double* xheu = (double*)malloc(ncols * sizeof(double));
	for (int i = 0; i < ncols; i++) xheu[i] = 0.0;								//initialize xheu to all zeros
	for (int i = 0; i < inst->nnodes; i++) xheu[xpos(i, succ[i], inst)] = 1.0;	//set solution edges to one
	int* ind = (int*)malloc(ncols * sizeof(int));
	for (int i = 0; i < ncols; i++) ind[i] = i;
	int effortlevel = CPX_MIPSTART_NOCHECK; //cplex don't check the incumbent solution
	int beg = 0;

	while (second() - inst->tstart < inst->time_limit) {

		if (VERBOSE > 40) printf("============================STARTING A HARD FIXING LOOP===============================\n");
		// 1 - post current heuristic solution to cplex
		if (CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");


		if (VERBOSE >= 40) printf("fixing random edges with probability %1.2lf\n", pfix);
		// 2 - fix random edges in the heuristic solution with probability pfix
		for (int i = 0; i < ncols; i++) {
			if (xheu[i] > 0.5 && random01() < pfix) {
				int index = i; // index of the variable in the matrix of cplex
				char lu = 'L';  //change lower bound
				double bd = 1.0; // set to 1 the lowe bound
				if (CPXchgbds(env, lp, 1, &index, &lu, &bd)) print_error("Error while fixing the edges");
			}
		}

		// 3 - call CPLEX callback with SMALL time limit
		// install a "lazyconstraint" callback to cut infeasible integer sol.s (found e.g. by heuristics) 
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE; // ... means lazyconstraints
		if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) print_error("CPXcallbacksetfunc() error");

		//CPXsetintparam(env, CPX_PARAM_THREADS, 1); 	// just for debugging, only one thread at time
		CPXmipopt(env, lp); 							// with the callback installed

		if (CPXgetx(env, lp, xheu, 0, ncols - 1)) { // optimal TSP sol. (or just feasible if time limit...), if any
			if (VERBOSE > 40) printf("we don't have any better solution\n");
		}
		else {
			build_sol(xheu, inst, succ, comp, &ncomp);
			current_cost = compute_solution_cost(succ, inst);
			if (current_cost < 0.9999 * best_cost) {
				best_cost = current_cost;
				for (int k = 0; k < inst->nnodes; k++) best_succ[k] = succ[k];
			}
		}

		if (VERBOSE >= 40) printf("Unfixing all the edges\n");
		// 4 - unfix all variables
		for (int i = 0; i < ncols; i++) {
			int index = i;
			char lu = 'L';
			double bd = 0.0;
			if (CPXchgbds(env, lp, 1, &index, &lu, &bd)) print_error("Error while unfixing the edges");
		}
	}
	free(ind);
	free(xheu);
	free(succ);
	free(comp);
	return best_cost;
}//hard_fixing

double local_branching(instance* inst, int* best_succ, int k, double cplex_tl) {
	double best_cost = INF;
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);  // you need to create your own build_model in order to populate env and lp

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 0);

	CPXsetdblparam(env, CPX_PARAM_TILIM, cplex_tl);				//SHORT time limit for callback

	inst->ncols = CPXgetnumcols(env, lp);
	int ncols = inst->ncols;
	int n = inst->nnodes;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));		// array with successors
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));		// array with component numbers for each node
	int ncomp;														// number of separate subtours

	//find a good starting heuristic solution
	double current_cost = multistart_grasp_nearest_neighbor(inst, succ, 0.1, 1);
	double* xheu = (double*)malloc(ncols * sizeof(double));
	for (int i = 0; i < ncols; i++) xheu[i] = 0.0;								//initialize xheu to all zeros
	for (int i = 0; i < inst->nnodes; i++) xheu[xpos(i, succ[i], inst)] = 1.0;	//set solution edges to one
	int* ind = (int*)malloc(ncols * sizeof(int));
	for (int i = 0; i < ncols; i++) ind[i] = i;
	int effortlevel = CPX_MIPSTART_NOCHECK;
	int beg = 0;

	while (second() - inst->tstart < inst->time_limit) {
		if (VERBOSE > 40) printf("============================STARTING A LOCAL BRANCHING LOOP===============================\n");

		// 1 - post current heuristic solution to cplex
		if (CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");

		//2 - add local branching contraint
		int nnz = 0;
		double rhs = (double)(n - k);
		char sense = 'G';
		int* index = (int*)malloc(ncols * sizeof(int));
		double* coeff = (double*)malloc(ncols * sizeof(double));
		int izero = 0;
		for (int i = 0; i < ncols; i++) {
			if (xheu[i] > 0.5) {
				index[nnz] = i;
				coeff[nnz] = 1.0;
				nnz++;
			}
		}
		char** rname = (char**)malloc(1 * sizeof(char*));		// (char **) required by cplex... row name
		rname[0] = (char*)malloc(100 * sizeof(char));
		sprintf(rname[0], "Local Branching constraint");

		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, coeff, NULL, &rname[0]))
			print_error("CPXaddrows(): error 1");
		if (VERBOSE > 50) printf("&&&&&&&&&&&&&&&&& Added local branching constraint with k = %d\n", k);
		free(index);
		free(coeff);
		free(rname[0]);
		free(rname);

		//3 - apply CPLEX callback for SHORT time limit
		// install a "lazyconstraint" callback to cut infeasible integer sol.s (found e.g. by heuristics) 
		CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE; // ... means lazyconstraints
		if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) print_error("CPXcallbacksetfunc() error");

		//CPXsetintparam(env, CPX_PARAM_THREADS, 1); 	// just for debugging, only one thread at time
		CPXmipopt(env, lp); 							// with the callback installed

		//4 - Get new solution, update if necessary and remove old local branching constraint
		if (CPXgetx(env, lp, xheu, 0, ncols - 1)) { // optimal TSP sol. (or just feasible if time limit...), if any
			if (VERBOSE > 40) printf("we don't have any better solution\n");
		}
		else {
			build_sol(xheu, inst, succ, comp, &ncomp);
			current_cost = compute_solution_cost(succ, inst);
			if (current_cost < 0.9999 * best_cost) {
				best_cost = current_cost;
				for (int k = 0; k < inst->nnodes; k++) best_succ[k] = succ[k];
			}//if
		}//else
		int numrows = CPXgetnumrows(env, lp);
		if (CPXdelrows(env, lp, numrows - 1, numrows - 1)) print_error("CPXdelrows(): error 1");
		if (VERBOSE > 50) printf("&&&&&&&&&&&&&&&&& deleted last local branching constraint\n");
	}//while

	free(ind);
	free(xheu);
	free(succ);
	free(comp);
	return best_cost;
}//local_branching

double rins_polishing(instance* inst, int* best_succ, int frequency, int t_p) {
	double best_cost = INF;
	// open CPLEX model
	int error;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);  // you need to create your own build_model in order to populate env and lp

	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	if (VERBOSE >= 60) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, 0);
	CPXsetdblparam(env, CPX_PARAM_TILIM, inst->time_limit);				//SHORT time limit for callback
	CPXsetintparam(env, CPXPARAM_MIP_Strategy_RINSHeur, frequency); // apply rins with a certain node frequency 
	CPXsetintparam(env, CPXPARAM_MIP_PolishAfter_Time, t_p); // apply rins with a certain node frequency 

	inst->ncols = CPXgetnumcols(env, lp);
	int ncols = inst->ncols;
	int n = inst->nnodes;
	int* succ = (int*)malloc(inst->nnodes * sizeof(int));		// array with successors
	int* comp = (int*)malloc(inst->nnodes * sizeof(int));		// array with component numbers for each node
	int ncomp;														// number of separate subtours

	//find a good starting heuristic solution
	double current_cost = multistart_grasp_nearest_neighbor(inst, succ, 0.1, 1);
	double* xheu = (double*)malloc(ncols * sizeof(double));
	for (int i = 0; i < ncols; i++) xheu[i] = 0.0;								//initialize xheu to all zeros
	for (int i = 0; i < inst->nnodes; i++) xheu[xpos(i, succ[i], inst)] = 1.0;	//set solution edges to one
	int* ind = (int*)malloc(ncols * sizeof(int));
	for (int i = 0; i < ncols; i++) ind[i] = i;
	int effortlevel = CPX_MIPSTART_NOCHECK;
	int beg = 0;

	if (VERBOSE > 40) printf("============================STARTING RINS_POLISHING===============================\n");

	// 1 - post current heuristic solution to cplex
	if (CPXaddmipstarts(env, lp, 1, ncols, &beg, ind, xheu, &effortlevel, NULL)) print_error("CPXaddmipstarts() error");

	//2 - apply CPLEX callback for SHORT time limit
	// install a "lazyconstraint" callback to cut infeasible integer sol.s (found e.g. by heuristics) 
	CPXLONG contextid = CPX_CALLBACKCONTEXT_CANDIDATE; // ... means lazyconstraints
	if (CPXcallbacksetfunc(env, lp, contextid, my_callback, inst)) print_error("CPXcallbacksetfunc() error");

	//CPXsetintparam(env, CPX_PARAM_THREADS, 1); 	// just for debugging, only one thread at time
	CPXmipopt(env, lp); 							// with the callback installed

	//3- Get new solution, update if necessary 
	if (CPXgetx(env, lp, xheu, 0, ncols - 1)) { // optimal TSP sol. (or just feasible if time limit...), if any
		if (VERBOSE > 40) printf("we don't have any better solution\n");
	}
	else {
		build_sol(xheu, inst, succ, comp, &ncomp);
		current_cost = compute_solution_cost(succ, inst);
		if (current_cost < 0.9999 * best_cost) {
			best_cost = current_cost;
			for (int k = 0; k < inst->nnodes; k++) best_succ[k] = succ[k];
		}//if
	}//if-else

	free(ind);
	free(xheu);
	free(succ);
	free(comp);
	return best_cost;
}//rins_polishing