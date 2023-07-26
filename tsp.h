#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <cplex.h> 
#include <stdbool.h>
#include <time.h>
#include <limits.h>
#define INF 1000000000.0

#define VERBOSE				    60		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)

//data structures

typedef struct {
	int* succ; //solution as succ
	int n;  // number of nodes 
	double cost; // cost of the solution
} solution;

typedef struct {

	//input data
	int nnodes;								//number of nodes
	double* xcoord;							// coordinate x
	double* ycoord;							// coordinate y
	double* cst;                            // all the distances between nodes

	// parameters 
	int randomseed;							// random instance
	char input_file[1000];		  			// input file
	int mode;								// solver mode selection variable
	double time_limit;                      // time limit
	double tstart;							// starting time of the solution computation
	int integer_costs; 						// = 1 for integer costs (rounded distances), 0 otherwise
	int p;									// = 1 print the nodes of the instance on screen
	int ncols;								// number of columns in the CPLEX model
	

	//global data
	double best_value;						//cost of the best solution
	int* best_sol;							//array containing the optimal sequence of nodes

} instance;

//function declarations

double second();

double random01();

/*!
* @function	dist
* @abstract	returns the distance between two points in the instance
* @param i	: index of node 1
* @param j	: index of node 2
* @param inst	: pointer to the instance in question
* @result	the double value of the distance
*/
double dist(int i, int j, instance* inst);

/*!
* @function	bubblesort
* @abstract sorts an array of integers in increasing order
* @param v	: pointer to the array in question
* @param n	: length of the array
*/
void bubblesort(int* v, int n);

/*!
* @function	bubblesort_solution
* @abstract sorts an array of solution structures in increasing order of cost
* @param v	: pointer to the array in question
* @param n	: length of the array
*/
void bubblesort_solution(solution* v, int n);

/*!
* @function	check_tour
* @abstract checks if a successor array consitutes a tour
* @param succ	: pointer to the successor array in question
* @param nnodes	: length of the array
* @result 1 if the succ array constitutes a tour, 0 otherwise
*/
int check_tour(int* succ, int nnodes);

/*!
* @function	plot_solution
* @abstract prints the arcs and nodes of the solution saved in the instance to a GNUplot window
* @param inst		: pointer to the instance structure in question
* @param plot_edges	: 1 for plotting the solution in inst->best sol, 0 if you only want to plot points
*/
void plot_solution(instance* inst, int plot_edges);

/*!
* @function	preprocessing
* @abstract precomputes the distances between all possible point couples in the graph
* @param inst	: pointer to the instance structure in question
* @param sol	: pointer to the solution structure in question
*/
void preprocessing(instance* inst, solution* sol);

/*!
* @function	update_incumbent
* @abstract updates the incumbent solution inside the instance structure with the values passed as parameters
* @param inst	: pointer to the instance structure in question
* @param cost	: double value of the solution cost
* @param succ	: successor array of the solution
*/
void update_incumbent(instance* inst, double cost, int* succ);

/*!
* @function	postprocessing
* @abstract updates the incumbent and plot the solution
* @param inst	: pointer to the instance structure in question
* @param sol	: pointer to the solution structure in question
*/
void postprocessing(instance* inst, solution* sol);

/*!
* @function succtopath
* @abstract takes the successor representation of a solution and converts it into a sequential path representation
*			starting from node 0
* @param succ	: pointer to the successor array in question
* @param path	: pointer to the path array to save the new representation into
* @param n	: length of both arrays
*/
void succtopath(int* succ, int* path, int n);

/*!
* @function pathtosucc
* @abstract takes the sequential path representation of a solution and converts it into a succ representation
* @param path	: pointer to the path array in question
* @param succ	: pointer to the succ array to save the new representation into
* @param n	: length of both arrays
*/
void pathtosucc(int* path, int* succ, int n);

/*!
* @function cost
* @abstract fetches the cost of the arc from node i to j from the instance structure
* @param i	: index of node 1
* @param j	: index of node 2
* @param inst	: pointer to the instance structure in question
* @result the double cost of arc (i,j)
*/
double cost(int i, int j, instance* inst);

/*!
* @function compute_solution_cost
* @abstract computes the cost of the solution saved in the successor array
* @param succ	: pointer to the successor array with the saved solution
* @param inst	: pointer to the instance structure to access the edge costs
* @resut the double cost of the solution in the successor array
*/
double compute_solution_cost(int* succ, instance* inst);

/*!
* @function	multistart_grasp_nearest_neighbor
* @abstract computes the GRASP Nearest Neighbor heuristic on every possible starting node
* @param inst		: pointer to the instance structure in question
* @param best_succ	: successor array where the best solution will be saved
* @param p			: probability to pick the second best move at every step of Nearest Neighbor
* @param do2opt		: 1 if you want 2opt to be applied to the best solution, 0 otherwise
* @result the cost of the best solution found by the algorithm
*/
double multistart_grasp_nearest_neighbor(instance* inst, int* best_succ, const double p, int do2opt);

/*!
* @function	grasp_nearest_neighbor
* @abstract computes the GRASP Nearest Neighbor heuristic on the given starting node
* @param inst		: pointer to the instance structure in question
* @param succ		: successor array where the solution will be saved
* @param startnode	: the index of the starting node for the heuristic (from 0 to nnodes-1)
* @param p			: probability to pick the second best move at every step of Nearest Neighbor
* @param do2opt		: 1 if you want 2opt to be applied to the solution, 0 otherwise
* @result the cost of the solution found by the algorithm
*/
double grasp_nearest_neighbor(instance* inst, int* succ, const int startnode, const double p, int do2opt);

/*!
* @function	farthest_grasp_extra_mileage
* @abstract computes the GRASP Extra Mileage heuristic starting from the subtour composed by 
			the 2 farthest nodes in the graph
* @param inst		: pointer to the instance structure in question
* @param best_succ	: successor array where the solution will be saved
* @param p			: probability to pick the second best move at every step of Extra Mileage
* @param do2opt		: 1 if you want 2opt to be applied to the solution, 0 otherwise
* @result the cost of the solution found by the algorithm
*/
double farthest_grasp_extra_mileage(instance* inst, int* best_succ, const double p, int do2opt);

/*!
* @function	multistart_grasp_extra_mileage
* @abstract computes the GRASP Extra Mileage heuristic from every possible starting subtour composed by 
			a starting node and the farthest node from it
* @param inst		: pointer to the instance structure in question
* @param best_succ	: successor array where the best solution will be saved
* @param p			: probability to pick the second best move at every step of Extra Mileage
* @param do2opt		: 1 if you want 2opt to be applied to the best solution, 0 otherwise
* @result the cost of the best solution found by the algorithm
*/
double multistart_grasp_extra_mileage(instance* inst, int* best_succ, const double p, int do2opt);

/*!
* @function	grasp_extra_mileage
* @abstract computes the GRASP Extra Mileage heuristic starting from the subtour composed by
			the 2 farthest nodes in the graph
* @param inst		: pointer to the instance stricture in question
* @param succ		: successor array where the solution will be saved
* @param p			: probability to pick the second best move at every step of Extra Mileage
* @param uncov		: pointer to the array of uncovered nodes
* @param n_uncov	: length of the array of uncovered nodes
* @param do2opt		: 1 if you want 2opt to be applied to the solution, 0 otherwise
* @result the cost of the solution found by the algorithm
*/
double grasp_extra_mileage(instance* inst, int* succ, double p, int* uncov, int n_uncov, int do2opt);

/*!
* @function invert_succ
* @abstract inverts the path from a starting node to an end node in the successor representation
* @param succ	: pointer to the succ array where the subath to invert is
* @param start	: index of the starting node from which to invert the path
* @param end	: index of the ending node up to which it will invert the path
*/
void invert_succ(int* succ, int start, int end);

/*!
* @function two_optimality_move
* @abstract performs one two optimality move on the provided solution
* @param inst	: pointer to the instance structure to access edge costs
* @param succ	: pointer to the successor array to refine with 2opt
* @result 1 if a move was performed successfully, 0 if there are no moves to perform
*/
int two_optimality_move(instance* inst, int* succ);

/*!
* @function two_optimality_move_tabu
* @abstract performs the tabu variant of the two optimality move on the provided solution, modifying the pair of edges
*			with the smallest cost change that aren't marked as tabu
* @param inst		: pointer to the instance structure to access edge costs
* @param succ		: pointer to the successor array to modify with 2opt
* @param tenure		: value of the tenure for the tabu list
* @param tabu_list	: pointer to the tabu list array
* @param curr_iter	: the current number of iterations of tabu search
*/
void two_optimality_move_tabu(instance* inst, int* succ, int tenure, int* tabu_list, int curr_iter);

/*!
* @function tabu_search
* @abstract computes the tabu search metaheuristic until the time limit is reached
* @param inst		: pointer to the instance structure to access edge costs
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @result the cost of the solution found by the algorithm
*/
double tabu_search(instance* inst, int* best_succ);

/*!
* @function five_opt_kick
* @abstract changes 5 random edges in the solution represented by the successor array, preserving a feasible cycle
* @param succ		: pointer to the successor array to modify with the 5opt kick
* @param nnodes		: length of the successor array
*/
void five_opt_kick(int* succ, int nnodes);

/*!
* @function vns
* @abstract computes the variable neighborhood search metaheuristic until the time limit is reached
* @param inst		: pointer to the instance structure to access edge costs
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @result the cost of the solution found by the algorithm
*/
double vns(instance* inst, int* best_succ);

/*!
* @function simulated_annealing
* @abstract computes the sinulated annealing metaheuristic until the time limit is reached
* @param inst		: pointer to the instance structure to access edge costs
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @param Tmax		: maximum temperature for the 2-opt move prbobaility computation
* @result the cost of the solution found by the algorithm
*/
double simulated_annealing(instance* inst, int* best_succ, double Tmax);

/*!
* @function create_child
* @abstract creates a feasible child solution out of 2 parent solutions and refines it with 2opt
* @param inst		: pointer to the instance structure to access edge costs
* @param sol		: pointer to the solution structure array where all solutions are saved
* @param parent_1	: index of the first parent in the solution array
* @param parent_2	: index of the second parent in the solution array
* @param child		: pointer to the solution structure representing the child
*/
void create_child(instance* inst, solution* sol, int parent_1, int parent_2, solution* child);

/*!
* @function mutate
* @abstract mutates a feasible solution from an existing one
* @param inst		: pointer to the instance structure to access edge costs
* @param sol		: pointer to the solution structure array where all solutions are saved
* @param parent		: index of the parent in the solution array
* @param child		: pointer to the solution structure representing the child
* @param n			: number of mutations to introduce
*/
void mutate(instance* inst, solution* sol, int parent, solution* child, int n);

/*!
* @function genetic
* @abstract computes the genetic search metaheuristic until the time limit is reached
* @param inst		: pointer to the instance structure to access edge costs
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @param num_sol	: size of the population
* @result the cost of the solution found by the algorithm
*/
double genetic(instance* inst, int* best_succ, int num_sol);

//CPLEX

/*!
* @function TSPopt_Benders
* @abstract computes the Benders' Loop solution until the time limit is reached
* @param inst		: pointer to the instance structure to solve
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @result the cost of the solution found by the algorithm
*/
double TSPopt_Benders(instance* inst, int* best_succ);

/*!
* @function TSPopt_branch_and_cut
* @abstract computes the CPLEX callback solution until the time limit is reached
* @param inst		: pointer to the instance structure to solve
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @result the cost of the solution found by the algorithm
*/
double TSPopt_branch_and_cut(instance* inst, int* best_succ);

/*!
* @function my_callback
* @abstract patches up the solution in the current b&c tree node and posts it to CPLEX to potentially improve bounds
* @param context		: CPXCALLBACKCONTEXTptr
* @param contextid		: CPXLONG
* @param userhandle		: void*
*/
static int CPXPUBLIC my_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void* userhandle);

/*!
* @function add_SEC
* @abstract adds a subtour elimination constraint to the CPLEX model for every connected component in the graph
* @param inst		: pointer to the instance structure in question
* @param env		: CPXENVptr
* @param lp			: CPXLPptr
* @param comp		: pointer to array marking what connected component each node is part of
* @param succ		: pointer to array storing the solution
* @param mycomp		: total number of conected components
* @param context	: CPXCALLBACKCONTEXTptr for use in callback, NULL otherwise
*/
void add_SEC(instance* inst, CPXENVptr env, CPXLPptr lp, int* comp, int* succ, int mycomp, CPXCALLBACKCONTEXTptr context);

/*!
* @function xpos
* @abstract returns the x position in the LP matrix of a n edge (i,j) in the graph
* @param i			: extreme point #1
* @param j			: extreme point #2
* @param inst		: pointer to the instance structure in question
* @result the x position of edge (i,j) in the LP matrix
*/
int xpos(int i, int j, instance* inst);

/*!
* @function build_model
* @abstract translates a graph into a CPLEX model
* @param inst		: pointer to the instance structure in question
* @param env		: CPXENVptr
* @param lp			: CPXLPptr
*/
void build_model(instance* inst, CPXENVptr env, CPXLPptr lp);

/*!
* @function build_sol
* @abstract translates a CPLEX solution into a successor vector format, counting the separate connected components
* @param xstar		: array containing the edges selected by the CPLEX solution
* @param inst		: pointer to the instance structure in question
* @param succ		: pointer to the array where the solution will be saved
* @param comp		: pointer to the array where the component idices will be saved
* @param ncomp		: pointer to the integer where the number of connected components will be saved
*/
void build_sol(const double* xstar, instance* inst, int* succ, int* comp, int* ncomp);

/*!
* @function	patchingHeuristic
* @abstract patches up a solution made up of different connected components with the least costly move, then 2-opts the final result, also adding SECs to the model
* @param inst		: pointer to the instance structure in question
* @param env		: CPXENVptr for adding SECs
* @param lp			: CPXLPptr for adding SECs
* @param succ		: pointer to the successor array containing the solution to patch up
* @param comp		: pointer to the array containing the component idices
* @param ncomp		: number of connected components in the solution before patching
*/
void patchingHeuristic(instance* inst, CPXENVptr env, CPXLPptr lp, int* succ, int* comp, const int ncomp);

/*!
* @function hard_fixing
* @abstract computes the Hard Fixing Matheuristic
* @param inst		: pointer to the instance structure to solve
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @param pfix		: probability of randomly fixing an edge of an intermediate heuristic solution
* @param cplex_tl	: time limit for the short CPLEX solver runs
* @result the cost of the solution found by the algorithm
*/
double hard_fixing(instance* inst, int* best_succ, double pfix, double cplex_tl);

/*!
* @function local_branching
* @abstract computes the Local Branching Matheuristic
* @param inst		: pointer to the instance structure to solve
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @param k			: size of the solution neighborhood to search in (number of different edges)
* @param cplex_tl	: time limit for the short CPLEX solver runs
* @result the cost of the solution found by the algorithm
*/
double local_branching(instance* inst, int* best_succ, int k, double cplex_tl);

/*!
* @function rins_polishing
* @abstract computes the CPLEX solution with RINS and Polishing heuristics
* @param inst		: pointer to the instance structure to solve
* @param best_succ	: pointer to the successor array where the best solution will be saved
* @param frequency	: frequency to use for RINS
* @param t_p		: time spent polishing during the algorithm
* @result the cost of the solution found by the algorithm
*/
double rins_polishing(instance* inst, int* best_succ, int frequency, int t_p);