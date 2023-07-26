#include "tsp.h"
void free_instance(instance* inst)
{
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->best_sol);
	free(inst->cst);
}

void print_error(const char* err) 
{ 
	printf("\n\n ERROR: %s \n\n", err); 
	fflush(NULL); 
	exit(1); 
}

void read_input(instance* inst) // simplified TSP parser, not all SECTIONs detected  
{

	FILE* fin = fopen(inst->input_file, "r");
	if (fin == NULL) printf(" input file not found!\n");

	inst->nnodes = -1;

	char line[180];
	char* par_name;
	char* token1;
	char* token2;

	int active_section = 0; // =1 NODE_COORD_SECTION 

	int do_print = (VERBOSE >= 1000);

	while (fgets(line, sizeof(line), fin) != NULL)
	{
		if (VERBOSE >= 2000) { printf("%s", line); fflush(NULL); }
		if (strlen(line) <= 1) continue; // skip empty lines
		par_name = strtok(line, " :");
		if (VERBOSE >= 3000) { printf("parameter \"%s\" ", par_name); fflush(NULL); }

		if (strncmp(par_name, "NAME", 4) == 0)
		{
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "COMMENT", 7) == 0)
		{
			active_section = 0;
			token1 = strtok(NULL, "");
			if (VERBOSE >= 10) printf(" ... solving instance %s\n\n", token1);
			continue;
		}

		if (strncmp(par_name, "TYPE", 4) == 0)
		{
			token1 = strtok(NULL, " :");
			if (strncmp(token1, "TSP", 3) != 0) print_error(" format error:  only TYPE == TSP implemented so far!!!!!!");
			active_section = 0;
			continue;
		}


		if (strncmp(par_name, "DIMENSION", 9) == 0)
		{
			if (inst->nnodes >= 0) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if (do_print) printf(" ... nnodes %d\n", inst->nnodes);
			inst->xcoord = (double*) malloc(inst->nnodes * sizeof(double));
			inst->ycoord = (double*) malloc(inst->nnodes * sizeof(double));
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "EDGE_WEIGHT_TYPE", 16) == 0)
		{
			token1 = strtok(NULL, " :");
			if ((strncmp(token1, "EUC_2D", 6) != 0) && (strncmp(token1, "ATT", 3) != 0))
				print_error(" format error:  only EDGE_WEIGHT_TYPE == EUC_2D or EDGE_WEIGHT_TYPE == ATT implemented so far!!!!!!");
			active_section = 0;
			continue;
		}

		if (strncmp(par_name, "NODE_COORD_SECTION", 18) == 0)
		{
			if (inst->nnodes <= 0) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;
			continue;
		}

		if (strncmp(par_name, "EOF", 3) == 0)
		{
			active_section = 0;
			break;
		}


		if (active_section == 1) // within NODE_COORD_SECTION
		{
			int i = atoi(par_name) - 1;
			if (i < 0 || i >= inst->nnodes) print_error(" ... unknown node in NODE_COORD_SECTION section");
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if (do_print) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i + 1, inst->xcoord[i], inst->ycoord[i]);
			continue;
		}

		printf(" final active section %d\n", active_section);
		print_error(" ... wrong format for the curren simplified parser!!!!!!!!!");

	}

	fclose(fin);

}

// read and analyze the arguments in the command line
void parse_command_line(int argc, char** argv, instance* inst)
{
	if (VERBOSE >= 100) printf(" running %s with %d parameters \n", argv[0], argc - 1);

	// default   
	strcpy(inst->input_file, "NULL");
	inst->randomseed = 0;
	inst->mode = 1;
	inst->time_limit = 5;
	inst->nnodes = 1000;
	inst->integer_costs = 0;
	inst->p = 0;


	int help = 0; if (argc < 1) help = 1;
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-file") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-input") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 			// input file
		if (strcmp(argv[i], "-f") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 				// input file
		if (strcmp(argv[i], "-seed") == 0) { inst->randomseed = abs(atoi(argv[++i])); continue; } 		// random seed
		if (strcmp(argv[i], "-mode") == 0) { inst->mode = atoi(argv[++i]); continue; } 					// solver mode
		if (strcmp(argv[i], "-tl") == 0) { inst->time_limit = atof(argv[++i]); continue; } 				// time limit
		if (strcmp(argv[i], "-n") == 0) { inst->nnodes = atoi(argv[++i]); continue; } 					// number of nodes of random instance
		if (strcmp(argv[i], "-i_c") == 0) { inst->integer_costs = 1; continue; } 						// set edge cost as integer
		if (strcmp(argv[i], "-p") == 0) { inst->p = 1; continue; } 										// set edge cost as integer
		if (strcmp(argv[i], "-help") == 0) { help = 1; continue; }											// program help
		/*
		mode = 1: GRASP Nearest neighbor heuristic + 2opt until no move is possible
		mode = 2: farthest couple start GRASP Extra mileage heuristic + 2opt until no move is possible
		mode = 3: multi start GRASP Extra mileage heuristic + 2opt until no move is possible (COMPUTATIONALLY INTENSIVE)
		mode = 4: Tabu search metaheuristic starting from GRASP Nearest neighbor solution
		mode = 5: VNS metaheuristic starting from GRASP Nearest neighbor solution
		mode = 6: Genetic metaheuristic from a population of 1000 solutions
		mode = 7: Cplex exact solution implementation 
		mode = 8: Cplex exact solution implementation with branch & cut
		mode = 9: Hard Fixing Matheuristic
		mode = 10: Local Branching Matheuristic
		mode = 11: Simulated Annealing heuristic
		mode = 12: Rins with solution polishing
		*/

		if (help || (VERBOSE >= 50))		// print current parameters
		{
			printf("\n\navailable parameters (vers. 27-aug-2022) --------------------------------------------------\n");
			printf("-file %s\n", inst->input_file);
			printf("-tl %lf\n", inst->time_limit);
			printf("-seed %d\n", inst->randomseed);
			printf("-mode %d\n", inst->mode);
			printf("-n %d\n", inst->nnodes);
			printf("-i_c %d\n", inst->integer_costs);
			printf("-p %d\n", inst->p);
			printf("\nenter -help or --help for help\n");
			printf("----------------------------------------------------------------------------------------------\n\n");
		}

		if (help) exit(1);
	}
}

//print out the number of nodes and their coordinates
void print_instance(instance* inst) 
{
	printf("Number of nodes: %d\n", inst->nnodes);
	int i = 0;
	for (i; i < inst->nnodes; i++) {
		printf("%d	%lf	%lf\n", i+1, inst->xcoord[i], inst->ycoord[i]);
	} //for
}

//generates a random TSP instance with 5 to 1000 nodes in a 1000x1000 square.
void random_instance(instance* inst)
{
	inst->xcoord = (double*)malloc(inst->nnodes * sizeof(double));
	inst->ycoord = (double*)malloc(inst->nnodes * sizeof(double));
	for (int i = 0; i < inst->nnodes; i++)
	{
		inst->xcoord[i] = 1000 * random01();
		inst->ycoord[i] = 1000 * random01();
	}
}

int main(int argc, char** argv) 
{
	if (argc < 4) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }
	if (VERBOSE >= 4) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	instance inst;
	solution current_sol;

	parse_command_line(argc, argv, &inst);

	if (strcmp(inst.input_file, "NULL") == 0) {
		random_instance(&inst);
	}
	else {
		read_input(&inst);
	}
	
	//initalize radnom number generator
	srand(7645321 + inst.randomseed);
	for (int k = 0; k < 1000; k++)
	{
		random01();
	}

	// preprocessing of the instance and current solution
	preprocessing(&inst, &current_sol);

	if (inst.p == 1)
	{
		//print the nodes of the instance on screen
		print_instance(&inst);
		plot_solution(&inst, 0);
	}
	if (inst.mode == 1)
	{
		//compute GRASP nearest neighbor heuristic for every starting node
		current_sol.cost = multistart_grasp_nearest_neighbor(&inst, current_sol.succ, 0.1, 1);
	}
	else if (inst.mode == 2)
	{
		//compute extra mileage heuristic starting from the farthest node couple
		current_sol.cost = farthest_grasp_extra_mileage(&inst, current_sol.succ, 0.1, 1);
	}
	else if (inst.mode == 3)
	{
		//compute extra mileage heuristic starting from each node and the farthest one from it
		current_sol.cost = multistart_grasp_extra_mileage(&inst, current_sol.succ, 0.1, 1);
	}
	
	else if (inst.mode == 4)
	{
		//compute the tabu search solution from a grasp nearest neighbor starting point
		current_sol.cost = tabu_search(&inst, current_sol.succ);
	}

	else if (inst.mode == 5)
	{
		//compute the variable neighborhood search solution from a grasp nearest neighbor starting point
		current_sol.cost = vns(&inst, current_sol.succ);
	}

	else if (inst.mode == 6)
	{
		//compute the genetic search solution
		current_sol.cost = genetic(&inst, current_sol.succ, 1000);
	}
	else if (inst.mode == 7)
	{
		//compute cplex exact solution with Benders loop method
		current_sol.cost = TSPopt_Benders(&inst, current_sol.succ);
	}
	else if (inst.mode == 8)
	{
		//compute cplex exact solution with branch & cut
		current_sol.cost = TSPopt_branch_and_cut(&inst, current_sol.succ);
	}
	else if (inst.mode == 9)
	{
		//compute hard fixing matheuristic solution
		current_sol.cost = hard_fixing(&inst, current_sol.succ, 0.3, 30);

	}
	else if (inst.mode == 10)
	{
		//compute local branching matheuristic solution
		current_sol.cost = local_branching(&inst, current_sol.succ, 10, 30);
	}

	else if (inst.mode == 11)
	{
		//compute local branching matheuristic solution
		current_sol.cost = simulated_annealing(&inst, current_sol.succ, 10);
	}

	else if (inst.mode == 12)
	{
		//compute local branching matheuristic solution
		current_sol.cost = rins_polishing(&inst, current_sol.succ, 30, 15, 30);
	}

	// update the incumbent and plot the solution
	postprocessing(&inst, &current_sol);

	// free the memory
	free_instance(&inst);
	free(current_sol.succ);
	return 0;
}