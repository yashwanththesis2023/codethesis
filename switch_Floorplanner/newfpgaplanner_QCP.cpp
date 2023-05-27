//============================================================================
// Name        : fpga-floorplanner.cpp
// Author      : Paul & Zakarya & Yashwanth 
// Version     :
// Copyright   :
// Description : exhaustive search using lp_solve to check all floorplans
//             : of rectangular soft module blocks
//============================================================================

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <climits>
#include <chrono>
#include <map>
#include <fstream>
#include <sstream>
#include <cstdlib>


#ifdef GUROBI_USE
#include "lib/gurobi_c++.h"
#endif


#ifdef LPSOLVE_USE
#include "lib/lp_lib.h"
#endif

#ifdef CPLEX_USE
#include <ilcplex/ilocplex.h>
ILOSTLBEGIN
typedef IloArray<IloNumVarArray> NumVarMatrix;
#endif

#define solver_use "define"


#if (CPLEX_USE==1 and LPSOLVE_USE+GUROBI_USE==0)
#define solver_use "cplex"
#endif

#if (GUROBI_USE==1 and LPSOLVE_USE+CPLEX_USE==0)
#define solver_use "gurobi"
#endif

#if (LPSOLVE_USE==1 and CPLEX_USE+GUROBI_USE==0)
#define solver_use "lpsolve"
#endif

#if (LPSOLVE_USE+CPLEX_USE+GUROBI_USE>=2)
#define solver_use "Undefined"
#endif

using namespace std;

/*
 Outline: <width> <height>

 <blockname> <area>
 <blockname> <area>
 <blockname> <x> <y> <block width> <block height>    <-- this is a preplaced block identifying FPGA areas that should not be used by the floorplanner

 */

typedef uint32_t pint;
typedef uint64_t xint;
typedef unsigned int uint;

int is_prime(xint);

inline int next_prime(pint p) {
	if (p == 2)
		return 3;
	for (p += 2; p > 1 && !is_prime(p); p += 2)
		;
	if (p == 1)
		return 0;
	return p;
}

int is_prime(xint n) {
#	define NCACHE 256
#	define S (sizeof(uint) * 2)
	static uint cache[NCACHE] = { 0 };

	pint p = 2;
	int ofs, bit = -1;

	if (n < NCACHE * S) {
		ofs = n / S;
		bit = 1 << ((n & (S - 1)) >> 1);
		if (cache[ofs] & bit)
			return 1;
	}

	do {
		if (n % p == 0)
			return 0;
		if (p * p > n)
			break;
	} while ((p = next_prime(p)));

	if (bit != -1)
		cache[ofs] |= bit;
	return 1;
}

int decompose(xint n, pint *out) {
	int i = 0;
	pint p = 2;
	while (n >= p * p) {
		while (n % p == 0) {
			out[i++] = p;
			n /= p;
		}
		if (!(p = next_prime(p)))
			break;
	}
	if (n > 1)
		out[i++] = n;
	return i;
}

int startswidth(const char *str, const char *prefix) {
	return strncmp(str, prefix, strlen(prefix)) == 0;
}

class block {

public:

	int area = 0, x = 0, y = 0, width = 0, height = 0, preplaced = 0;
	char *name = NULL;
	pint *primes = NULL;
	pint *exponents = NULL;
	pint *tmp_exponents = NULL;
	int primes_len = 0;

	block() {
	}

	~block() {
		if (name)
			free(name);
		if (primes)
			free(primes);
		if (exponents)
			free(exponents);
		if (tmp_exponents)
			free(tmp_exponents);
	}

	void prime_factor_decompose_area() {
		if (preplaced) {
			primes = exponents = NULL;
			primes_len = 0;
			return;
		}

		pint out[100];
		primes_len = decompose(area, out);
		pint lprimes[100];
		pint lexponents[100];

		lprimes[0] = out[0];
		lexponents[0] = 1;
		int j = 0;
		for (int i = 1; i < primes_len; i++)
			if (lprimes[j] == out[i])
				lexponents[j]++;
			else {
				lprimes[++j] = out[i];
				lexponents[j] = 1;
			}
		primes_len = j + 1;

		primes = (pint*) malloc(primes_len * sizeof(pint));
		memcpy(primes, lprimes, primes_len * sizeof(pint));

		exponents = (pint*) malloc(primes_len * sizeof(pint));
		memcpy(exponents, lexponents, primes_len * sizeof(pint));

		tmp_exponents = (pint*) malloc(primes_len * sizeof(pint));
	}

	void print() {
		printf("name: %s area: %d x: %d y: %d width: %d height: %d pre-placed: %d\nprime factor decomposition: ", name, area, x, y, width, height, preplaced);
		for (int i = 0; i < primes_len; i++)
			printf("%u^%u ", primes[i], exponents[i]);
		printf("\nprime factor decomposit[2]: ");
		for (int i = 0; i < primes_len; i++)
			printf("%u^%u ", primes[i], tmp_exponents[i]);
		printf("\n");
	}

	/**
	 * compute the width and height of the tmp_exponents configuration
	 * return the width
	 */
	int compute_width_height() {
		if (!preplaced) {
			width = 1;
			for (int i = 0; i < primes_len; i++)
				width *= (int) pow(primes[i], tmp_exponents[i]);
			height = area / width;
		}
		return width;
	}

	void print_block_configurations() {
		compute_width_height();
		printf("%s w=%d h=%d ", name, width, height);
	}
};

class fpga_floorplan {

public:

	int fpga_width;
	int fpga_height;
	int AreaBounding;
	vector<block*> *block_list;
    map<string, int> solver_map;

	fpga_floorplan()
    {
        fpga_width = fpga_height = 0;
        block_list = new vector<block *>();
        solver_map["cplex"] = 1;
        solver_map["CPLEX"] = 1;
        solver_map["gurobi"] = 2;
        solver_map["GUROBI"] = 2;
        solver_map["lp_solve"] = 3;
        solver_map["lpsolve"] = 3;
        solver_map["LPSOLVE"] = 3;
        solver_map["LP_SOLVE"] = 3;
    }

	~fpga_floorplan() {
		for (auto &v : *block_list)
			free(v);
		free(block_list);
	}

	void load(char *filename) {

		printf("in load");

		FILE * f;
		if (filename == NULL)
			f = stdin;
		else if (!(f = fopen(filename, "r"))) {
			printf("can't open %s\n", filename);
			exit(1);
		}

		char buffer1[25501];
		char *line = buffer1;
		char *ptr;

		while (fgets(line, 25500, f)) {

			if (startswidth(line, "#"))
				continue;
			else if (strlen(line) <= 1)
				continue;
			else if (startswidth(line, "Outline:")) {
				ptr = strtok(line, " ");
				ptr = strtok(NULL, " ");
				fpga_width = atoi(ptr);
				ptr = strtok(NULL, " ");
				fpga_height = atoi(ptr);
			} else { // read block name and 1 to 4 integers
				block *b = new block();

				// copy block name; malloc space for it
				char *name = strtok(line, " ");
				int name_size = strlen(name);
				b->name = (char*) malloc((name_size + 1) * sizeof(char));
				strcpy(b->name, name);

				// continue with area
				ptr = strtok(NULL, " ");
				b->area = atoi(ptr);

				// regular module with area specified only?
				if ((ptr = strtok(NULL, " "))) {
					// pre-placed module
					b->x = b->area;
					b->area = 0;
					b->preplaced = 1;
					b->y = atoi(ptr);
					ptr = strtok(NULL, " ");
					b->width = atoi(ptr);
					ptr = strtok(NULL, " ");
					b->height = atoi(ptr);
					b->area = b->width * b->height;
				}
				block_list->push_back(b);

			}
		} // end while
		fclose(f);
	}

	/**
	 * prime factor decomposition of all block areas
	 */
	void decompose() {
		for (auto const& v : *block_list)
			v->prime_factor_decompose_area();
	}

	/**
	 * go through all blocks and within the blocks through all prime factors
	 * and create all possible block <height, width> configurations
	 * all symmetrical configurations can also be created, i.e. <w,h> and <h,w>
	 * but this is suppressed at the moment
	 */
	bool generate_exhaustive(int block_index, int prime_index, fpga_floorplan* fp, int max_unchanged_iterations,int area) {

		static int unchanged_iterations = 0;  // Counter variable for unchanged iterations
		
		// Set all block configurations? Then print and return
		if (block_index >= fp->block_list->size()) {
			double objval = (fp->*fp->solver)();
			printf("current floorplan area: %f\n", objval);
			static int best_height = INT_MAX;

			int new_height = (int)objval;
			if (new_height < best_height) {
				best_height = min(best_height, new_height);
				fp->AreaBounding = best_height;
				cout << "Area bounding in gen" << AreaBounding << endl;
				unchanged_iterations = 0;
			} else if (best_height < area) {
				return true;  // Exit the function and continue with the main code
			} else {
				unchanged_iterations++;
			}

			cout << "unchanged_iterations = " << unchanged_iterations << endl;
			printf("best bounding area found: %d\n", best_height);

			if (unchanged_iterations >= max_unchanged_iterations) {
				printf("Exiting the process due to no improvement in best height.\n");
				return true;  // Exit the function and continue with the main code
			}

			return false;
		}

		block* b = fp->block_list->at(block_index);

		// If the block is pre-placed, then only one geometry is possible
		if (b->preplaced) {
			if (generate_exhaustive(block_index + 1, 0, fp, max_unchanged_iterations,area))
				return true;  // Exit the function and continue with the main code
			return false;
		}

		// If all prime exponents in the current block set, continue with prime 0 from the next block
		if (prime_index >= b->primes_len) {
			if (!(b->compute_width_height() > (((int)sqrt(b->area)) + 1))) {
				// Skip symmetric case
				if (generate_exhaustive(block_index + 1, 0, fp, max_unchanged_iterations,area))
					return true;  // Exit the function and continue with the main code
			}
			return false;
    }

    // If the block is not pre-placed and not all prime exponents are set
    for (int i = 0; i <= b->exponents[prime_index]; i++) {
        b->tmp_exponents[prime_index] = i;
        if (generate_exhaustive(block_index, prime_index + 1, fp, max_unchanged_iterations,area))
            return true;  // Exit the function and continue with the main code
    }

    return false;
	}
    double getBoundingRectAreaNEW() {
    int min_x = INT_MAX, max_x = INT_MIN, min_y = INT_MAX, max_y = INT_MIN;
    for (auto const &v : *block_list) { 
        //cout<< endl<< v->name <<endl;        
        int x = v->x, y = v->y, width = v->width, height = v->height;
        //cout<<"updated in mbr"<<endl;
        //printf("%s, x = %d, y  = %d, width = %d, height = %d\n",v->name,x,y,width,height);
        
        min_x = min(min_x, x);
        max_x = max(max_x, x + width);
        min_y = min(min_y, y);
        max_y = max(max_y, y + height);
    }
    cout<<endl;
    printf("min_x = %d, max_x  = %d, min_y = %d, max_y = %d\n",min_x,max_x,min_y,max_y);
    return (max_x - min_x) * (max_y - min_y);
    }


	void print() {
		printf("Outline/FPGA width: %d height: %d\n", fpga_width, fpga_height);
		for (auto const& v : *block_list) {
			v->print();
			printf("\n");
		}
	}

	void print_block_configurations() {
		for (auto const& v : *block_list)
			v->print_block_configurations();
		printf("\n");
	}

	void print_lp(FILE *f) {
		block *b;

		fprintf(f, "// %10s %10s %10s %10s %10s %10s\n", "name", "idx", "width", "height", "x", "y");
		// print list of comments with block names and geometries
		for (int i = 0; i < block_list->size(); i++) {
			b = block_list->at(i);
			if (b->preplaced)
				fprintf(f, "// %10s %10d %10d %10d %10d %10d\n", b->name, i, b->width, b->height, b->x, b->y);
			else
				fprintf(f, "// %10s %10d %10d %10d\n", b->name, i, b->width, b->height);
		}

		fprintf(f, "\nmin: y;\n\n");

		// compute maximum over widths and heights of all modules
		// to have an upper bound on FPGA area height
		int M = 0;
		for (auto const& v : *block_list)
			M += max(v->width, v->height);

		int eq = 1; // equation number
		int hi, hj, wi, wj;
		// xi + hi ri + wi - wi ri <= W for 1 <= i <= n
		// yi + wi ri + hi - hi ri <= y for 1 <= i <= n
		for (int i = 1; i <= block_list->size(); i++) {
			hi = block_list->at(i - 1)->height;
			wi = block_list->at(i - 1)->width;
			fprintf(f, "c%d:\t x%-3d + %-3d r%-3d + %-3d - %-3d r%-3d <= %-3d;\n", eq++, i, hi, i, wi, wi, i, fpga_width);
			fprintf(f, "c%d:\t y%-3d + %-3d r%-3d + %-3d - %-3d r%-3d <= y; \n", eq++, i, wi, i, hi, hi, i);
		}

		// xi + hi ri + wi - wi ri <= xj +       M pij + M qij for 1 <= i < j <= n
		// yi + wi ri + hi - hi ri <= yj +   M + M pij - M qij for 1 <= i < j <= n
		// xi - hj rj - wj + wj rj >= xj -   M + M pij - M qij for 1 <= i < j <= n
		// yi - wj rj - hj + hj rj >= yj - 2 M + M pij + M qij for 1 <= i < j <= n
		fprintf(f, "\n");
		for (int i = 1; i <= block_list->size(); i++) {
			for (int j = i + 1; j <= block_list->size(); j++) {
				hi = block_list->at(i - 1)->height;
				hj = block_list->at(j - 1)->height;
				wi = block_list->at(i - 1)->width;
				wj = block_list->at(j - 1)->width;
				fprintf(f, "c%d:\t x%-3d + %-3d r%-3d + %-3d - %-3d r%-3d <= x%-3d +       %-3d p%d%d + %-3d q%d%d;\n", eq++, i, hi, i, wi, wi, i, j, M, i, j, M, i, j);
				fprintf(f, "c%d:\t y%-3d + %-3d r%-3d + %-3d - %-3d r%-3d <= y%-3d + %-3d + %-3d p%d%d - %-3d q%d%d;\n", eq++, i, wi, i, hi, hi, i, j, M, M, i, j, M, i, j);
				fprintf(f, "c%d:\t x%-3d - %-3d r%-3d - %-3d + %-3d r%-3d >= x%-3d - %-3d + %-3d p%d%d - %-3d q%d%d;\n", eq++, i, hj, j, wj, wj, j, j, M, M, i, j, M, i, j);
				fprintf(f, "c%d:\t y%-3d - %-3d r%-3d - %-3d + %-3d r%-3d >= y%-3d - %-3d + %-3d p%d%d + %-3d q%d%d;\n", eq++, i, wj, j, hj, hj, j, j, (2 * M), M, i, j, M, i, j);
				fprintf(f, "\n");
			}
		}

		// print rotational constraints
		for (int i = 1; i <= block_list->size(); i++)
			fprintf(f, "r%d <= 1;\n", i);
		for (int i = 1; i <= block_list->size(); i++)
			for (int j = i + 1; j <= block_list->size(); j++)
				fprintf(f, "p%d%d <= 1;\nq%d%d <= 1;\n", i, j, i, j);

		// print x y pairs of preplaced blocks and overwrite in this case rotational constraints to non-rotational
		fprintf(f, "\n");
		for (int i = 1; i <= block_list->size(); i++)
			if (block_list->at(i - 1)->preplaced) {
				fprintf(f, "x%d = %d;\n", i, block_list->at(i - 1)->x);
				fprintf(f, "y%d = %d;\n", i, block_list->at(i - 1)->y);
				fprintf(f, "r%d = 0;\n", i);
			}

		fprintf(f, "\nint ");
		for (int i = 1; i <= block_list->size(); i++)
			fprintf(f, "r%d, ", i);
		for (int i = 1; i <= block_list->size(); i++)
			for (int j = i + 1; j <= block_list->size(); j++)
				fprintf(f, "p%d%d, q%d%d%c ", i, j, i, j, (i == block_list->size() - 1) && (j == block_list->size()) ? ';' : ',');
		fprintf(f, "\n");

	}
#ifdef GUROBI_USE
    double gurobi()
    {
        block *b;
        GRBEnv env = GRBEnv();
        GRBModel m1 = GRBModel(env);
        // compute maximum over widths and heights of all modules
        // to have an upper bound on FPGA area height
        int M = max(fpga_width,fpga_height);
        double bounding_area = 0.0;
        GRBVar* vars = NULL;
        // for (auto const &v : *block_list)
        //     M += max(v->width, v->height);

        int eq1 = 1; // equation number
        int hi, hj, wi, wj;
        
        int sizeB = block_list->size() + 1;
        GRBVar yh, x[sizeB], y[sizeB], r[sizeB], p[sizeB][sizeB], q[sizeB][sizeB];
        int j;
        yh = m1.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y");
        for (int i = 1; i <= block_list->size(); i++)
        {
            x[i] = m1.addVar(0, GRB_INFINITY, 0, GRB_INTEGER, "x" + to_string(i));
            r[i] = m1.addVar(0, 1, 0, GRB_INTEGER, "r" + to_string(i));
            y[i] = m1.addVar(0, GRB_INFINITY, 0, GRB_INTEGER, "y" + to_string(i));
            for (j = i + 1; j <= block_list->size(); j++)
            {
                p[i][j] = m1.addVar(0, 1, 0, GRB_BINARY, "p" + to_string(i) + to_string(j));
                q[i][j] = m1.addVar(0, 1, 0, GRB_BINARY, "q" + to_string(i) + to_string(j));
            }
        }

        // Subject To :
        for (int i = 1; i <= block_list->size(); i++)
        {
            hi = block_list->at(i - 1)->height;
            wi = block_list->at(i - 1)->width;
            m1.addConstr(x[i] + hi * r[i] - wi * r[i] <= fpga_width - wi, "c" + to_string(eq1++));
            m1.addConstr(y[i] - yh <= -hi, "c" + to_string(eq1++));
        }
        
        for (int i = 1; i <= block_list->size(); i++)
        {
            for (int j = i + 1; j <= block_list->size(); j++)
            {
                hi = block_list->at(i - 1)->height;
                hj = block_list->at(j - 1)->height;
                wi = block_list->at(i - 1)->width;
                wj = block_list->at(j - 1)->width;
                m1.addConstr(x[i] + hi * r[i] - wi * r[i] - x[j] - M * p[i][j] - M * q[i][j] <= -wi, "c" + to_string(eq1++));
                m1.addConstr(y[i] + wi * r[i] - hi * r[i] - y[j] - M * p[i][j] + M * q[i][j] <= M - hi, "c" + to_string(eq1++));
                m1.addConstr(x[i] - hj * r[j] + wj * r[j] - x[j] - M * p[i][j] + M * q[i][j] >= -M + wj, "c" + to_string(eq1++));
                m1.addConstr(y[i] - wj * r[j] + hj * r[j] - y[j] - M * p[i][j] - M * q[i][j] >= -2*M + hj, "c" + to_string(eq1++));
            }
        }

        // print x y pairs of preplaced blocks and overwrite in this case rotational constraints to non-rotational
        for (int i = 1; i <= block_list->size(); i++)
            if (block_list->at(i - 1)->preplaced)
            {
                x[i].set(GRB_DoubleAttr_LB, block_list->at(i - 1)->x);
                x[i].set(GRB_DoubleAttr_UB, block_list->at(i - 1)->x);
                y[i].set(GRB_DoubleAttr_LB, block_list->at(i - 1)->y);
                y[i].set(GRB_DoubleAttr_UB, block_list->at(i - 1)->y);
                r[i].set(GRB_DoubleAttr_LB, 0);
                r[i].set(GRB_DoubleAttr_UB, 0);
            }
        m1.setObjective(yh + 0 * x[1]);

        m1.update();

        for (int i = 1; i <= block_list->size(); i++)
         {
             hi = block_list->at(i - 1)->height;
             wi = block_list->at(i - 1)->width;
             //printf("%s: w=%d h=%d \n", block_list->at(i - 1)->name, wi,hi);
         }
        
        m1.write("model_gurobi.lp");
        m1.optimize();
        if (m1.get(GRB_IntAttr_Status) == 2){
        //m1.computeIIS();
        // if(m1.get(GRB_IntAttr_Status) != 2)
        // {return 1000;}
        m1.write("model_gurobi.sol");

        cout << "--------------" << __FUNCTION__ << "---------------" << endl;
        for (int i = 1; i <= block_list->size(); i++){
            block_list->at(i - 1)->x = x[i].get(GRB_DoubleAttr_X);
            block_list->at(i - 1)->y = y[i].get(GRB_DoubleAttr_X);
            // cout << "hi" << block_list->at(i - 1)->height << endl;
            // cout << "wi" << block_list->at(i - 1)->width << endl;
        }
        
        bounding_area = getBoundingRectAreaNEW();
        cout<< "current bounding area  " << bounding_area << endl;
        cout<< "height  "<< m1.get(GRB_DoubleAttr_ObjVal)<<endl;
        //print_block_configurations();
        return bounding_area;}
        else{return 100000;}

    }
#endif

void solverCheck(string name)
    {
        switch(solver_map[name])
        {
#ifdef CPLEX_USE
        case 1:
            this->solver = &fpga_floorplan::cplex;
            break;
#endif
#ifdef GUROBI_USE
        case 2:
            this->solver = &fpga_floorplan::gurobi;
            break;
#endif
#ifdef LPSOLVE_USE 
        case 3:
            this->solver = &fpga_floorplan::lp_solve;
            break;
#endif
        default:
            cout << "solver not found or not defined" << endl;
            exit(0);
        }

    }

    double (fpga_floorplan::*solver)();

int readAreaBounding(char *filename){
	
	std::ifstream file("final_MBLA_SA_optimization_req_results.csv"); 

    if (file.is_open()) {
        std::string line;
        std::getline(file, line); // Skip the header line if needed
		std::string filename1=filename;
        int AreaBounding_value;
        while (std::getline(file, line)) {
            std::istringstream iss(line);
            std::string name,AreaBounding;
			//printf("after reading the file\n");

            // Read the values using a delimiter (e.g., comma)
            std::getline(iss, name, ';'); // Read the first value
			//cout << name << endl;
			//cout << filename << endl;
			cout <<"name " << filename1 << endl;
            if(filename1.find(name) != std::string::npos){
            std::getline(iss, AreaBounding, ';'); // Read the second value
            AreaBounding_value = std::stoi(AreaBounding);
            std::cout << "AreaBounding value: " << AreaBounding_value << std::endl;
            return AreaBounding_value;
            }
            
        }	
}
}

};


int main(int argcc, char** argv) {
    int value = 0;
    int width = 0;

    // Record start time
    auto start = chrono::high_resolution_clock::now();

    fpga_floorplan* fp = new fpga_floorplan();

    int area;

    if (argcc == 1)
        fp->load(NULL);
    else
        fp->load(argv[1]);

    area = fp->readAreaBounding(argv[1]);
	cout << area << endl;
    printf("loading done\n");

    if (argcc == 3)
        fp->solverCheck(argv[2]);
    else
        fp->solverCheck(solver_use);

    fp->decompose();

    printf("decomposition done\n");

    int max_unchanged_iterations = 30;

    while (true) {
        bool stop_generating = fp->generate_exhaustive(0, 0, fp, max_unchanged_iterations, area);
		if (stop_generating) {
			// Continue with the remaining code in the loop
			value = fp->AreaBounding;
			width = fp->fpga_width;
			printf("bounding area = %d \n", value);

			if (area < fp->AreaBounding) {
				fp->fpga_width = floor(area / fp->AreaBounding);
				printf("width = %d\n", width);
			} else {
				break;  // Exit the loop and continue with the main code
			}
		}
	}
    

    // Record end time
    auto finish = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed = finish - start;

    cout << "Elapsed time: " << elapsed.count() << endl;

    return 0;
}