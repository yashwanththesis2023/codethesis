//============================================================================
// Name        : fpga-floorplanner.cpp
// Author      : Paul
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
	vector<block*> *block_list;

	fpga_floorplan() {
		fpga_width = fpga_height = 0;
		block_list = new vector<block*>();
	}

	~fpga_floorplan() {
		for (auto &v : *block_list)
			free(v);
		free(block_list);
	}

	void load(char *filename) {

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
	void generate_exhaustive(int block_index, int prime_index) {

		// Set all block configurations? Then print and return
		if (block_index >= block_list->size()) {
		    print_block_configurations();

			// printf("\n start evaluating the following program\n\n");
			// print_lp(stdout);

			// printf("\nstarting lpsolve\n\n");
			// FILE *tmp = fopen("/tmp/lp.txt", "w");
			// print_lp(tmp);
			// fclose(tmp);

			// start lp_solve
			// FILE *pipe;
			// char buf[1000];
			// pipe = popen("/bin/lp_solve /tmp/lp.txt", "r");
			// print_lp(pipe);

			// // find best floorplan (width the smalles height)
			// static int best_height = INT_MAX;
			// char *ptr;
			// while (fgets(buf, sizeof buf, pipe))
			// 	if (startswidth(buf, "Value of objective function:")) {
			// 		printf("%s", buf);
			// 		ptr = strtok(buf, " ");
			// 		ptr = strtok(NULL, " ");
			// 		ptr = strtok(NULL, " ");
			// 		ptr = strtok(NULL, " ");
			// 		ptr = strtok(NULL, " ");
			// 		int new_height = atoi(ptr);
			// 		if (new_height < best_height) {
			// 			print_block_configurations();
			// 			best_height = min(best_height, atoi(ptr));
			// 		}
			// 	}
			// pclose(pipe);

			// printf("best floor plan found with the minimal height of %d\n", best_height);
			 return;
		}

		block *b = block_list->at(block_index);

		// is the block pre-placed, then only one geometry is possible
		if (b->preplaced) {
			generate_exhaustive(block_index + 1, 0);
			return;
		}

		// all prime exponents in current block set, continue with prime 0 from the next block
		if (prime_index >= b->primes_len) {
			if (!(b->compute_width_height() > (((int) sqrt(b->area)) + 1))) // skip symmetric case
				generate_exhaustive(block_index + 1, 0);
			return;
		}

		// not pre-placed module and not all prime exponents set
		for (int i = 0; i <= b->exponents[prime_index]; i++) {
			b->tmp_exponents[prime_index] = i;
			generate_exhaustive(block_index, prime_index + 1);
		}
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
};

int main(int argcc, char** argv) {

	fpga_floorplan *fp = new fpga_floorplan();
	if (argcc == 1)
		fp->load(NULL);
	else
		fp->load(argv[1]);

	printf("loading done\n");

	fp->decompose();

	printf("decomposition done\n");

	fp->generate_exhaustive(0, 0);

	// return OK
	return (0);

}
