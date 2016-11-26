#include <cstdio>
#include "treenode.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <queue>
#include <limits>
#include <bitset>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>

#include "densemap.h"
#include "netnode.h"
#include "example.h"

// Manually check some derivatives.

std::vector<TreeNode> geneNodes;
auto gene = createGene(geneNodes);

auto taxa = getTaxa(gene);
auto events = getEvents(gene, taxa);
std::vector<NetNode> results;
auto species = createSpeciesWithIntro(results, 0.0);

int main() {


	std::vector<double> result;

	calcProbability(species, gene, &result);
	for (unsigned int i = 0; i < result.size() ;i++) {
		printf("%d %g\n", i, result[i]);
	}

	double dx = 0.0000000001;
	std::vector<NetNode> temp1;
	std::vector<NetNode> temp2;
	double blah = calcProbability(createSpeciesWithIntro(temp1, dx), gene, nullptr) - calcProbability(createSpeciesWithIntro(temp1, -dx), gene, nullptr);

	printf("Prob: %g\n", (blah/(2 * dx)));

}