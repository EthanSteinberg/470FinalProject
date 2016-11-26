#pragma once

#include <memory>
#include "netnode.h"
#include "treenode.h"

/**
 * Create a simple species network.
 */
inline NetNode& createSpecies(std::vector<NetNode>& results, double dx = 0) {
	results.reserve(13);

	results.emplace_back("A");
	NetNode& A = results.back();

	results.emplace_back("B");
	NetNode& B = results.back();

	results.emplace_back("C");
	NetNode& C = results.back();

	results.emplace_back("D");
	NetNode& D = results.back();

	results.emplace_back("E");
	NetNode& E = results.back();

	results.emplace_back("F");
	NetNode& F = results.back();

	results.emplace_back("G");
	NetNode& G = results.back();

	results.emplace_back("one", Edge<NetNode>(0, B, 1.0), Edge<NetNode>(1, C, 1.0));
	NetNode& one = results.back();

	results.emplace_back("two", Edge<NetNode>(2, A, 1.0), Edge<NetNode>(11, one, 1.0));
	NetNode& two = results.back();

	results.emplace_back("three", Edge<NetNode>(3, D, 1.0), Edge<NetNode>(10, E, 1.0));
	NetNode& three = results.back();

	results.emplace_back("four", Edge<NetNode>(4, F, 1.0), Edge<NetNode>(9, G, 1.0));
	NetNode& four = results.back();

	results.emplace_back("five", Edge<NetNode>(5, three, 1.0), Edge<NetNode>(8, four, 1.0));
	NetNode& five = results.back();

	results.emplace_back("six", Edge<NetNode>(6, two, 1.0 + dx), Edge<NetNode>(7, five, 1.0));
	NetNode& six = results.back();

	return six;
}


/**
 * Create a simple species network with a trivial introgression.
 */
inline NetNode& createSpeciesWithTrivialIntro(std::vector<NetNode>& results) {
	results.reserve(15);

	results.emplace_back("A");
	NetNode& A = results.back();

	results.emplace_back("B");
	NetNode& B = results.back();

	results.emplace_back("C");
	NetNode& C = results.back();

	results.emplace_back("D");
	NetNode& D = results.back();

	results.emplace_back("E");
	NetNode& E = results.back();

	results.emplace_back("F");
	NetNode& F = results.back();

	results.emplace_back("G");
	NetNode& G = results.back();

	results.emplace_back("one", Edge<NetNode>(0, B, 1.0), Edge<NetNode>(1, C, 1.0));
	NetNode& one = results.back();

	results.emplace_back("oneIntrogressed", Edge<NetNode>(2, one, 0.0), 1.0, 15);
	NetNode& oneIntrogressed = results.back();

	results.emplace_back("two", Edge<NetNode>(3, A, 1.0), Edge<NetNode>(5, oneIntrogressed, 0.0, EdgeType::LEFT));
	NetNode& two = results.back();

	results.emplace_back("three", Edge<NetNode>(4, D, 1.0), Edge<NetNode>(14, E, 1.0));
	NetNode& three = results.back();

	results.emplace_back("superThree", Edge<NetNode>(6, three, 1.0), Edge<NetNode>(7, oneIntrogressed, 0.0, EdgeType::RIGHT));
	NetNode& superThree = results.back();

	results.emplace_back("four", Edge<NetNode>(8, F, 1.0), Edge<NetNode>(9, G, 1.0));
	NetNode& four = results.back();

	results.emplace_back("five", Edge<NetNode>(10, superThree, 1.0), Edge<NetNode>(11, four, 1.0));
	NetNode& five = results.back();

	results.emplace_back("six", Edge<NetNode>(12, two, 1.0), Edge<NetNode>(13, five, 1.0));
	NetNode& six = results.back();

	return six;
}


/**
 * Create a simple species network with a sophisticated introgression.
 */
inline NetNode& createSpeciesWithIntro(std::vector<NetNode>& results, double dx = 0) {
	results.reserve(15);

	results.emplace_back("A");
	NetNode& A = results.back();

	results.emplace_back("B");
	NetNode& B = results.back();

	results.emplace_back("C");
	NetNode& C = results.back();

	results.emplace_back("D");
	NetNode& D = results.back();

	results.emplace_back("E");
	NetNode& E = results.back();

	results.emplace_back("F");
	NetNode& F = results.back();

	results.emplace_back("G");
	NetNode& G = results.back();

	results.emplace_back("one", Edge<NetNode>(0, B, 1.0), Edge<NetNode>(1, C, 1.0));
	NetNode& one = results.back();

	results.emplace_back("oneIntrogressed", Edge<NetNode>(2, one, 0.0), 0.25 + dx, 15);
	NetNode& oneIntrogressed = results.back();

	results.emplace_back("two", Edge<NetNode>(3, A, 1.0), Edge<NetNode>(4, oneIntrogressed, 0.0, EdgeType::LEFT));
	NetNode& two = results.back();

	results.emplace_back("three", Edge<NetNode>(5, D, 1.0), Edge<NetNode>(6, E, 1.0));
	NetNode& three = results.back();

	results.emplace_back("superThree", Edge<NetNode>(7, three, 0.0), Edge<NetNode>(8, oneIntrogressed, 0.0, EdgeType::RIGHT));
	NetNode& superThree = results.back();

	results.emplace_back("four", Edge<NetNode>(9, F, 1.0), Edge<NetNode>(10, G, 1.0));
	NetNode& four = results.back();

	results.emplace_back("five", Edge<NetNode>(11, superThree, 1.0), Edge<NetNode>(12, four, 1.0));
	NetNode& five = results.back();

	results.emplace_back("six", Edge<NetNode>(13, two, 1.0), Edge<NetNode>(14, five, 1.0));
	NetNode& six = results.back();

	return six;
}


/**
 * Create a simple gene tree.
 */
inline TreeNode& createGene(std::vector<TreeNode>& results) {
	results.reserve(13);

	results.emplace_back("A");
	TreeNode& A = results.back();

	results.emplace_back("B");
	TreeNode& B = results.back();

	results.emplace_back("C");
	TreeNode& C = results.back();

	results.emplace_back("D");
	TreeNode& D = results.back();

	results.emplace_back("E");
	TreeNode& E = results.back();

	results.emplace_back("F");
	TreeNode& F = results.back();

	results.emplace_back("G");
	TreeNode& G = results.back();

	results.emplace_back("one", A, B);
	TreeNode& one = results.back();

	results.emplace_back("two", D, F);
	TreeNode& two = results.back();

	results.emplace_back("three", E, G);
	TreeNode& three = results.back();

	results.emplace_back("four", two, three);
	TreeNode& four = results.back();

	results.emplace_back("five", C, four);
	TreeNode& five = results.back();

	results.emplace_back("six", one, five);
	TreeNode& six = results.back();

	return six;
}


/**
 * Create a simple species network.
 */
inline NetNode& createSimpleSpecies(std::vector<NetNode>& results, double* params) {
	results.reserve(7);

	results.emplace_back("A");
	NetNode& A = results.back();

	results.emplace_back("B");
	NetNode& B = results.back();

	results.emplace_back("C");
	NetNode& C = results.back();

	results.emplace_back("BIntrogressed", Edge<NetNode>(0, B, params[0]), params[7], 7);
	NetNode& BIntrogressed = results.back();

	results.emplace_back("one", Edge<NetNode>(1, A, params[1]), Edge<NetNode>(2, BIntrogressed, params[2], EdgeType::LEFT));
	NetNode& one = results.back();

	results.emplace_back("two", Edge<NetNode>(3, BIntrogressed, params[3], EdgeType::RIGHT), Edge<NetNode>(4, C, params[4]));
	NetNode& two = results.back();

	results.emplace_back("three", Edge<NetNode>(5, one, params[5]), Edge<NetNode>(6, two, params[6]));
	NetNode& three = results.back();

	return three;
}


/**
 * Create an even simpler gene tree.
 */
inline TreeNode& createSimpleGene(std::vector<TreeNode>& results) {
	results.reserve(5);

	results.emplace_back("A");
	TreeNode& A = results.back();

	results.emplace_back("B");
	TreeNode& B = results.back();

	results.emplace_back("C");
	TreeNode& C = results.back();

	results.emplace_back("one", A, B);
	TreeNode& one = results.back();

	results.emplace_back("two", one, C);
	TreeNode& two = results.back();

	return two;
}

/**
 * Create an even simpler gene tree.
 */
inline TreeNode& createSimpleGeneTwo(std::vector<TreeNode>& results) {
	results.reserve(5);

	results.emplace_back("A");
	TreeNode& A = results.back();

	results.emplace_back("B");
	TreeNode& B = results.back();

	results.emplace_back("C");
	TreeNode& C = results.back();

	results.emplace_back("one", A, C);
	TreeNode& one = results.back();

	results.emplace_back("two", one, B);
	TreeNode& two = results.back();

	return two;
}

/**
 * Create an even simpler gene tree.
 */
inline TreeNode& createSimpleGeneThree(std::vector<TreeNode>& results) {
	results.reserve(5);

	results.emplace_back("A");
	TreeNode& A = results.back();

	results.emplace_back("B");
	TreeNode& B = results.back();

	results.emplace_back("C");
	TreeNode& C = results.back();

	results.emplace_back("one", B, C);
	TreeNode& one = results.back();

	results.emplace_back("two", one, A);
	TreeNode& two = results.back();

	return two;
}