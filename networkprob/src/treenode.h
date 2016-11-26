#pragma once

#include <string>
#include <memory>
#include <vector>
#include <map>
#include <experimental/optional>
#include <iostream>

template<class T>
using optional = std::experimental::optional<T>;

/**
 * A node in a tree.
 */
class TreeNode {
public:
	/**
	 * A leaf node.
	 */
	TreeNode(std::string name) {
		isLeaf = true;
		this->name = name;
	}

	/**
	 * A tree node.
	 */
	TreeNode(std::string name, const TreeNode& a_leftChild, const TreeNode& a_rightChild): leftChild(&a_leftChild), rightChild(&a_rightChild) {
		isLeaf = false;
		this->name = name;
	}

	/**
	 * Print the tree.
	 */
	void print() const {
		if (isLeaf) {
			std::cout<<name;
		} else {
			std::cout<<'(';
			leftChild->print();
			std::cout<<',';
			rightChild->print();
			std::cout<<')'<<name;
		}

	}

	bool isLeaf;
	std::string name;

	const TreeNode* leftChild;
	const TreeNode* rightChild;

};

/**
 * Get all the names of all the nodes.
 */
inline void processTaxa(const TreeNode& node, std::vector<std::string>& result, bool leafOnly) {
	if (node.isLeaf == leafOnly) {
		result.push_back(node.name);
	}

	if (!node.isLeaf) {
		processTaxa(*(node.leftChild), result, leafOnly);
		processTaxa(*(node.rightChild), result, leafOnly);
	}
}

/**
 * Get a mapping for the taxa of the tree.
 */
inline std::map<std::string, int> getTaxa(const TreeNode& gene) {
	std::vector<std::string> temp;
	processTaxa(gene, temp, false);
	temp.resize(6, "invalid event");
	processTaxa(gene, temp, true);


	std::map<std::string, int> result;
	for (unsigned int i = 0; i < temp.size(); i++) {
		result[temp[i]] = i;
	}
	return result;
}

/**
 * Compute all the events in the tree.
 */
inline void processEvents(const TreeNode& node, const std::map<std::string, int>& taxa, std::vector<int>& result) {
	if (!node.isLeaf) {
		int left = taxa.find(node.leftChild->name)->second;
		int right = taxa.find(node.rightChild->name)->second;

		int current = 0;
		current |= (1 << left);
		current |= (1 << right);
		result.push_back(current);

		processEvents(*node.leftChild, taxa, result);
		processEvents(*node.rightChild, taxa, result);
	}
}

/**
 * Get all the events in the tree.
 */
inline std::vector<int> getEvents(const TreeNode& gene, const std::map<std::string, int>& taxa){
	std::vector<int> current;
	processEvents(gene, taxa, current);
	return current;
}