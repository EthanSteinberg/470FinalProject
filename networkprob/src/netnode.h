#pragma once

#include <string>
#include <memory>
#include <map>
#include <bitset>
#include <experimental/optional>
#include <vector>
#include <limits>

#include "densemap.h"
#include "mathutils.h"
#include "treenode.h"

template<class T>
using optional = std::experimental::optional<T>;

const bool debug = false;

enum class NodeType {
	LEAF = 0,
	TREE = 1,
	NETWORK = 2,
};

enum class EdgeType {
	NORMAL = 0,
	LEFT = 1,
	RIGHT = 2,
};

/**
 * An edge in a network.
 */
template<typename Node>
class Edge {
public:
	/**
	 * The constructor for an edge.
	 */
	Edge(unsigned int a_id, Node& a_toNode, double a_distance, EdgeType a_type = EdgeType::NORMAL) : id(a_id), toNode(a_toNode), distance(a_distance), type(a_type) {}

	/**
	 * Get the data associated with the edge.
	 */
	std::vector<densemap> getData(const std::map<std::string, int>& netNodes, const std::map<std::string, int>& taxa, const std::vector<int>& events, int numDerivativeParams) {
		auto result = update(toNode.getData(type, netNodes, taxa, events, numDerivativeParams), events, distance);

		if (debug) {
			std::cout<<"------------------------------------"<<std::endl;
			std::cout<<"From node: "<<toNode.name<<std::endl;

			for (auto&& map : result) {
				std::cout<<"next: "<<std::bitset<8>(map.getTaxaBits()>>6)<<' ';
				for (unsigned int i =0; i < map.choices.size(); i++) {
					std::cout<<map.choices[i]<<' ';
				}
				std::cout<<std::endl;
				for (int i = 0; i < 1<<6; i++){
					double thingy = map.getHistory(i);

					if (thingy != 0) {
						std::cout<<std::bitset<8>(i)<<' '<<thingy<<std::endl;
					}
				}
			}
		}

		return result;
	}

	/**
	 * Get the derivatives of the data for the edge.
	 */
	std::vector<std::vector<densemap>> getDataDerivative(const std::map<std::string, int>& netNodes, const std::map<std::string, int>& taxa, const std::vector<int>& events, int numDerivativeParams) {
		std::vector<std::vector<densemap>> result;
		auto derivative = toNode.getDataDerivative(type, netNodes, taxa, events, numDerivativeParams);

		for (unsigned int i = 0; i < derivative.size(); i++) {
			if (i == id) {
				// std::cout<<"Got one! at " <<i <<std::endl;
				// That means that I need to originate the derivative
				result.push_back(derivativeUpdate(toNode.getData(type, netNodes, taxa, events, numDerivativeParams), events, distance));
			} else {
				// This means that the derivative is hopefully farther down the line
				result.push_back(update(derivative[i], events, distance));
			}
		}

		if (debug) {
			std::cout<<"------------------------------------"<<std::endl;
			std::cout<<"From node: "<<toNode.name<<std::endl;

			for (unsigned int i = 0; i < result.size() ;i++) {
				for (auto&& map : result[i]) {
					std::cout<<"next: "<<std::bitset<8>(map.getTaxaBits()>>6)<<' ';
					for (unsigned int i =0; i < map.choices.size(); i++) {
						std::cout<<map.choices[i]<<' ';
					}
					std::cout<<std::endl;
					for (int i = 0; i < 1<<6; i++){
						double thingy = map.getHistory(i);

						if (thingy != 0) {
							std::cout<<std::bitset<8>(i)<<' '<<thingy<<std::endl;
						}
					}
				}
			}
		}

		return result;
	}

	/**
	 * Print the edge.
	 */
	void print() {
		std::cout<<toNode.name<<' '<<distance<<' '<<(int)type<<std::endl;
	}

	/**
	 * Print the children of hte edge.
	 */
	void printChildren() {
		toNode.print();
	}

	/**
	 * Get the maximum parameter index.
	 */
	int getMaximumParamId() const {
		return std::max((int)id, toNode.getMaximumParamId());
	}

	/**
	 * Set the parameters.
	 */
	void setParams(double* params) {
		distance = params[id];
		toNode.setParams(params);
	}

	unsigned int id; // The index for the edge.
	Node& toNode; // The node it points to.
	double distance; // The length of the edge.
	EdgeType type; // THe type of edge.
};

/**
 * A network node.
 */
class NetNode {
public:
	/**
	 * A leaf node.
	 */
	NetNode(std::string name) {
		type = NodeType::LEAF;
		this->name = name;
		this->initialized = false;
	}

	/**
	 * A tree node.
	 */
	NetNode(std::string name, Edge<NetNode> a_leftEdge, Edge<NetNode> a_rightEdge): leftEdge(a_leftEdge), rightEdge(a_rightEdge) {
		type = NodeType::TREE;
		this->name = name;
		this->initialized = false;
	}

	/**
	 * A network node.
	 */
	NetNode(std::string name, Edge<NetNode> a_childEdge, double leftProbability, unsigned int introgressionId): childEdge(a_childEdge) {
		type = NodeType::NETWORK;
		this->name = name;

		this->leftProbability = leftProbability;
		this->initialized = false;
		this->introgressionId = introgressionId;
	}

	/**
	 * Set the parameters.
	 */
	void setParams(double* params) {
		initialized = false;

		switch (type) {
			case NodeType::LEAF:
				break;

			case NodeType::TREE:
				leftEdge->setParams(params);
				rightEdge->setParams(params);
				break;

			case NodeType::NETWORK:
				leftProbability = params[introgressionId];
				childEdge->setParams(params);
				break;

			default:
				std::cerr<<"Unknown type"<<std::endl;
				exit(-1);
		}
	}

	/**
	 * Get the maximum parameter index.
	 */
	int getMaximumParamId() const {
		switch (type) {
			case NodeType::LEAF:
				return -1;
			case NodeType::TREE:
				return std::max(leftEdge->getMaximumParamId(), rightEdge->getMaximumParamId());
			case NodeType::NETWORK:
				return std::max(childEdge->getMaximumParamId(), (int)introgressionId);

			default:
				std::cerr<<"Unknown type"<<std::endl;
				exit(-1);
		}


	}

	/**
	 * Get the data for this node.
	 */
	const std::vector<densemap>& getData(EdgeType type, const std::map<std::string, int>& netNodes, const std::map<std::string, int>& taxa, const std::vector<int>& events, int numDerivativeParams) {
		if (!initialized) {
			initialized = true;
			computeDenseMap(netNodes, taxa, events, numDerivativeParams);
		}

		switch (type) {
			case EdgeType::NORMAL:
				return currentData;
			case EdgeType::LEFT:
				return leftData;
			case EdgeType::RIGHT:
				return rightData;

			default:
				std::cerr<<"Unknown type"<<std::endl;
				exit(-1);
		}
	}

	/**
	 * Get the derivative of the data for this node.
	 */
	const std::vector<std::vector<densemap>> & getDataDerivative(EdgeType type, const std::map<std::string, int>& netNodes, const std::map<std::string, int>& taxa, const std::vector<int>& events, int numDerivativeParams) {
		if (!initialized) {
			initialized = true;
			computeDenseMap(netNodes, taxa, events, numDerivativeParams);
		}

		switch (type) {
			case EdgeType::NORMAL:
				return derivatives;
			case EdgeType::LEFT:
				return leftDerivatives;
			case EdgeType::RIGHT:
				return rightDerivatives;

			default:
				std::cerr<<"Unknown type"<<std::endl;
				exit(-1);
		}
	}

	/**
	 * Print the node and children.
	 */
	void print() {
		if (type == NodeType::LEAF) {
			std::cout<<"Leaf: "<<name<<std::endl;
		} else if (type == NodeType::TREE) {
			std::cout<<"Tree: "<<name<<std::endl;
			leftEdge->print();
			rightEdge->print();
			leftEdge->printChildren();
			rightEdge->printChildren();
		} else if (type == NodeType::NETWORK) {
			std::cout<<"Net: "<<name<<' '<<leftProbability<<std::endl;
			childEdge->print();
			childEdge->printChildren();
		}
	}

	/**
	 * Compute the values for this node.
	 */
	void computeDenseMap(const std::map<std::string, int>& netNodes, const std::map<std::string, int>& taxa, const std::vector<int>& events, int numDerivativeParams) {
		if (type == NodeType::LEAF) {

			int id = taxa.find(name)->second;
			currentData.resize(1);
			std::vector<int64_t> choices(netNodes.size(), -1);
			currentData[0].init(1 << id, choices);
			currentData[0].setHistory(0, 1.0);

			derivatives.resize(numDerivativeParams);
			for (int i = 0 ; i < numDerivativeParams;i++) {
				std::vector<densemap> nextMap;
				nextMap.resize(1);
				nextMap[0].init(1 << id, choices);
				derivatives[i] = nextMap;
			}
		} else if (type == NodeType::TREE) {
			std::vector<densemap> left = leftEdge->getData(netNodes, taxa, events, numDerivativeParams);
			std::vector<densemap> right = rightEdge->getData(netNodes, taxa, events, numDerivativeParams);

			currentData = combine(left, right);

			auto leftDerivatives = leftEdge->getDataDerivative(netNodes, taxa, events, numDerivativeParams);
			auto rightDerivatives = rightEdge->getDataDerivative(netNodes, taxa, events, numDerivativeParams);

			derivatives.resize(leftDerivatives.size());
			for (unsigned int i = 0; i < leftDerivatives.size(); i++) {
				derivatives[i] = combineDerivatives(left, leftDerivatives[i], right, rightDerivatives[i]);
			}

			if (debug) {
				std::cout<<"---------------------------------"<<std::endl;
				std::cout<<"Computed node: "<<name<<std::endl;

				for (auto&& map : currentData) {
					std::cout<<"next: "<<std::bitset<8>(map.getTaxaBits()>>6)<<' ';
					for (unsigned int i =0; i < map.choices.size(); i++) {
						std::cout<<map.choices[i]<<' ';
					}
					std::cout<<std::endl;
					for (int i = 0; i < 1<<6; i++){
						double thingy = map.getHistory(i);

						if (thingy != 0) {
							std::cout<<std::bitset<8>(i)<<' '<<thingy<<std::endl;
						}
					}
				}
			}

		} else if (type == NodeType::NETWORK) {
			std::vector<densemap> child = childEdge->getData(netNodes, taxa, events, numDerivativeParams);

			std::vector<std::vector<densemap>> childDerivatives = childEdge->getDataDerivative(netNodes, taxa, events, numDerivativeParams);

			int netNodeId = netNodes.find(name)->second;

			std::tie(leftData, rightData) = split(child, netNodeId, events, leftProbability);

			leftDerivatives.resize(childDerivatives.size());
			rightDerivatives.resize(childDerivatives.size());

			for (unsigned int i = 0; i < childDerivatives.size() ; i++) {
				if (i == introgressionId) {
					std::tie(leftDerivatives[i], rightDerivatives[i]) = splitDerivativeHere(child, netNodeId, events, leftProbability);
				} else {
					std::tie(leftDerivatives[i], rightDerivatives[i]) = splitDerivatives(childDerivatives[i], child, netNodeId, events, leftProbability);
				}
			}

			if (debug) {
				std::cout<<"---------------------------------"<<std::endl;
				std::cout<<"Computed node: "<<name<<std::endl;

				std::cout<<"Left"<<std::endl;
				for (auto&& map : leftData) {
					std::cout<<"next: "<<std::bitset<8>(map.getTaxaBits()>>6)<<' ';
					for (unsigned int i =0; i < map.choices.size(); i++) {
						std::cout<<std::bitset<64>(map.choices[i])<<' ';
					}
					std::cout<<std::endl;
					for (int i = 0; i < 1<<6; i++){
						double thingy = map.getHistory(i);

						if (thingy != 0) {
							std::cout<<std::bitset<8>(i)<<' '<<thingy<<std::endl;
						}
					}
				}

				std::cout<<"Right"<<std::endl;
				for (auto&& map : rightData) {
					std::cout<<"next: "<<std::bitset<8>(map.getTaxaBits()>>6)<<' ';
					for (unsigned int i =0; i < map.choices.size(); i++) {
						std::cout<<map.choices[i]<<' ';
					}
					std::cout<<std::endl;
					for (int i = 0; i < 1<<6; i++){
						double thingy = map.getHistory(i);

						if (thingy != 0) {
							std::cout<<std::bitset<8>(i)<<' '<<thingy<<std::endl;
						}
					}
				}
			}
		}
	}

	NodeType type;
	std::string name;
	bool initialized;

	std::vector<densemap> currentData;
	std::vector<std::vector<densemap>> derivatives;

	std::vector<densemap> leftData;
	std::vector<densemap> rightData;

	std::vector<std::vector<densemap>> leftDerivatives;
	std::vector<std::vector<densemap>> rightDerivatives;

	// Tree node properties
	optional<Edge<NetNode>> leftEdge;
	optional<Edge<NetNode>> rightEdge;

	// Network node properties
	optional<Edge<NetNode>> childEdge;
	double leftProbability;

	unsigned int introgressionId;
};

/**
 * Get all the network nodes into result.
 */
inline void processNetNodes(const NetNode& species, std::vector<std::string>& result) {

	switch (species.type) {
		case NodeType::LEAF:
			return;
		case NodeType::TREE:
			processNetNodes(species.leftEdge->toNode, result);
			processNetNodes(species.rightEdge->toNode, result);
			return;
		case NodeType::NETWORK:
			if (std::find(result.begin(), result.end(), species.name) == result.end())
				result.push_back(species.name);
			processNetNodes(species.childEdge->toNode, result);
			return;
	}
}

/**
 * Get a mapping of network nodes to indices.
 */
inline std::map<std::string, int> getNetNodes(const NetNode& species) {
	std::vector<std::string> temp;
	processNetNodes(species, temp);

	std::map<std::string, int> result;
	for (unsigned int i = 0; i < temp.size(); i++) {
		result[temp[i]] = i;
	}
	return result;
}

/**
 * Compute the probability of a gene tree given a species tree.
 * Also computes the derivatives if it is not nullptr.
 */
inline double calcProbability(NetNode& species, TreeNode& geneTree, std::vector<double>* derivatives = nullptr) {
	auto taxa = getTaxa(geneTree);
	auto events = getEvents(geneTree, taxa);

	if (debug) {
		for (int event : events) {
			std::cout<<std::bitset<16>(event)<<',';
		}
		std::cout<<std::endl;
	}

	auto netNodes = getNetNodes(species);
	auto rootEdge = Edge<NetNode>(-1, species, std::numeric_limits<double>::infinity());

	int numDerivativeParams = rootEdge.getMaximumParamId() + 1;

	const std::vector<densemap>& root = rootEdge.getData(netNodes, taxa, events, (derivatives != nullptr) ?  numDerivativeParams : 0 );

	int maxTaxa = 0;
	for (const auto& entry: taxa) {
		maxTaxa = std::max(maxTaxa, entry.second);
	}

	uint16_t targetTaxaBits = 0;
	for (int i = 6; i <= maxTaxa; i++) {
		targetTaxaBits |= 1 << i;
	}

	double probability = 0.0;

	for(auto&& map: root) {
		if (map.getTaxaBits() == targetTaxaBits) {
			probability += map.getHistory(map.getMaxHistory(events.size()) - 1);
		}
	}

	if (derivatives != nullptr) {
		auto derivativeRoot = rootEdge.getDataDerivative(netNodes, taxa, events, numDerivativeParams);

		for (int i = 0; i < numDerivativeParams; i++) {
			double nextVal = 0.0;
			for(auto&& map: derivativeRoot[i]) {
				if (map.getTaxaBits() == targetTaxaBits) {
					nextVal += map.getHistory(map.getMaxHistory(events.size()) - 1);
				}
			}
			derivatives->push_back(nextVal);
		}
	}

	return probability;
}