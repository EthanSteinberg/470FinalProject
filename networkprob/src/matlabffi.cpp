#include "matlabffi.h"

#include "densemap.h"
#include "netnode.h"
#include "example.h"

struct NetworkBuffer {
    std::vector<NetNode> data;
};

struct NetworkBuffer* allocNetworkBuffer(int size) {
    NetworkBuffer* result = new NetworkBuffer();
    result->data.reserve(size);
    return result;
}

void freeNetworkBuffer(struct NetworkBuffer* buffer) {
    delete buffer;
}

int createLeafNetNode(struct NetworkBuffer* buffer, const char* name) {
    buffer->data.emplace_back(name);
    return buffer->data.size() - 1;
}

Edge<NetNode> createEdgeFromStruct(struct NetworkBuffer* buffer, NetworkEdge edge) {
    return Edge<NetNode>(edge.index, buffer->data[edge.sourceNode], edge.distance, (EdgeType)edge.type);
}

int createTreeNetNode(struct NetworkBuffer* buffer, const char* name, struct NetworkEdge leftEdge, struct NetworkEdge rightEdge) {
    buffer->data.emplace_back(name, createEdgeFromStruct(buffer, leftEdge), createEdgeFromStruct(buffer, rightEdge));
    return buffer->data.size() - 1;
}

int createNetworkNetNode(struct NetworkBuffer* buffer, const char* name,  struct NetworkEdge bottomEdge, double leftProbability, int introgressionId) {
    buffer->data.emplace_back(name, createEdgeFromStruct(buffer, bottomEdge), leftProbability, introgressionId);
    return buffer->data.size() - 1;
}

struct TreeBuffer {
    std::vector<TreeNode> data;
};

struct TreeBuffer* allocTreeBuffer(int size) {
    TreeBuffer* result = new TreeBuffer();
    result->data.reserve(size);
    return result;
}

void freeTreeBuffer(struct TreeBuffer* buffer) {
    delete buffer;
}

int createLeafTreeNode(struct TreeBuffer* buffer, const char* name) {
    buffer->data.emplace_back(name);
    return buffer->data.size() - 1;
}

int createTreeTreeNode(struct TreeBuffer* buffer, const char* name, int leftDescendant, int rightDescendant) {
    buffer->data.emplace_back(name, buffer->data[leftDescendant], buffer->data[rightDescendant]);
    return buffer->data.size() - 1;
}

void changeParams(Network net, double* params) {
    net.buffer->data[net.rootNode].setParams(params);
}

double computeProbability(Network net, Tree tree, double* derivatives) {
    if (derivatives == nullptr) {
        return calcProbability(net.buffer->data[net.rootNode], tree.buffer->data[tree.rootNode], nullptr);
    } else {
        std::vector<double> derivativeResults;
        double prob = calcProbability(net.buffer->data[net.rootNode], tree.buffer->data[tree.rootNode], &derivativeResults);

        for (unsigned int i = 0; i < derivativeResults.size(); i++) {
            derivatives[i] = derivativeResults[i];
        }

        return prob;
    }

}