#ifndef MATLAB_FFI_INCLUDED
#define MATLAB_FFI_INCLUDED

#ifdef __cplusplus
extern "C" {
#endif

    /**
     * A network buffer stores network nodes.
     */
    struct NetworkBuffer;
    struct NetworkBuffer* allocNetworkBuffer(int size);
    void freeNetworkBuffer(struct NetworkBuffer* buffer);

    enum NetworkEdgeType {
        NORMAL = 0,
        LEFT = 1,
        RIGHT = 2,
    };

    struct NetworkEdge {
        int index; // The index in the params vector.
        int sourceNode;
        double distance;
        enum NetworkEdgeType type;
    };

    /**
     * Functions for creating network nodes.
     */
    int createLeafNetNode(struct NetworkBuffer* buffer, const char* name);
    int createTreeNetNode(struct NetworkBuffer* buffer, const char* name, struct NetworkEdge leftEdge, struct NetworkEdge rightEdge);
    int createNetworkNetNode(struct NetworkBuffer* buffer, const char* name,  struct NetworkEdge bottomEdge, double leftProbability, int introgressionId);

    /**
     * A tree buffer stores tree nodes.
     */
    struct TreeBuffer;
    struct TreeBuffer* allocTreeBuffer(int size);
    void freeTreeBuffer(struct TreeBuffer* buffer);

    /**
     * Functions for creating tree nodes.
     */
    int createLeafTreeNode(struct TreeBuffer* buffer, const char* name);
    int createTreeTreeNode(struct TreeBuffer* buffer, const char* name, int leftDescendant, int rightDescendant);

    struct Network {
        struct NetworkBuffer* buffer;
        int rootNode;
    };

    struct Tree {
        struct TreeBuffer* buffer;
        int rootNode;
    };

    /**
     * Change all the params for the given tree.
     */
    void changeParams(struct Network net, double* params);

    /**
     * Compute the probability of a gene tree given a network.
     * If derivatives is non-null, then also computes the derivatives.
     */
    double computeProbability(struct Network net, struct Tree tree, double* derivatives);

#ifdef __cplusplus
}
#endif

#endif