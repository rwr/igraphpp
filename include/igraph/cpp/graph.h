/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_GRAPH_H
#define IGRAPHPP_GRAPH_H

#include <cassert>
#include <memory>
#include <igraph/igraph_constructors.h>
#include <igraph/igraph_datatype.h>
#include <igraph/igraph_foreign.h>
#include <igraph/igraph_interface.h>
#include <igraph/igraph_conversion.h>
#include <igraph/cpp/attributes.h>
#include <stdexcept>
#include <vector>

#include <igraph/cpp/error.h>
#include <igraph/cpp/types.h>
#include <igraph/cpp/vector.h>
#include <unordered_set>

namespace igraph {

class Edge;
class EdgeSelector;
class Vertex;
class VertexSelector;

/// C++-style wrapper around an igraph_t object
class Graph {
private:
    /// The igraph_t instance encapsulated by the wrapper
    igraph_t* m_pGraph;
    /// Set of fast deleted edges. Works only on graphs without multiedges.
    std::unordered_set<integer_t> deleted_edges;

public:
    /*****************************/
    /* Constructors, destructors */
    /*****************************/

    /// Constructs an empty graph
    Graph(long numVertices = 0, bool directed = false) : m_pGraph(new igraph_t) {
        IGRAPH_TRY(igraph_empty(m_pGraph, numVertices, directed));
    }

    /// Constructs a wrapper that wraps the given igraph_t instance
    /**
     * The ownership of the wrapped instance is stolen by the wrapper.
     * The caller should not destroy the graph on its own, ever;
     * the wrapper should be destroyed instead.
     *
     * This function never throws an exception.
     */
    Graph(igraph_t* graph) throw() : m_pGraph(graph) {
        assert(graph);
    }

    /// Constructs a wrapper that wraps the given igraph_t instance
    /**
     * The ownership of the wrapped instance is stolen by the wrapper.
     * The caller should not destroy the graph on its own, ever;
     * the wrapper should be destroyed instead.
     *
     * This function never throws an exception.
     */
    Graph(std::unique_ptr<igraph_t> graph) throw() : m_pGraph(graph.release()) {
    }

    /// Copy constructor
    Graph(const Graph& other) : m_pGraph(new igraph_t) {
        IGRAPH_TRY(igraph_copy(m_pGraph, other.m_pGraph));
    }

    /// Destroys the graph
    ~Graph();

    /******************/
    /* Static methods */
    /******************/

    /********************/
    /* Instance methods */
    /********************/

    /// Adds a single edge to the graph
    void addEdge(integer_t source, integer_t target);

    /// Adds a single edge to the graph (the fast version, works only for simple graphs)
    void addEdge_fast(igraph_integer_t source, igraph_integer_t target);

    /// Adds a list of edges to the graph
    void addEdges(const Vector& edges);

    /// Adds a list of edges to the graph (the fast version, works only for simple graphs)
    void addEdges_fast(const Vector& edges);

    /// Adds a single vertex to the graph
    void addVertex();

    /// Adds the given number of vertices to the graph
    void addVertices(long numVertices);

    /// Returns whether the two vertices are connected
    bool areConnected(long u, long v) const;

    /// Returns whether the two vertices are connected (the fast version, works only for simple graphs)
    bool areConnected_fast(long int u, long int v) const;

    /// Extracts a pointer to the encapsulated graph object
    igraph_t* c_graph() { return m_pGraph; }

    /// Extracts a pointer to the encapsulated graph object (const)
    const igraph_t* c_graph() const { return m_pGraph; }

    /// Returns the degrees of some vertices
    Vector degree(const VertexSelector& vids,
            NeighborMode mode = IGRAPH_ALL, bool loops = false) const;
    /// Returns the degrees of some vertices
    void degree(Vector* result, const VertexSelector& vids,
            NeighborMode mode = IGRAPH_ALL, bool loops = false) const;

    /// Returns the degree of the given vertex if the graph is directed and the degree otherwise.
    /// The fast version, works only for simple graphs.
    integer_t degree_fast(long int v, igraph_neimode_t mode = IGRAPH_ALL){
        Vector eids;
        igraph_incident(m_pGraph, eids.c_vector(), v, mode);
        size_t count = 0;
        for (auto eid : eids)
            if (! deleted_edges.count(eid))
                count++;
        return count;
    }

    /// Returns the in-degree of the given vertex. If the graph is directed, this is the sum of in-degrees and out-degrees,
    /// otherwise, just the degree of the vertex.
    /// The fast version, works only for simple graphs.
    integer_t in_degree_fast(long int v){
        return degree_fast(v, IGRAPH_IN);
    }

    /// Returns the out-degree of the given vertex. If the graph is directed, this is the sum of in-degrees and out-degrees,
    /// otherwise, just the degree of the vertex.
    /// The fast version, works only for simple graphs.
    integer_t out_degree_fast(long int v){
        return degree_fast(v, IGRAPH_OUT);
    }

    /// Returns the degree of the given vertex. If the graph is directed, this is the sum of in-degrees and out-degrees,
    /// otherwise, just the degree of the vertex.
    /// The fast version, works only for simple graphs.
    integer_t degree_fast(long int v){
        return degree_fast(v, IGRAPH_ALL);
    }

    /// Deletes some edges from the graph
    void deleteEdges(const EdgeSelector& es);

    /// Deletes one edge given by the endpoints.
    void deleteEdge(long u, long v);

    /// Deletes one edge given by the endpoints. If the edge does not exist, does nothing.
    /// The fast version, works only for simple graphs.
    void deleteEdgeFast(long u, long v);

    /// Returns the head and tail vertices of an edge
    void edge(integer_t eid, integer_t* from, integer_t* to) const;

    /// Returns the edge with the given index
    Edge edge(integer_t eid);

    /// Returns the number of edges in the graph
    integer_t ecount() const { return igraph_ecount(m_pGraph); }

    /// Returns the number of edges in the graph
    /// The fast version, works only for simple graphs.
    integer_t ecount_fast() const { return igraph_ecount(m_pGraph) - deleted_edges.size(); }


    /// Returns a copy of the value of the given graph attribute
    AttributeValue getAttribute(const std::string& attribute) const;

    /// Returns the edge list of the graph
    Vector getEdgelist(bool bycol=false) const;

    /// Returns the edge list of the graph.
    /// The fast version, works only for simple graphs.
    Vector getEdgelist_fast(bool bycol=false) const;

    /// Returns the edge list of the graph
    void getEdgelist(Vector* result, bool bycol=false) const;

    /// Returns the ID of an arbitrary edge between the two given nodes
    integer_t getEid(integer_t source, integer_t target, bool directed=true,
            bool error=false) const;

    /// Tests whether an edge given by source and target is fast deleted.
    bool isDeleted(integer_t source, integer_t target) const{
        integer_t eid = getEid(source, target);
        return deleted_edges.count(eid);
    }

    /// Tests whether an edge given by edge id is fast deleted.
    bool isDeleted(integer_t eid) const{
        return deleted_edges.count(eid);
    }

    /// Returns whether the graph has the given graph attribute
    bool hasAttribute(const std::string& attribute) const;

    /// Returns the edges incident on a given vertex
    void incident(Vector* result, long int vertex, NeighborMode mode = IGRAPH_OUT) const;
    /// Returns the edges incident on a given vertex
    Vector incident(long int vertex, NeighborMode mode = IGRAPH_OUT) const;

    /// Returns the edges incident on a given vertex
    /// The fast version, works only for simple graphs.
    Vector incident_fast(long int vertex, NeighborMode mode = IGRAPH_OUT) const{
        Vector eids_prelim, eids;
        igraph_incident(m_pGraph, eids_prelim.c_vector(), vertex, mode);
        for (size_t i=0; i < eids_prelim.size(); i++)
            if (! isDeleted(eids_prelim[i]))
                eids.push_back(eids_prelim[i]);
        return eids;
    };

    /// Returns whether the graph is directed
    bool isDirected() const;

    /// Checks whether the graph is a simple graph
    bool isSimple() const;

    /// Checks whether the graph is a simple graph accounting edges deleted fast
    int isSimple_fast(const igraph_t *graph, igraph_bool_t *res);

    /// Returns the neighbors of a vertex
    void neighbors(Vector* result, long int vertex, NeighborMode mode = IGRAPH_OUT) const;

    /// Returns the neighbors of a vertex
    Vector neighbors(long int vertex, NeighborMode mode = IGRAPH_OUT) const;

    /// Returns the in_neighbors of a vertex
    /// The fast version, works only for simple graphs.
    Vector in_neighbors_fast(long int vertex) const {
        if (!isDirected())
            return neighbors_fast(vertex);

        Vector eids = incident_fast(vertex, IGRAPH_IN);
        Vector neighbors;
        for (size_t i = 0; i < eids.size(); i++)
            if (!isDeleted(eids[i]))
                neighbors.push_back(IGRAPH_FROM(m_pGraph, eids[i]));
        return neighbors;
    }



    /// Returns the neighbors of a vertex. If the graph is directed, in-neighbors and out-neighbors are returned.
    /// The fast version, works only for simple graphs.
    Vector neighbors_fast(integer_t vertex, NeighborMode mode = IGRAPH_ALL) const{
        Vector neighbors;
        Vector eids;
        igraph_incident(m_pGraph, eids.c_vector(), vertex, mode);

        for (auto eid : eids)
            if (! deleted_edges.count(eid)){
                integer_t from;
                edge(eid, &from, &vertex);
                neighbors.push_back(from);
            }
        return neighbors;
    }

    /// Returns the in-neighbors of a vertex in a directed graph or all neighbors in an undirected graph
    Vector in_neighbors(integer_t vertex) const{
        Vector eids;
        igraph_incident(m_pGraph, eids.c_vector(), vertex, IGRAPH_IN);
        return eids;
    };

    /// Returns the in-neighbors of a vertex in a directed graph or all neighbors in an undirected graph.
    /// The fast version, works only for simple graphs.
    Vector in_neighbors_fast(integer_t vertex) const{
        Vector eids_prelim, eids;
        igraph_incident(m_pGraph, eids_prelim.c_vector(), vertex, IGRAPH_IN);
        for (size_t i=0; i<eids_prelim.size(); i++)
            if (! isDeleted(eids_prelim[i]))
                eids.push_back(eids_prelim[i]);
        return eids;
    };

    /// Returns the out-neighbors of a vertex in a directed graph or all neighbors in an undirected graph
    Vector out_neighbors(long int vertex) const{
        Vector eids;
        igraph_incident(m_pGraph, eids.c_vector(), vertex, IGRAPH_OUT);
        return eids;
    };

    /// Returns the in-neighbors of a vertex in a directed graph or all neighbors in an undirected graph.
    /// The fast version, works only for simple graphs.
    Vector out_neighbors_fast(integer_t vertex) const{
        Vector eids_prelim, eids;
        igraph_incident(m_pGraph, eids_prelim.c_vector(), vertex, IGRAPH_OUT);
        for (size_t i=0; i<eids_prelim.size(); i++)
            if (! isDeleted(eids_prelim[i]))
                eids.push_back(eids_prelim[i]);
        return eids;
    };


    /// Sets the value of the given graph attribute
    void setAttribute(const std::string& attribute, const AttributeValue& value);

    /// Removes loop and/or multiple edges from the graph
    void simplify(bool multiple=true, bool loops=true);

    /// Returns the number of vertices in the graph
    integer_t vcount() const { return igraph_vcount(m_pGraph); }

    /// Returns the vertex with the given index in the graph
    Vertex vertex(integer_t vid);

    Graph induced_subgraph(Vector const& allowed_vertices) const;

    /// Makes the undirected graph bidirected.
    /// If the graph is already directed, does nothing.
    void make_bidirected(){
        if (isDirected())
            return;
        IGRAPH_TRY(igraph_to_directed(m_pGraph, IGRAPH_TO_DIRECTED_MUTUAL));
    }
    /*************/
    /* Operators */
    /*************/

    /// Disjoint union of two graphs
    Graph operator+(const Graph& other) const;

    /// Assignment operator
    Graph& operator=(const Graph& other) {
        if (&other == this)
            return *this;

        igraph_t new_graph;
        IGRAPH_TRY(igraph_copy(&new_graph, other.m_pGraph));

        if (m_pGraph) {
            igraph_destroy(m_pGraph);
        } else {
            m_pGraph = new igraph_t;
        }
        *m_pGraph = new_graph;

        deleted_edges = other.deleted_edges; // Changed for fast deleting edges.

        return *this;
    }

    /// Retrieves the value of the given graph attribute
    /**
     * This method works similar to \c std::map<>.operator[]: if the attribute
     * is found, its value is returned; if the attribute is not found, a new
     * attribute will be created with the default constructor of \c AttributeValue
     * and this will be returned. Therefore, this operator won't work on const
     * graphs.
     */
    AttributeValue& operator[](const std::string& attribute);

    /// Prints the edge list to cout.
    void print_edges() const;

    /// Prints the edge list to cout.
    /// The fast version.
    void print_edges_fast() const{
        Vector edges = getEdgelist_fast();
        for (auto source_idx = 0; source_idx < edges.size(); source_idx += 2)
            std::cout << edges[source_idx] << "," << edges[source_idx+1] << "  ";
        std::cout << std::endl;
    };

private:
    /// Returns a pointer to the attribute holder of the graph
    AttributeHolder* getAttributeHolder();

    /// Returns a pointer to the attribute holder of the graph (const)
    const AttributeHolder* getAttributeHolder() const;

    friend class Edge;
    friend class EdgeSelector;
    friend class Vertex;
    friend class VertexSelector;
};

}       // end of namespaces

#endif  // IGRAPHPP_GRAPH_H
