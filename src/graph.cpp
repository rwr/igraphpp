/* vim:set ts=4 sw=4 sts=4 et: */

#include <igraph/igraph_conversion.h>
#include <igraph/igraph_foreign.h>
#include <igraph/igraph_games.h>
#include <igraph/igraph_operators.h>
#include <igraph/igraph_structural.h>

#include <igraph/cpp/edge.h>
#include <igraph/cpp/edge_selector.h>
#include <igraph/cpp/graph.h>
#include <igraph/cpp/vertex.h>
#include <igraph/cpp/vertex_selector.h>
#include <igraph/cpp/util/allocation_guard.h>
#include <memory>

namespace igraph {

/// Destructor
Graph::~Graph() {
    if (m_pGraph) {
        igraph_destroy(m_pGraph);
        m_pGraph = nullptr;
    }
}

/******************/
/* Static methods */
/******************/

/********************/
/* Instance methods */
/********************/

void Graph::addEdge(igraph_integer_t source, igraph_integer_t target) {
    Vector edge(2);
    edge[0] = source; edge[1] = target;
    addEdges(edge);
}

void Graph::addEdges(const Vector& edges) {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_add_edges(m_pGraph, edges.c_vector(), 0));
}

void Graph::addEdge_fast(igraph_integer_t source, igraph_integer_t target) {
    if (isDeleted(source, target))
        deleted_edges.erase(getEid(source, target));
    Vector edge(2);
    edge[0] = source; edge[1] = target;
    addEdges(edge);
}

void Graph::addEdges_fast(const Vector& edges) {
    assert(m_pGraph);
    Vector new_edges;
    for (int i = 0; i < edges.size(); i += 2) {
        integer_t eid = getEid(edges[i], edges[i + 1]);
        if (isDeleted(eid))
            deleted_edges.erase(eid);
        else {
            new_edges.push_back(edges[i]);
            new_edges.push_back(edges[i+1]);
        }
    }
    IGRAPH_TRY(igraph_add_edges(m_pGraph, new_edges.c_vector(), 0));
}

void Graph::addVertex() { addVertices(1); }

void Graph::addVertices(long int numVertices) {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_add_vertices(m_pGraph, numVertices, 0));
}

bool Graph::areConnected(long int u, long int v) const {
    igraph_bool_t result;
    assert(m_pGraph);
    IGRAPH_TRY(igraph_are_connected(m_pGraph, u, v, &result));

    return result;
}

bool Graph::areConnected_fast(long int u, long int v) const {
    igraph_bool_t result = false;
    assert(m_pGraph);
    if (isDeleted(u, v))
        return result;
    IGRAPH_TRY(igraph_are_connected(m_pGraph, u, v, &result));
    return result;
}


Vector Graph::degree(const VertexSelector& vids, NeighborMode mode, bool loops) const {
    Vector result;
    degree(&result, vids, mode, loops);
    return result;
}

void Graph::degree(Vector* result, const VertexSelector& vids,
                   NeighborMode mode, bool loops) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_degree(m_pGraph, result->c_vector(), *vids.c_vs(),
                mode, loops));
}

void Graph::deleteEdges(const EdgeSelector& es) {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_delete_edges(m_pGraph, *es.c_es()));
}

void Graph::deleteEdge(long u, long v){
    assert(m_pGraph);
    igraph_es_t c_es;
    igraph_es_pairs_small(&c_es, isDirected(), u, v, -1);
    IGRAPH_TRY(igraph_delete_edges(m_pGraph, c_es));
}

void Graph::deleteEdgeFast(long u, long v){
    if (areConnected_fast(u, v))
        deleted_edges.insert(getEid(u, v));
}

void Graph::edge(integer_t eid, integer_t* from, integer_t* to) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_edge(m_pGraph, eid, from, to));
}

Edge Graph::edge(integer_t eid) {
    assert(m_pGraph);
    return Edge(this, eid);
}

AttributeValue Graph::getAttribute(const std::string& attribute) const {
    return getAttributeHolder()->getGraphAttribute(attribute);
}

void Graph::getEdgelist(Vector* result, bool bycol) const {
    assert(m_pGraph);
    Vector prelim_result;
    IGRAPH_TRY(igraph_get_edgelist(m_pGraph, prelim_result.c_vector(), bycol));
}

Vector Graph::getEdgelist_fast(bool bycol) const {
    Vector result_prelim, result;
    getEdgelist(&result_prelim, bycol);
    for (int i = 0; i < result_prelim.size(); i++)
        if (! isDeleted(result_prelim[i], result_prelim[i+1])) {
            result.push_back(result_prelim[i]);
            result.push_back(result_prelim[i+1]);
        }
    return result;
}

integer_t Graph::getEid(integer_t source, integer_t target,
        bool directed, bool error) const {
    integer_t eid;
    igraph_get_eid(m_pGraph, &eid, source, target, isDirected(), true);
    if (deleted_edges.count(eid))
        return -1;
    IGRAPH_TRY(igraph_get_eid(m_pGraph, &eid, source, target, directed, error));
    return eid;
}

bool Graph::isDirected() const {
    return igraph_is_directed(m_pGraph);
}

bool Graph::isSimple() const {// TODO account edges deleted fast
    bool_t result;
    IGRAPH_TRY(igraph_is_simple(m_pGraph, &result));
    return result;
}

bool Graph::hasAttribute(const std::string& attribute) const {
    return getAttributeHolder()->hasGraphAttribute(attribute);
}

void Graph::incident(Vector* result, long int vertex, NeighborMode mode) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_incident(m_pGraph, result->c_vector(), vertex, mode));
}

Vector Graph::incident(long int vertex, NeighborMode mode) const {
    Vector result;
    incident(&result, vertex, mode);
    return result;
}

void Graph::neighbors(Vector* result, long int vertex, NeighborMode mode) const {
    assert(m_pGraph);
    IGRAPH_TRY(igraph_neighbors(m_pGraph, result->c_vector(), vertex, mode));
}

Vector Graph::neighbors(long int vertex, NeighborMode mode) const {
    Vector result;
    neighbors(&result, vertex, mode);
    return result;
}

Vector Graph::in_neighbors(long int vertex) const {
    Vector result;
    neighbors(&result, vertex, IGRAPH_IN);
    return result;
}

Vector Graph::out_neighbors(long int vertex) const {
    Vector result;
    neighbors(&result, vertex, IGRAPH_OUT);
    return result;
}

void Graph::setAttribute(const std::string& attribute, const AttributeValue& value) {
    return getAttributeHolder()->setGraphAttribute(attribute, value);
}

void Graph::simplify(bool multiple, bool loops) {
    // TODO: last argument (attribute combination)
    assert(m_pGraph);
    IGRAPH_TRY(igraph_simplify(m_pGraph, multiple, loops, 0));
}

Vertex Graph::vertex(integer_t vid) {
    assert(m_pGraph);
    return Vertex(this, vid);
}

Graph Graph::induced_subgraph(Vector const& allowed_vertices) const {
    // create selector from Vector
    igraph_vs_t* selector;
    IGRAPH_TRY(igraph_vs_vector(selector, allowed_vertices.c_vector()));
    AllocationGuard delete_vs([&]{
        igraph_vs_destroy(selector);
    });

    // use selector to create subgraph
    igraph_t* induced;
    IGRAPH_TRY(igraph_induced_subgraph(m_pGraph, induced, *selector, IGRAPH_SUBGRAPH_AUTO));

    return Graph(induced);
}

/*************/
/* Operators */
/*************/

Graph Graph::operator+(const Graph& other) const {
    std::unique_ptr<igraph_t> result(new igraph_t);
    IGRAPH_TRY(igraph_disjoint_union(result.get(), m_pGraph, other.m_pGraph));
    return Graph(result.release());
}

AttributeValue& Graph::operator[](const std::string& attribute) {
    return getAttributeHolder()->getGraphAttributeReference(attribute);
}

/***************************************************************************/

AttributeHolder* Graph::getAttributeHolder() {
    return static_cast<AttributeHolder*>(m_pGraph->attr);
}

const AttributeHolder* Graph::getAttributeHolder() const {
    return static_cast<AttributeHolder*>(m_pGraph->attr);
}

void Graph::print_edges() const{
    Vector edges = getEdgelist();
    for (auto source_idx = 0; source_idx < edges.size(); source_idx += 2)
        std::cout << edges[source_idx] << "," << edges[source_idx+1] << "  ";
    std::cout << std::endl;
}

    int Graph::isSimple_fast(const igraph_t *graph, igraph_bool_t *res) {
        long int vc = igraph_vcount(graph);
        long int ec = igraph_ecount(graph);

        if (vc == 0 || ec == 0) {
            *res = 1;
        } else {
            igraph_vector_t neis;
            long int i, j, n;
            igraph_bool_t found = 0;
            IGRAPH_VECTOR_INIT_FINALLY(&neis, 0);
            for (i = 0; i < vc; i++) {
                IGRAPH_CHECK(igraph_neighbors(graph, &neis, (igraph_integer_t) i, IGRAPH_OUT));
                n = igraph_vector_size(&neis);
                for (j = 0; j < n; j++) {
                    // changed compared to igraph_is_simple
                    integer_t eid;
                    igraph_get_eid(m_pGraph, &eid, i, j, isDirected(), true);
                    if (deleted_edges.count(eid))
                        continue;
                    // end change
                    if (VECTOR(neis)[j] == i) {
                        found = 1; break;
                    }
                    if (j > 0 && VECTOR(neis)[j - 1] == VECTOR(neis)[j]) {
                        found = 1; break;
                    }
                }
            }
            *res = !found;
            igraph_vector_destroy(&neis);
            IGRAPH_FINALLY_CLEAN(1);
        }

        return 0;
    }

}         // end of namespaces

