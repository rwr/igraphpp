/* vim:set ts=4 sw=4 sts=4 et: */

#ifndef IGRAPHPP_ISOMORPHISM_H
#define IGRAPHPP_ISOMORPHISM_H

namespace igraph {

class Graph;

/// Checks if two given (di)graphs are isomorphic. A wrapper around isomorphic(...).
bool isomorphic(const Graph& graph1, const Graph& graph2);

}         // end of namespace

#endif    // IGRAPHPP_ISOMORPHISM_H
