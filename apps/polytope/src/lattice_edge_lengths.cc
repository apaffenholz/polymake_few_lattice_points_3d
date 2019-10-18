#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Graph.h"
#include "polymake/Rational.h"
#include "polymake/permutations.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/common/lattice_tools.h"
#include "polymake/graph/compare.h"


namespace polymake { namespace polytope {

Map<int,int> lattice_edge_lengths(const perl::Object& p) {

    Map<int,int> ret;
    Graph<> G = p.give("GRAPH.ADJACENCY");
    Matrix<Rational> V = p.give("VERTICES");

    for ( auto e = entire(edges(G)); !e.at_end(); ++e ) {
        int v1 = e.from_node();
        int v2 = e.to_node();
        Vector<Integer> v = convert_to<Vector<Integer>>(V[v1]-V[v2]);
        ret[convert_to<int>(gcd(v))]++;
    }

    return ret;
} 


Function4perl(&lattice_edge_lengths, "lattice_edge_lengths(Polytope<Rational>)");

}}