#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/permutations.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/common/lattice_tools.h"
#include "polymake/graph/compare.h"


namespace polymake { namespace polytope {

perl::Object coning(const perl::Object& p, const Vector<Rational>& v) {
    
    Set<int> new_verts;
    int n_verts = p.give("N_VERTICES");
    IncidenceMatrix<> ftv = p.give("FACETS_THRU_VERTICES");
    Matrix<Rational> F = p.give("FACETS");
    Matrix<Rational> V = p.give("VERTICES");
    for ( int s=0; s < n_verts; ++s ) {
        for ( auto f : ftv[s] ) {
            if ( F[f]*v > 0 ) {
                new_verts += s;
                break;
            }
        }
    }

    Matrix<Rational> Vn = V.minor(new_verts,All)/v;
    perl::Object q("Polytope<Rational>");
    q.take("VERTICES") << Vn;
    return q;
} 


Function4perl(&coning, "coning(Polytope<Rational>,Vector<Rational>)");

}}