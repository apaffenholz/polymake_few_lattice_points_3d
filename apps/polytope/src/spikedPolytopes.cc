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


std::vector<perl::Object> spiked_polytopes(int n) {

    std::vector<perl::Object> spiked_polytopes;

    return spiked_polytopes;
}

Function4perl(&spiked_polytopes, "spiked_polytopes($)");

}}