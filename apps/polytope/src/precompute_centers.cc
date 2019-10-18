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
namespace {

    std::vector<Matrix<Rational>> find_all_equivalences(const perl::Object& q) {

        std::vector<Matrix<Rational>> automorphisms;
        int dim = q.give("CONE_DIM");
        Matrix<Rational> V = q.give("VERTICES");
        Set<int> basis = basis_rows(V);
        Matrix<Rational> basis_mat = V.minor(basis,All);
        Set<Vector<Rational>> q_vertex_set(rows(V));

        Subsets_of_k<const Series<int,true> > subsets = all_subsets_of_k(range(0,V.rows()-1),dim);
        AllPermutations<> perms(dim);

        for ( auto s = entire(subsets); !s.at_end(); ++s ) {
            Array<int> sa(*s);
            for (AllPermutations<>::const_iterator perm=entire(perms);  !perm.at_end();  ++perm) {
                Matrix<Rational> C;
                for ( auto e : *perm ) {
                    C /= V[sa[e]];
                }
                if ( det(C) == 0 ) {
                    continue;
                }
                Matrix<Rational> t = inv(basis_mat)*C;
                if ( t == polymake::common::primitive(t) && abs(det(t)) == 1 ) {
                    perl::Object p("Polytope<Rational>");
                    p.take("POINTS") << V*t;
                    Matrix<Rational> Vp = p.give("VERTICES");
                    Set<Vector<Rational>> p_vertex_set(rows(Vp));
                    if ( incl(q_vertex_set, p_vertex_set) == 0 ) {
                        automorphisms.push_back(t);
                    } 
                }
            }
        }

        return automorphisms;
    }

    bool check_lattice_equivalence ( const perl::Object& p, const perl::Object& q, Matrix<Rational>& trafo) {

        int dim = p.give("CONE_DIM");
        int qdim = q.give("CONE_DIM");
        if ( dim != qdim ) return false;

        int p_n_verts = p.give("N_VERTICES");
        int q_n_verts = q.give("N_VERTICES");
        if ( q_n_verts != p_n_verts ) return false;

        int p_n_facets = p.give("N_FACETS");
        int q_n_facets = q.give("N_FACETS");
        if ( p_n_facets != q_n_facets ) return false;

        int p_n_edges = p.give("N_EDGES");
        int q_n_edges = q.give("N_EDGES");
        if ( p_n_edges != q_n_edges ) return false;

        int p_lv = p.give("LATTICE_VOLUME");
        int q_lv = q.give("LATTICE_VOLUME");
        if ( p_lv != q_lv ) return false;

        int p_ilp = p.give("N_INTERIOR_LATTICE_POINTS");
        int q_ilp = q.give("N_INTERIOR_LATTICE_POINTS");
        if ( p_ilp != q_ilp ) return false;

        const IncidenceMatrix<> Mp=p.give("RAYS_IN_FACETS"), Mq=q.give("RAYS_IN_FACETS");
        if ( !graph::isomorphic(Mp,Mq) ) return false;

        Matrix<Rational> V = p.give("VERTICES");
        Set<int> basis = basis_rows(V);
        Matrix<Rational> basis_mat = V.minor(basis,All);
        Set<Vector<Rational>> p_vertex_set(rows(V));

        Subsets_of_k<const Series<int,true> > subsets = all_subsets_of_k(range(0,V.rows()-1),dim);
        AllPermutations<> perms(dim);

        for ( auto s = entire(subsets); !s.at_end(); ++s ) {
            Array<int> sa(*s);

            for (AllPermutations<>::const_iterator perm=entire(perms);  !perm.at_end();  ++perm) {
                Matrix<Rational> Vq = q.give("VERTICES");
                Matrix<Rational> C;

                for ( auto e : *perm ) { 
                    C /= Vq[sa[e]];
                }

                if ( det(C) == 0 ) {
                    continue;
                }

                if ( det(C) == 0 ) continue;
                Matrix<Rational> t = inv(C)*basis_mat;
                if ( t == polymake::common::primitive(t) && abs(det(t)) == 1 ) {
                    Set<Vector<Rational>> q_vertex_set(rows(Vq*t));

                    if ( incl(q_vertex_set, p_vertex_set) == 0 ) {
                        trafo = t;
                        return true;
                    } 
                }
            }
        }

        return false;
    }

    std::pair<int,Matrix<Rational>> find_in_centers( perl::Object& pred, std::vector<perl::Object>& centers, std::vector<std::vector<Matrix<Rational>>>& center_automorphisms ) {

        std::pair<int,Matrix<Rational>> ret;
        int dim = pred.give("CONE_DIM");

        Matrix<Rational> trafo;
        for ( unsigned int i = 0; i < centers.size(); i++ ) {
            if ( check_lattice_equivalence(pred,centers[i],trafo) ) {
                ret.first = i;
                ret.second = inv(trafo);
                return ret;
            }
        }
        centers.push_back(pred);
        center_automorphisms.push_back(find_all_equivalences(pred));
        ret.second = unit_matrix<Rational>(dim);
        ret.first = centers.size()-1;
    
        return ret;
    } 


    bool is_unit_matrix ( const Matrix<Rational>& m ) {
        int c = m.cols();
        if ( m.rows() != c ) return false;
        for ( int i = 0; i < c; i++ ) {
            for ( int j = 0; j < c; j++ ) {
                if ( ( i != j && m(i,j) != 0 ) || ( i == j && m(i,j) != 1 ) ) return false;
            }
        }
        return true;
    }

}

/*
 * a: a list of all lattice polytopes with n-1 lattice points, where we want to compute those with n
 * centers: an array of all possible centers, up to lattice equivalence
 * polys: an array that records a map for each polytope. 
 *      The map gives for each vertex the corresponding center and the lattice transformation needed to make it equal to the center
 *      We reduce up to symmetry of the polytope, so not all vertices appear in the hash
 * center_automorphisms: Automorphisms for each center
 */

perl::ListReturn precompute_centers( const Array<perl::Object>& a ) {

    std::vector<perl::Object> centers;
    std::vector<Map<int,std::pair<int,Matrix<Rational>>>> polys;
    std::vector<std::vector<Matrix<Rational>>> center_automorphisms;
    int asize = a.size();
    int index = 0;

    for ( const perl::Object p : a ) {
        cout << "working on " << index++ << " of " << asize << endl;
        int dim = p.give("CONE_DIM");
        int n_verts = p.give("N_VERTICES");
        Matrix<Rational> V = p.give("VERTICES");
        Set<Vector<Rational>> p_vertex_set(rows(V));

        Matrix<Rational> plp(0,dim);
        Array<Matrix<Rational>> lp = p.give("LATTICE_POINTS_GENERATORS");
        for (auto v : rows(lp[0]) ) {
            if ( !p_vertex_set.contains(v)) {
                plp /= v;
            }
        }

        Map<int,std::pair<int,Matrix<Rational>>> pmap;

        for ( int pv = 0; pv < n_verts; pv++ ) {
            Matrix<Rational> V_pred = V.minor(~range(pv,pv),All);
            if ( plp.rows() ) {
                V_pred /= plp; 
            }
            perl::Object pred("Polytope<Rational>");
            pred.take("POINTS") << V_pred;
            int dim_pred = pred.give("CONE_DIM");
            if ( dim_pred < dim ) {
                continue;
            }
            std::pair<int,Matrix<Rational>> ret = find_in_centers(pred, centers, center_automorphisms);
            if ( !is_unit_matrix(ret.second) && V[pv]*ret.second == V[pv] ) {
                cout << "skipped one" << endl;
            }
            pmap[pv] = ret;
        }

        polys.push_back(pmap);
    }
  perl::ListReturn result;
  result << centers << polys << center_automorphisms;
  return result;
}


bool lattice_equivalent(const perl::Object& p, const perl::Object&q ) {
    Matrix<Rational>t;

    return check_lattice_equivalence(p,q,t);
}

Array<Matrix<Rational>> lattice_automorphisms ( const perl::Object& p)  {
    return Array<Matrix<Rational>>(find_all_equivalences(p));
}

Function4perl(&precompute_centers, "precompute_centers(Array<Polytope>)");

Function4perl(&lattice_equivalent, "lattice_equivalent(Polytope,Polytope)");

Function4perl(&lattice_automorphisms, "lattice_automorphisms(Polytope)");

} }