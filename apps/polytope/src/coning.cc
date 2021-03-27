/* Copyright (c) 2020-2021
   Andreas Paffenholz
   Technische Universität Berlin, Germany
   https://www2.mathematik.tu-darmstadt.de/~paffenholz

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
--------------------------------------------------------------------------------
*/

#include "polymake/client.h"
#include "polymake/Array.h"
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/linalg.h"


namespace polymake { namespace polytope {

void coning_impl(const Matrix<Rational>& V, const Matrix<Rational>& F, const Vector<Rational>& v, Matrix<Rational>& Vnew, Matrix<Rational>& Fnew) {
    
    Set<Int> new_verts;
    IncidenceMatrix<> ftv(V.rows(), F.rows(), attach_operation(product(rows(V), rows(F), operations::mul()), operations::is_zero()).begin());
    Int n_verts = V.rows();
    for ( Int s=0; s < n_verts; ++s ) {
        for ( auto f : ftv[s] ) {
            if ( F[f]*v > 0 ) {
                new_verts += s;
                break;
            }
        }
    }

    Set<Int> remaining_facets;
    Set<Int> remaining_facets_containing_v;
    Set<Int> discarded_facets;
    for ( Int i = 0; i < F.rows(); ++i ) {
        if ( F[i]*v >= 0 ) {
            if ( F[i] * v > 0 ) {
                remaining_facets += i;
            } else {
                remaining_facets_containing_v += i;
            }
        } else {
            discarded_facets += i;
        }
    }

    Vector<Rational> interior_point = accumulate(rows(V),operations::add())/n_verts;
    Matrix<Rational> new_facets(0,F.cols());
    for ( const auto fd: discarded_facets ) {
        for ( const auto fr: remaining_facets ) {
            Set<Int> ridge;
            for ( Int s = 0; s < n_verts; ++ s ) {
                if ( F[fd] * V[s] == 0 && F[fr] * V[s] == 0 ) {
                    ridge += s;
                }
            }
            Matrix<Rational> verts = V.minor(ridge,All)/v;
            Matrix<Rational> ns = null_space(verts);
            if ( ns.rows() == 1 ) {
                Vector<Rational> nf(ns[0]);
                new_facets /= ns[0] * interior_point > 0 ? nf : -nf;
            }
        }
    }

    Vnew = V.minor(new_verts,All)/v;
    Fnew = F.minor(remaining_facets+remaining_facets_containing_v,All)/new_facets;

} 

BigObject coning(const BigObject& p, const Vector<Rational>& v) {
    
    Matrix<Rational> F = p.give("FACETS");
    Matrix<Rational> V = p.give("VERTICES");

    Matrix<Rational> Vnew;
    Matrix<Rational> Fnew;
    coning_impl(V,F,v,Vnew,Fnew);

    BigObject q("Polytope<Rational>");
    q.take("VERTICES") << Vnew;
    q.take("FACETS") << Fnew;
    return q;
}

ListReturn coning_from_vertices_facets(const Matrix<Rational>& V, const Matrix<Rational>& F, const Vector<Rational>& v) {

    Matrix<Rational> Vnew;
    Matrix<Rational> Fnew;
    coning_impl(V,F,v,Vnew,Fnew);

    ListReturn ret;
    ret << Vnew;
    ret << Fnew;
    return ret;  
}


Function4perl(&coning, "coning(Polytope<Rational>,Vector<Rational>)");
Function4perl(&coning_from_vertices_facets, "coning_from_vertices_facets(Matrix<Rational>,Matrix<Rational>,Vector<Rational>)");
}}