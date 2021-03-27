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
#include "polymake/Vector.h"
#include "polymake/Matrix.h"
#include "polymake/Integer.h"
#include "polymake/Rational.h"


namespace polymake { namespace polytope {

bool lattice_point_count(const Matrix<Rational>& V, const Matrix<Rational>& F, Int s) {

    Vector<Integer> max(V.cols());
    Vector<Integer> min(V.cols()); 
     const Int d = V.cols()-1;
    min[0] = 1;
    max[0] = 1;
    for (Int i = 1; i < V.cols(); ++i ) {
        min[i] = floor(V[0][i]);
        max[i] = ceil(V[0][i]);
    }

    for ( const auto& v: rows(V) ) {
        for (Int i = 1; i < V.cols(); ++i ) {
            if ( v[i] < min[i] ) { min[i] = floor(v[i]); }
            if ( v[i] > max[i] ) { max[i] = ceil(v[i]); }
        }
    }

    Int count = 0;
    Vector<Integer> cur(min);  

    for (;;) {
        bool valid = true;

        if (valid) {
            for (Int i = 0; i < F.rows(); ++i) {
                if (cur * F.row(i) < 0) {
                    valid = false;
                    break;
                }
            }
        }

        if (valid) {
           ++count;
          }

        if ( count > s ) {
            return false;
        }

        Int curd = 1;
        while (curd <= d && cur[curd] >= max[curd]) {
             cur[curd] = min[curd];
            ++curd;
        }

        if (curd <= d)
            cur[curd] += 1;
        else
            break;
    }

    if ( count < s ) {
        return false;
    }
    return true;
}


Function4perl(&lattice_point_count, "lattice_point_count(Matrix<Rational>,Matrix<Rational>,Int)");

}}