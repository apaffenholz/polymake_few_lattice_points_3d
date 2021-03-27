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
#include "polymake/permutations.h"
#include "polymake/linalg.h"
#include "polymake/PowerSet.h"


namespace polymake { namespace polytope {

    Array<Matrix<Integer>> matrix_automorphisms(const Matrix<Rational>& m) {

        std::vector<Matrix<Integer> > a;

        Matrix<Rational> b = m.minor(basis_rows(m),All);
        Set<Vector<Rational> > ms(rows(m));
        for (auto s = entire(all_subsets_of_k(sequence(0, m.rows()), m.cols())); !s.at_end(); ++s) {
           Matrix<Rational> n = m.minor(*s,All);
            if ( det(n) == 0 ) {
                continue;
            }
            for (AllPermutations<>::const_iterator perm=entire(all_permutations(m.cols()));  !perm.at_end();  ++perm) {
                Matrix<Rational> t = inv(b)*n.minor(Array<Int>(*perm), All);
                bool is_integral = true;
                for ( const auto& r : rows(t) ) {
                    for ( const auto& e : r ) {
                        if ( abs(denominator(e)) != 1 ) {
                            is_integral = false;
                            break;
                        }
                        if ( !is_integral ) {
                            break;
                        }
                    }
                }
                if ( !is_integral ) {
                    continue;
                }
                if ( abs(det(t)) != 1 ) {
                    continue;
                }
                Set<Vector<Rational> > ms_t(rows(m*t));
                if ( ms_t == ms ) {
                    a.push_back(Matrix<Integer>(t));
                }
            }
        }

        return Array<Matrix<Integer> >(a);
    }

    Function4perl(&matrix_automorphisms, "matrix_automorphisms(Matrix<Rational>)");

}}