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
#include "polymake/PowerSet.h"
#include "polymake/linalg.h"

namespace polymake { namespace polytope {
namespace {
    void smaller_volume_vector(const Matrix<Rational>& m, Vector<Rational>& v, const Array<Int>& perm, const Subsets_of_k<const Series<Int,true> >& subsets, const Map<Array<Int>,Rational>& dets) {

        int i = 0;
        bool undecided = true;
        bool first = true;
        int sign = 1;
        for ( auto s = entire(subsets); !s.at_end(); ++s, ++i ) {
            Array<Int> sel(m.cols());
            int j = 0;
            for ( auto k = entire(*s); !k.at_end(); ++k, ++j ) {
                sel[j] = perm[*k];
            }
            int psign = 1;
            for ( int l1 = 0; l1 < m.cols()-1; ++l1 ) {
                for ( int l2 = 0; l2 < m.cols()-1; ++l2 ) {
                    if ( sel[l2] > sel[l2+1] ) {
                        Int temp = sel[l2];
                        sel[l2] = sel[l2+1];
                        sel[l2+1] = temp;
                        psign *= -1;
                    }
                }
            }


            Rational d = psign * dets[sel];
            if ( first && v[0] >= -d ) {
                sign = -1;
            }
            first = false;
            if ( v[i] == sign*d ) {
                continue;
            }
            if ( undecided ) {
                if ( sign * d > v[i] ) {
                    return;
                } else {
                    undecided = false;
                }
            }
            v[i] = sign *d;
        }
    }


    Vector<Rational> volume_vector( const Matrix<Rational>& m) {

        Map<Array<Int>,Rational> dets;
        Subsets_of_k<const Series<Int,true> > subsets = all_subsets_of_k(range(0,m.rows()-1),m.cols());
        Integer n = Integer::binom(m.rows(),m.cols());
        Vector<Rational> v(convert_to<int>(n));
        int i = 0;
        for ( auto s = entire(subsets); !s.at_end(); ++s, ++i ) {
            Array<Int> sa(*s);
            dets[sa] = det(m.minor(*s,All));
            v[i] = dets[sa];
        }

        AllPermutations<> perms(m.rows());
        for (AllPermutations<>::const_iterator perm=entire(perms);  !perm.at_end();  ++perm) {
            smaller_volume_vector(m,v,*perm,subsets,dets);
        }

        return v;
    }

}

Function4perl(&volume_vector, "volume_vector");

} }