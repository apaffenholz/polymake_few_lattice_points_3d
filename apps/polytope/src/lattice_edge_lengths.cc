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
#include "polymake/Graph.h"
#include "polymake/common/lattice_tools.h"

namespace polymake { namespace polytope {

namespace {
    Map<Int,Int> lattice_edge_lengths(const BigObject& p) {

        Map<Int,Int> ret;
        Graph<> G = p.give("GRAPH.ADJACENCY");
        Matrix<Rational> V = p.give("VERTICES");

        for ( auto e = entire(edges(G)); !e.at_end(); ++e ) {
            Int v1 = e.from_node();
            Int v2 = e.to_node();
            Vector<Integer> v = convert_to<Vector<Integer>>(V[v1]-V[v2]);
            ret[convert_to<int>(gcd(v))]++;
        }

        return ret;
    }

}

Function4perl(&lattice_edge_lengths, "lattice_edge_lengths(Polytope<Rational>)");

}}