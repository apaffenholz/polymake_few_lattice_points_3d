#  Copyright (c) 2020-2021
#  Andreas Paffenholz
#  Technische Universität Darmstadt, Germany
#  https://www2.mathematik.tu-darmstadt.de/~paffenholz
#
#  This program is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation; either version 2, or (at your option) any
#  later version: http://www.gnu.org/licenses/gpl.txt.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#-------------------------------------------------------------------------------


use application "polytope";

script("boxed_and_spiked_polytopes.pl");
script("merge.pl");


# all polytopes of given size and width at least 2
# based on the methods in
# Mónica Blanco and Francisco Santos
# On the enumeration of lattice 3-polytopes
# arxiv: 1601.02577
#
# $size: the size to be computed
# $col: a list of polytopes of size $size-1
sub generate_polytopes_of_size {
    my ($size, $coll) = @_;

    if ( $size == 5 ) {
        return size_5_polytopes();
    }

    if ( $size == 6 ) {
        return size_6_polytopes();
    }

    my $polys_in = [];
    foreach my $p (@$coll) {
        push @$polys_in, $p;
    }

    print "merged polytopes ... ";
    my $polys = merge($polys_in, $options);
    print "... done\n spiked polytopes ... ";
    my $spiked = spiked_polytopes($size);
    foreach my $p (@$spiked) {
        push @$polys, $p;
    }
    print "... done\nboxed quasiminial polytopes ... ";
    my $boxed = boxed_quasi_minimal_polytopes($size);
    foreach my $p (@$boxed) {
        push @$polys, $p;
    }
    print " ... done\n";

    return $polys;
}