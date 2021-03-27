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

sub push_if_new {

    my ($list, $list_hash, $m, $v) = @_;

    if ( exists($list_hash->{$v} ) ) {
        if ( exists($list_hash->{$v}->{$m}) ) {
            return false;   
        }
    } else {
        $list_hash->{$v} = {};
    }

    $list_hash->{$v}->{$m} = 1;
    push @$list, $m;

    return true;
}

# computes all polytopes we may obtain
# from a center
# and stores it if it has not been seen so far
# $c: center information (list of points leading to center with vertices/facets, a list of points we can use for extension, representatives of the cosets of the symmetry of the center vs the center with one point)
sub polytopes_from_center {
    my ($c,$size,$list,$list_hash) = @_;

    foreach my $i (0..scalar(@{$c->{'vertex'}}-1)) {
        foreach my $j ($i..scalar(@{$c->{'vertex'}}-1)) {
            foreach my $auto (@{$c->{'vertex_automorphism_coset'}->[$i]}) {
                my @ret = coning_from_vertices_facets(convert_to<Rational>($c->{'vertices'}->[$i]), convert_to<Rational>($c->{'facets'}->[$i]) , new Vector<Rational>($c->{'vertex'}->[$j]*$auto));
                next unless lattice_point_count($ret[0],$ret[1],$size);

                my $V = new Matrix<Integer>($ret[0]);
                my $F = primitive($ret[1]);
                my $pairing_matrix = $F * transpose($V);

                my @nf = affine_normal_form_from_vertices_and_pairing_matrix(dehomogenize($V), $pairing_matrix, 0);
                my $m = join ",", @{$nf[0]};

                # my $vv = new Vector<Integer>(volume_vector($r->VERTICES));
                my $vv = $V->rows()."-".$F->rows();

                push_if_new($list,$list_hash,$m, $vv);
            }
        }   
    }
}