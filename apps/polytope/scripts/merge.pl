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

# compute all merged polytopes given a list of centers
# we run over all centers and try all possible extensions
# we use normal form to test for equivalence with a polytope 
# that we already found

use application "polytope";
use Term::ANSIColor;
use Time::HiRes qw( time );
use DateTime;

script("polytope_centers.pl");
script("polytopes_from_center.pl");

sub merge_impl {
    my ($centers,$size) = @_;
    my $total_centers = scalar(@$centers);

    my $list = [];
    my $list_hash = {};

    foreach my $c (@$centers) {
        print "starting center ".$c->{'index'}."/".$total_centers."\n";
        polytopes_from_center($c,$size,$list,$list_hash);
    }
    my $polytope_list = [];
    my $count = 1;
    foreach my $ms (@$list) {
        my @entries = split / /, $ms;
        my $rows = new Int(scalar(@entries)/3);
        my $k = 0;
        my $m = new Matrix<Rational>($rows,3);
        foreach my $i (0..$rows-1) {
            foreach my $j (0..2) {
                $m->elem($i,$j) = $entries[$k++];
            }
        }
        my $p = new Polytope(VERTICES=>(ones_vector<Rational>($rows))|$m);
        my $name = sprintf "merged.s%02d.%0*d", $size, $size-4, $count;
        $p->name = $name;
        $count++;
        push @$polytope_list, $p;
    }
    return $polytope_list;
}


sub merge {
    my ($a_in) = @_;

    my $total_start = DateTime->now();

    my $asize = scalar(@$a_in)-1;
    my $s = new Int($a_in->[0]->N_LATTICE_POINTS+1);

    my $listmap = new Map<Vector<Integer>,Set<Int>>;
    
    my $centers = [];
    print "preparing centers ...  \n";
    $centers = centers($a_in);
    print "   ... done\n";
    my $total_centers = scalar(@$centers);
    print "working on ".$total_centers." centers\n";

    my $polytope_list = [];
    $polytope_list = merge_impl($centers,$s);

    print "\n", $total_start, "\n", DateTime->now(), "\n";
    return $polytope_list;
}
