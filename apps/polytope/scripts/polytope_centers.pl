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

# computes for a list of polytopes of the same size s
# a list of centers c in affine normal form
# together with
# a list of vertices v that extend a center to a polytope P_v = conv(c,v) of size s, up to lattice equivalence
# lists of vertices/facets for P_v
# a list of vertices we can add to P_v (vertices of some polytope of size s that lead to the same center if removed)
# a representative of each coset of equivalences of the c vs P_v
#
# two versions: one single threaded, one with simple parallelization

use application "polytope";

sub centers_list_of_polytopes {
    my $a = shift; 
    return if !scalar(@$a);

    my $centers = [];
    my $center_hash = {};
    my $count = 0;

    foreach my $p (@$a) {
        my $dim = $p->DIM;
        my $p_n_verts = $p->N_VERTICES;

        # the vertices of $p as set, to check containment
        my $p_vertex_set = new Set<Vector<Rational>>(@{$p->VERTICES});

        # the non-vertex lattice points
        my $non_vertex_lattice_points = new Matrix<Rational>(0,$dim+1);
        foreach (@{$p->LATTICE_POINTS}) {
            if ( !$p_vertex_set->contains($_)) {
                $non_vertex_lattice_points /= $_;
            }
        }

        my $V = $p->VERTICES;

        foreach my $pv (0..$p_n_verts-1) { 
            my $pred = $non_vertex_lattice_points->rows ? 
                            new Polytope(POINTS=>$p->VERTICES->minor(~range($pv,$pv),All)/$non_vertex_lattice_points) : 
                            new Polytope(POINTS=>$p->VERTICES->minor(~range($pv,$pv),All));

            next if ( $pred->DIM < 3 );

            my @nf_with_info = affine_normal_form($pred);

            # the hermite normal form
            my $center =  ones_vector<Integer>($nf_with_info[0]->rows())|$nf_with_info[0];
            my $center_string = "";
            for my $r (@$center) {
                for my $e (@$r) {
                    $center_string .= $e." ";
                }
            }
            # the transformation leading to this hnf
            my $transformation = unit_vector<Integer>(4,0)/(zero_vector<Integer>(3)|$nf_with_info[1]);
            # the shifted and transformed additional vertex
            my $v = new Vector<Integer>(($p->VERTICES->[$pv] + (0|$nf_with_info[2])) * $transformation);

            # the shifted and transformed full list of vertices
            my $V = new Matrix<Rational>( translate($p, new Vector<Rational>($nf_with_info[2]))->VERTICES * $transformation );

            # the shifted and transformed full list of facets
            my $F = primitive(new Matrix<Rational>(  translate($p, new Vector<Rational>($nf_with_info[2]))->FACETS * transpose(inv($transformation)) ));

            my $center_automorphisms = matrix_automorphisms(new Matrix<Rational>($center));
            my $center_automorphisms_set = new Set<Matrix<Integer>>;
            foreach my $ca (@$center_automorphisms) {
                $center_automorphisms_set += $ca;
            }

            if ( exists($center_hash->{$center_string}) ) {

                my $index = $center_hash->{$center_string};
                my $seen = 0;
                foreach my $w (@{$centers->[$index]->{'vertex'}}) {
                    if ( $v == $w ) {
                        $seen = 1;
                    }
                }
                next if $seen == 1;

                my $center_automorphisms_subgroup = matrix_automorphisms(new Matrix<Rational>($center/$v));
                my $center_automorphisms_set_copy = new Set<Matrix<Integer>>($center_automorphisms_set);
                my $center_automorphisms_coset = [];
                while ($center_automorphisms_set_copy->size() != 0 ) {
                    my $auto = $center_automorphisms_set_copy->[0];
                    push @$center_automorphisms_coset, $auto;
                    foreach my $cas (@$center_automorphisms_subgroup) {
                        $center_automorphisms_set_copy -= $auto*$cas;
                    }
                }
                push @{$centers->[$index]->{'vertex'}}, $v;
                push @{$centers->[$index]->{'vertices'}}, $V;
                push @{$centers->[$index]->{'facets'}}, $F;
                push @{$centers->[$index]->{'vertex_automorphism_coset'}}, $center_automorphisms_coset;
            } else {
                my $index = $count++;
                $center_hash->{$center_string} = $index;
                $centers->[$index]->{'center_string'} = $center_string;
                $centers->[$index]->{'index'} = $index;
                $centers->[$index]->{'center'} = $center;
                $centers->[$index]->{'automorphisms'} = $center_automorphisms;

                my $center_automorphisms_subgroup = matrix_automorphisms(new Matrix<Rational>($center/$v));
                my $center_automorphisms_set_copy = new Set<Matrix<Integer>>($center_automorphisms_set);
                my $center_automorphisms_coset = [];
                while ($center_automorphisms_set_copy->size() != 0 ) {
                    my $auto = $center_automorphisms_set_copy->[0];
                    push @$center_automorphisms_coset, $auto;
                    foreach my $cas (@$center_automorphisms_subgroup) {
                        $center_automorphisms_set_copy -= $auto*$cas;
                    }
                }  
                $centers->[$index]->{'vertex'} = [];
                push @{$centers->[$index]->{'vertex'}}, $v;
                $centers->[$index]->{'vertices'} = [];
                push @{$centers->[$index]->{'vertices'}}, $V;
                $centers->[$index]->{'facets'} = [];
                push @{$centers->[$index]->{'facets'}}, $F;
                $centers->[$index]->{'vertex_automorphism_coset'} = [];
                push @{$centers->[$index]->{'vertex_automorphism_coset'}}, $center_automorphisms_coset;
            }
        }
    }
    foreach my $center (@$centers) {
        delete($center->{'automorphisms'});
        delete($center->{'center'});
    }
    return $centers;
}


sub centers_to_file {
    my $in_name = shift;
    my $out_name = $in_name =~ s/merged/partial_centers/r;

    my $list = load_data($in_name);
    my $centers = centers_list_of_polytopes($list);
    save_data($centers,$out_name);
}

sub merge_centers_from_file {
    my ($dir,$prefix) = @_;
    my @files = glob( $dir . '/'.$prefix.'*' );
    
    my $centers = [];
    my $center_hash = {};
    my $count = 0;
    my $i = 0;
    foreach my $file (@files) {
        my $new_centers = load_data($file);
        print "merging ".$i."\n"; $i++;
        foreach my $center_data (@$new_centers) {

            my $center_string = $center_data->{'center_string'};            
            
            if ( exists($center_hash->{$center_string}) ) {
                my $index = $center_hash->{$center_string};

                foreach my $j (0..scalar(@{$center_data->{'vertex'}})-1 ) {
                    my $seen = 0;
                    my $vp = new Vector<Integer>($center_data->{'vertex'}->[$j]);
                    foreach my $w (@{$centers->[$index]->{'vertex'}}) {
                        if ( $vp == $w ) {
                            $seen = 1;
                        }
                    }
                    next if $seen == 1;

                    push @{$centers->[$index]->{'vertex'}}, $vp;
                    push @{$centers->[$index]->{'vertices'}}, new Matrix<Integer>($center_data->{'vertices'}->[$j]);
                    push @{$centers->[$index]->{'facets'}}, new Matrix<Rational>($center_data->{'facets'}->[$j]);
                    push @{$centers->[$index]->{'vertex_automorphism_coset'}}, new Array<Matrix<Integer>>($center_data->{'vertex_automorphism_coset'}->[$j]);
                }
            } else {
                my $index = $count++;
                $center_hash->{$center_string} = $index;
                $centers->[$index]->{'index'} = $index;

                $centers->[$index]->{'vertex'} = [];
                $centers->[$index]->{'vertices'} = [];
                $centers->[$index]->{'facets'} = [];
                $centers->[$index]->{'vertex_automorphism_coset'} = [];
                foreach my $j (0..scalar(@{$center_data->{'vertex'}})-1 ) {
                    push @{$centers->[$index]->{'vertex'}}, new Vector<Integer>($center_data->{'vertex'}->[$j]);
                    push @{$centers->[$index]->{'vertices'}}, new Matrix<Integer>($center_data->{'vertices'}->[$j]);
                    push @{$centers->[$index]->{'facets'}}, new Matrix<Rational>($center_data->{'facets'}->[$j]);
                    push @{$centers->[$index]->{'vertex_automorphism_coset'}}, new Array<Matrix<Integer>>($center_data->{'vertex_automorphism_coset'}->[$j]);
                }   
            }
        }
    }

    my $out_name = $prefix =~ s/partial_centers/centers/r;
    save_data($centers,$dir.'/'.$prefix.".pdata");
}

sub centers {

    my $a = shift;  # the polytopes with n-1 lattice points
    return if !scalar(@$a);

    my $centers = [];
    my $center_hash = {};
    my $count = 0;

    foreach my $p (@$a) {
        my $localcount = 0;
        my $dim = $p->DIM;
        my $p_n_verts = $p->N_VERTICES;

        # the vertices of $p as set, to check containment
        my $p_vertex_set = new Set<Vector<Rational>>(@{$p->VERTICES});

        # the non-vertex lattice points
        my $non_vertex_lattice_points = new Matrix<Rational>(0,$dim+1);
        foreach (@{$p->LATTICE_POINTS}) {
            if ( !$p_vertex_set->contains($_)) {
                $non_vertex_lattice_points /= $_;
            }
        }

        my $V = $p->VERTICES;

        foreach my $pv (0..$p_n_verts-1) { 
            my $pred = $non_vertex_lattice_points->rows ? 
                            new Polytope(POINTS=>$p->VERTICES->minor(~range($pv,$pv),All)/$non_vertex_lattice_points) : 
                            new Polytope(POINTS=>$p->VERTICES->minor(~range($pv,$pv),All));

            next if ( $pred->DIM < 3 );

            my @nf_with_info = affine_normal_form($pred);

            # the hermite normal form
            my $center =  ones_vector<Integer>($nf_with_info[0]->rows())|$nf_with_info[0];
            my $center_string = "";
            for my $r (@$center) {
                for my $e (@$r) {
                    $center_string .= $e." ";
                }
            }
            # the transformation leading to this hnf
            my $transformation = unit_vector<Integer>(4,0)/(zero_vector<Integer>(3)|$nf_with_info[1]);
            # the shifted and transformed additional vertex
            my $v = new Vector<Integer>(($p->VERTICES->[$pv] + (0|$nf_with_info[2])) * $transformation);

            # the shifted and transformed full list of vertices
            my $V = new Matrix<Rational>( translate($p, new Vector<Rational>($nf_with_info[2]))->VERTICES * $transformation );

            # the shifted and transformed full list of facets
            my $F = primitive(new Matrix<Rational>(  translate($p, new Vector<Rational>($nf_with_info[2]))->FACETS * transpose(inv($transformation)) ));


            my $center_automorphisms = matrix_automorphisms(new Matrix<Rational>($center));
            my $center_automorphisms_set = new Set<Matrix<Integer>>;
            foreach my $ca (@$center_automorphisms) {
                $center_automorphisms_set += $ca;
            }

            if ( exists($center_hash->{$center_string}) ) {

                my $index = $center_hash->{$center_string};
                my $seen = 0;
                foreach my $w (@{$centers->[$index]->{'vertex'}}) {
                    if ( $v == $w ) {
                        $seen = 1;
                    }
                }
                next if $seen == 1;

                $localcount++;
                my $center_automorphisms_subgroup = matrix_automorphisms(new Matrix<Rational>($center/$v));
                my $center_automorphisms_set_copy = new Set<Matrix<Integer>>($center_automorphisms_set);
                my $center_automorphisms_coset = [];
                while ($center_automorphisms_set_copy->size() != 0 ) {
                    my $auto = $center_automorphisms_set_copy->[0];
                    push @$center_automorphisms_coset, $auto;
                    foreach my $cas (@$center_automorphisms_subgroup) {
                        $center_automorphisms_set_copy -= $auto*$cas;
                    }
                }
                push @{$centers->[$index]->{'vertex'}}, $v;
                push @{$centers->[$index]->{'vertices'}}, $V;
                push @{$centers->[$index]->{'facets'}}, $F;
                push @{$centers->[$index]->{'vertex_automorphism_coset'}}, $center_automorphisms_coset;
            } else {
                my $index = $count++;
                $localcount++;
                $center_hash->{$center_string} = $index;
                $centers->[$index]->{'center_string'} = $center_string;
                $centers->[$index]->{'index'} = $index;
                $centers->[$index]->{'center'} = $center;
                $centers->[$index]->{'automorphisms'} = $center_automorphisms;

                my $center_automorphisms_subgroup = matrix_automorphisms(new Matrix<Rational>($center/$v));
                my $center_automorphisms_set_copy = new Set<Matrix<Integer>>($center_automorphisms_set);
                my $center_automorphisms_coset = [];
                while ($center_automorphisms_set_copy->size() != 0 ) {
                    my $auto = $center_automorphisms_set_copy->[0];
                    push @$center_automorphisms_coset, $auto;
                    foreach my $cas (@$center_automorphisms_subgroup) {
                        $center_automorphisms_set_copy -= $auto*$cas;
                    }
                }
                $centers->[$index]->{'vertex'} = [];
                push @{$centers->[$index]->{'vertex'}}, $v;
                $centers->[$index]->{'vertices'} = [];
                push @{$centers->[$index]->{'vertices'}}, $V;
                $centers->[$index]->{'facets'} = [];
                push @{$centers->[$index]->{'facets'}}, $F;
                $centers->[$index]->{'vertex_automorphism_coset'} = [];
                push @{$centers->[$index]->{'vertex_automorphism_coset'}}, $center_automorphisms_coset;
            }
        }
        print $p->name." contributed ".$localcount." entries\n";
    }
    return $centers;
}

sub vector_to_perl { my $v = shift; my $vp = []; foreach (@$v) { push @$vp, convert_to<Int>($_); } return $vp; }
sub matrix_to_perl { my $m = shift; my $mp = []; foreach (@$m) { push @$mp, \@{vector_to_perl($_)}; } return $mp; }




