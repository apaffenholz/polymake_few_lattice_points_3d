use application "polytope";
use Term::ANSIColor;
use Time::HiRes qw( time );
use DateTime;

sub find_all_equivalences {
    my ($q) = @_;

    my $automorphisms = ();
    my $dim = $q->DIM;
    my $basis = basis_rows($q->VERTICES);
    my $basis_mat = $q->VERTICES->minor($basis,All);
    my $q_vertex_set = new Set<Vector<Rational>>(@{$q->VERTICES});

    foreach  my $subset (@{all_subsets_of_k(range(0,$q->N_VERTICES-1),$dim+1)}) {
        foreach my $perm (@{all_permutations($dim+1)}) {
            my @subset_arr = @$subset[@$perm];
            my $target = $q->VERTICES->minor(\@subset_arr, All);
            next if det($target) == 0;
            my $t = inv($basis_mat)*$target;
            if ( $t == primitive($t) && abs(det($t)) == 1) {
                my $r=new Polytope(POINTS=>$q->VERTICES*$t);
                my $r_vertex_set = new Set<Vector<Rational>>(@{$r->VERTICES});
                if ( incl($q_vertex_set,$r_vertex_set) == 0 ) {
                    push @$automorphisms, $t;
                }
            }
        }
    }
    return $automorphisms;
}

sub find_in_centers {

    my ($p,$centers,$center_automorphisms) = @_;
    my $ret =();
    foreach (0..scalar(@$centers)-1) {
        my $r = check_center_lattice_equivalence($p,$centers->[$_]);
        if ( $r->{equiv} ) {
            $ret->{index} = $_;
            $ret->{trafo} = $r->{trafo};
            return $ret;
        }
    }
    push @$centers, $p;
    push @$center_automorphisms, find_all_equivalences($p);
    $ret->{trafo} = unit_matrix($p->DIM+1);
    $ret->{index} = scalar(@$centers)-1;
    return $ret;
}

sub prepare_centers {

    my $a = shift;

    my $polys = [];
    my $centers = [];
    my $center_automorphisms = [];

    foreach my $p (@$a) {
        my $dim = $p->DIM;
        my $p_n_verts = $p->N_VERTICES;

        # the vertices of $p as set, to check containment
        my $p_vertex_set = new Set<Vector<Rational>>(@{$p->VERTICES});

        # the non-vertex lattice points
        my $plp = new Matrix<Rational>(0,$dim+1);
        foreach (@{$p->LATTICE_POINTS}) {
            if ( !$p_vertex_set->contains($_)) {
                $plp /= $_;
            }
        }
        my $pmap = {};

        # loop over all vertices and define the polytope obtained by removing one vertex
        foreach my $pv (0..$p_n_verts-1) {
            my $pred = $plp->rows ? new Polytope(POINTS=>$p->VERTICES->minor(~range($pv,$pv),All)/$plp) : new Polytope(POINTS=>$p->VERTICES->minor(~range($pv,$pv),All));
            next if $pred->DIM < $dim;
            my $ret = find_in_centers($pred, $centers, $center_automorphisms);
            if ( $ret->{trafo} != unit_matrix($dim+1) && $p->VERTICES->[$pv]*$ret->{trafo} == $p->VERTICES->[$pv] ) {
                print "skipped one\n";
            }
            $pmap->{$pv} = $ret;
        }
        push @$polys, $pmap;
    }
    return ($polys, $centers, $center_automorphisms);
}

sub check_center_lattice_equivalence {
    my ($p,$q) = @_;

    my $ret;

    if ( $p->N_VERTICES != $q->N_VERTICES ||
         $p->F_VECTOR != $q->F_VECTOR ||
         $p->LATTICE_VOLUME != $q->LATTICE_VOLUME ||
         $p->N_INTERIOR_LATTICE_POINTS != $q->N_INTERIOR_LATTICE_POINTS
         ||
        !isomorphic($p,$q)
        ) {
             $ret->{equiv} = false;
             return $ret;
         }

    my $r = $p->CONE_DIM;
    if ($p->DIM != $q->DIM ) {
        croak "polytopes of different dimension\n";
    }
    my $m = $p->VERTICES->minor(basis_rows($p->VERTICES),All);
    my $sp = new Set<Vector<Rational>>(@{$p->VERTICES});

    foreach my $s (@{all_subsets_of_k(range(0,$q->N_VERTICES-1),$r)}) {
        foreach my $perm (@{all_permutations($r)}) {
            my $n = $q->VERTICES->minor($s,All)->minor($perm,All);
            next if det($n) == 0;
            my $t = inv($n) * $m;
            if ( $t == primitive($t) && abs(det($t)) == 1 ) {
                my $sq = new Set<Vector<Rational>>(@{$q->VERTICES*$t});
                if ( incl($sq,$sp) == 0 ) {
                    $ret->{trafo} = inv($t);
                    $ret->{equiv} = true;
                    return $ret;
                }
            }
        }
    }
    $ret->{equiv} = false;
    return $ret;
}

sub check_lattice_equivalence {
    my ($p,$q) = @_;

#    my $pe = lattice_edge_lengths($p);
#    my $qe = lattice_edge_lengths($q);
#    if ( $pe != $qe ) {
#      print "1";
#      return false;
#    }

    print ".";
    my $r = $p->DIM+1;
    my $m = $p->VERTICES->minor(basis_rows($p->VERTICES),All);
    my $sp = new Set<Vector<Rational>>(@{$p->VERTICES});

    foreach my $s (@{all_subsets_of_k(range(0,$q->N_VERTICES-1),$r)}) {
        foreach my $perm (@{all_permutations($r)}) {
            my $n = $q->VERTICES->minor($s,All)->minor($perm,All);
            next if det($n) == 0;
            my $t = inv($n) * $m;
            if ( $t == primitive($t) && abs(det($t)) == 1 ) {
                my $sq = new Set<Vector<Rational>>(@{$q->VERTICES*$t});
                if ( incl($sq,$sp) == 0 ) {
                    return true;
                }
            }
        }
    }
    return false;
}

sub add_to_list {

    my ($listmap, $list, $lel_list, $gcd_list, $p, $v) = @_;
    my $lel = undef;
    if ( defined($listmap->{$v})) {
#        print "checking ", scalar(@{$listmap->{$v}}), " possible equivalences\n";
        return false if $gcd_list->{$v} == 1;
        foreach my $i (reverse @{$listmap->{$v}}) {
            if ( $p->N_INTERIOR_LATTICE_POINTS == $list->[$i]->N_INTERIOR_LATTICE_POINTS ) {
              if ( !defined($lel) ) {
                $lel = lattice_edge_lengths($p);
              }
              if ( $lel_list->[$i] == $lel ) {
                return false if check_lattice_equivalence($p,$list->[$i]);
              }
            }
        }
    }

    if ( !defined($lel) ) {
      $lel = lattice_edge_lengths($p);
    }
    push @$list, $p;
    push @$lel_list, $lel;
    $gcd_list->{$v} = gcd($v);
    $listmap->{$v} += scalar(@$list)-1;
    return true;
}

sub check_for_common_face {

    my ($p,$i,$j) = @_;
    foreach (@{$p->VERTICES_IN_FACETS}) {
        if ( $_->contains($i) && $_->contains($j) ) {
            return true;
        }
    }
    return false;
}


sub merge {

    my ($a_in,$a_red,$centers,$center_automorphisms) = @_;
    my $total_start = DateTime->now();
    my $s = $a_in->[0]->N_LATTICE_POINTS;
    my $listmap = new Map<Vector<Integer>,Set<Int>>;
    my $list = [];
    my $lel_list = [];
    my $gcd_list = new Map<Vector<Integer>,Integer>;
    my $times = [];
    my ($start, $end);
    $PolyDB::default::db_section_name="Polytopes.Lattice";
    $PolyDB::default::db_collection_name="FewLatticePoints";
    my $asize = scalar(@$a_in)-1;
    if ( $a_red == undef ) {
        print "preparing centers ...\n";
        ($a_red, $centers, $center_automorphisms) = prepare_centers($a_in);
        print " ---done\n";
    }
    my $c=cube($a_in->[0]->DIM);
    my $v=2*ones_vector($a_in->[0]->DIM+1);
    $v->[0]--;
    my $d = coning($c,$v);
    my $sched=$d->get_schedule("VERTICES", "CONE_DIM", "FACETS", "N_LATTICE_POINTS", "N_INTERIOR_LATTICE_POINTS", "F_VECTOR");

    my $listsize = 0;
    foreach my $i (0..$asize) {
        print $i, " -- ", scalar @$list, "\n";
        my $stat = {};
        foreach my $id (keys %{$listmap}) {
            $stat->{scalar(@{$listmap->{$id}})}++;
        }
        foreach my $id (sort { $a <=> $b } keys %$stat) {
            print "(", $id, ": ", $stat->{$id}, ")  ";
        }
        print "\n";
        print join " ", @$times;
        print "\n";
        foreach my $ik (keys %{$a_red->[$i]}) {
            my $pt = transform($a_in->[$i],dense($a_red->[$i]->{$ik}->{trafo}));
            foreach my $j ($i..$asize) {
                foreach my $jk (keys %{$a_red->[$j]}) {
                    next if $a_red->[$i]->{$ik}->{index} != $a_red->[$j]->{$jk}->{index};

                    my $v2 = $a_in->[$j]->VERTICES->[$jk]*$a_red->[$j]->{$jk}->{trafo};

                    foreach my $aq ( @{$center_automorphisms->[$a_red->[$j]->{$jk}->{index}]} ) {
                        $start = time();
                        my $r = coning($pt,$v2*$aq);
                        $sched->apply($r);
                        $end = time();
                        $times->[0] += $end-$start;

                        next if $r->N_LATTICE_POINTS != $s+1;

                        $start = time();
#                        my $vv = lex_min_volume_vector($r);
                        my $vv = new Vector<Integer>(volume_vector($r->VERTICES));
                        $end = time();
                        $times->[1] += $end-$start;

                        $start=time();
                        $listsize += 1 if add_to_list($listmap,$list,$lel_list,$gcd_list,$r,$vv);
                        $end = time();
                        $times->[2] += $end-$start;
                    }
                }
            }
        }
    }
    print "\n", $total_start, "\n", DateTime->now(), "\n";
    return $list;
}


sub is_quasi_minimal {

    my $is_qm = 2;
    my $rem_points = new Set<Int>;
    my $p = shift;
    my $p_vertex_set = new Set<Vector<Rational>>(@{$p->VERTICES});
    my $plp = new Matrix<Rational>(0,$p->DIM+1);
    foreach (@{$p->LATTICE_POINTS}) {
        if ( !$p_vertex_set->contains($_)) {
            $plp /= $_;
        }
    }
    foreach my $i (0..$p->N_VERTICES-1) {
        my $q = new Polytope(POINTS=>$p->VERTICES->minor(~[$i],All)/$plp);
        if ( $q->LATTICE_WIDTH > 1 ) {
            $is_qm--;
            $rem_points += $i;
        }
        if ( $is_qm <= 0 ) { return $rem_points; }
    }
    return true;
}

sub delete_vertex {

    my ($p,$i) = @_;

    my $p_vertex_set = new Set<Vector<Rational>>(@{$p->VERTICES});
    my $plp = new Matrix<Rational>(0,$p->DIM+1);
    foreach (@{$p->LATTICE_POINTS}) {
        if ( !$p_vertex_set->contains($_)) {
            $plp /= $_;
        }
    }

    return new Polytope(POINTS=>$p->VERTICES->minor(~[$i],All)/$plp);
}

sub delete_two_vertices {

    my ($p,$i,$j) = @_;

    my $p_vertex_set = new Set<Vector<Rational>>(@{$p->VERTICES});
    my $plp = new Matrix<Rational>(0,$p->DIM+1);
    foreach (@{$p->LATTICE_POINTS}) {
        if ( !$p_vertex_set->contains($_)) {
            $plp /= $_;
        }
    }

    return new Polytope(POINTS=>$p->VERTICES->minor(~[$i,$j],All)/$plp);
}


sub initial_volume_vector {

    my ($m,$subsets) = @_;

    my $v = new Vector<Integer>(new Int(binomial($m->rows,$m->cols)));
    my $i = 0;
    foreach my $s (@$subsets) {
        $v->[$i++] = det($m->minor($s,All));
    }
    return $v;
}

sub smaller_volume_vector {

    my ($m,$v, $perm, $subsets, $dets) = @_;
    my $i = 0;
    my $undecided = true;
    my $sign = 1;
    my $first = true;
    foreach my $s (@$subsets) {
        my @ps = @$perm[@$s];
        my $psign = 1;
        foreach (0..$m->cols-2) {
            foreach my $id (0..$m->cols-2) {
                if ( $ps[$id] > $ps[$id+1] ) {
                    my $temp = $ps[$id];
                    $ps[$id] = $ps[$id+1];
                    $ps[$id+1] = $temp;
                    $psign *= -1;
                }
            }
        }
        my $d = $psign * $dets->{\@ps};
        if ( $first && $v->[0] >= -$d ) {
            $sign = -1;
        }
        $first = false;
        if ( $v->[$i] == $sign*$d ) {
            $i++;
            next;
        }
        if ( $undecided )  {
            if ( $sign*$d > $v->[$i] ) {
                return;
            } else {
                $undecided = false;
            }
        }
        $v->[$i++] = $sign*$d;
    }
}

sub lex_comparison {
    my ($v,$w) = @_;
    foreach my $i (0..$v->dim-1) {
        if ( $v->[$i] > $w->[$i] ) {
            return 1;
        } elsif( $v->[$i] < $w->[$i] ) {
            return -1;
        }
    }
    return 0;
}

sub lex_min_volume_vector {

    my $p = shift;
    my $m = $p->VERTICES;
    my $mr=$m->rows;
    my $all_perms = all_permutations($mr);
    my $subsets = new Array<Array<Int>>(all_subsets_of_k(range(0,$mr-1),$m->cols));
    my $dets = new Map<Array<Int>,Int>;
    foreach my $s (@$subsets) {
        $dets->{$s} = det($m->minor($s,All));
    }
    my $v = initial_volume_vector($m,$subsets);
    foreach my $perm  (@$all_perms) {
        smaller_volume_vector($m, $v, $perm, $subsets, $dets);
    }
    return $v;
}

sub transform_a_red {

    my $d = shift; # Array<Map<Int,Pair<Int,Matrix<Rational>>>>
    my $p = [];

    foreach my $e (@$d) { # Map<Int,Pair,Int,Matrix<Rational>>>
        my $s = ();
        foreach my $k (keys %$e ) { 
            my $r = ();
            $r->{index} = $e->{$k}->first;
            $r->{trafo} = $e->{$k}->second;
            $s->{$k} = $r;
        }
        push @$p, $s;
    }
    return $p;
}