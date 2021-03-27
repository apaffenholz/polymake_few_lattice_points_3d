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

script("generate_polytopes.pl");

my ($basepath,$prefix,$n) = @ARGV;
my @files = glob( $basepath . '/'.$prefix.'*' );

$Verbose::credits=0;
my $hlist = ();
my $full_count = 0;
my $added_count;
foreach my $file (@files) {
    my $polys = load_data($file);
    foreach my $p (@$polys) {
        $full_count++;
        if ( !defined($hlist->{$p}) ) {
            $added_count++;
        }
        $hlist->{$p} = 1;
    } 
}

printf "read ".$full_count." merged polytopes and found ".$added_count." different ones\n";

# get size
my $ms_size = (keys %$hlist)[0];
my @entries_size = split /,/, $ms_size;
my $marr_size = []; 
foreach my $s (@entries_size) { 
    my @r = split / /, $s; 
    push @$marr_size, \@r; 
} 
my $m_size = new Matrix<Rational>($marr_size);
my $p_size = new Polytope(VERTICES=>(ones_vector<Rational>($m_size->rows()))|$m_size);
my $size = new Int($p_size->N_LATTICE_POINTS);

my $other_polytopes = spiked_polytopes($size);
my $boxed = boxed_quasi_minimal_polytopes($size);
push @$other_polytopes, @$boxed;

$added_count += scalar(@$other_polytopes);

print "total number of polytopes of size ".$size." is ".$added_count."\n";

# write polytopes to files
my $t = $added_count;
my $digits = 1;
while ( $t > 10 ) {
    $digits++;
    $t = $t/10;
}
$t = $n;
while ( $t > 10 ) {
    $digits--;
    $t = $t/10;
}

my $polytope_list = [];
my $count = 0;
my $file = 0;
foreach my $ms (keys %$hlist) {
    my $name = sprintf "merged.s%02d.%0*d", $size, $size-4, $count;
    my $mat_name_string = $name.'---'.$ms;
    push @$polytope_list, $mat_name_string;
    $count++;
    if ( $count % $n == 0 ) {
        my $num = sprintf "%0*d", $digits, $file;
        my $name = $basepath."/string.".$size."v.".$num.".pdata";
        save_data($polytope_list,$name);
        $file++;
        $polytope_list = [];
        $count = 0;
    }
}
foreach my $p (@$other_polytopes) {
    my $V = dehomogenize($p->VERTICES);
    my $name = $p->name;
    my $ms = join ",", @$V;
    my $mat_name_string = $name.'---'.$ms;
    push @$polytope_list, $mat_name_string;
    $count++;
    if ( $count % $n == 0 ) {
        my $num = sprintf "%0*d", $digits, $file;
        my $name = $basepath."/string.".$size."v.".$num.".pdata";
        save_data($polytope_list,$name);
        $file++;
        $polytope_list = [];
        $count = 0;
    }
}


if ( $count > 0 ) {
    my $num = sprintf "%0*d", $digits, $file;
    my $name = $basepath."/string.".$size."v.".$num.".pdata";
    save_data($polytope_list,$name);
}
