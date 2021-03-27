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

my ($in_name,$basepath) = @ARGV;

use File::Basename;

$Verbose::credits=0;
my($filename, $dir, $suffix) = fileparse($in_name);
if ( !defined($basepath) ) {
    $basepath = $dir;
}
$in_name = $filename;

my $mat_strings = load_data($basepath.'/'.$in_name);

my $ms_size = $mat_strings->[0];
my @name_mat_size = split /---/, $ms_size;
my @entries_size = split /,/, $name_mat_size[1];
my $marr_size = []; 
foreach my $s (@entries_size) { 
    my @r = split / /, $s; 
    push @$marr_size, \@r; 
} 
my $m_size = new Matrix<Rational>($marr_size);
my $p_size = new Polytope(VERTICES=>((ones_vector<Rational>($m_size->rows()))|$m_size));
my $size = new Int($p_size->N_LATTICE_POINTS);

my ($out_name) = $in_name =~ m/([0-9]+.pdata)/i;
$out_name = $size."v.".$out_name;

my $polytope_list = [];
foreach my $ms (@$mat_strings) {
    my @name_mat = split /---/, $ms;
    my $name = $name_mat[0];
    my @entries = split /,/, $name_mat[1];
    my $marr = []; 
    foreach my $s (@entries) { 
        my @r = split / /, $s; 
        push @$marr, \@r; 
    } 
    my $m = new Matrix<Rational>($marr);
    my $p = new Polytope(VERTICES=>((ones_vector<Rational>($m->rows()))|$m), FEASIBLE=>true);
    $p->name = $name;
    push @$polytope_list, $p;
}
save_data($polytope_list,$basepath.'/'.$out_name);

1;
