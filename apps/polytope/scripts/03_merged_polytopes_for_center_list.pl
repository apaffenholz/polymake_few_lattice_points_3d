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

use File::Basename;

script("polytopes_from_center.pl");

my ($in_name, $basepath) = @ARGV;

$Verbose::credits=0;
my($filename, $dir, $suffix) = fileparse($in_name);
if ( !defined($basepath) ) {
    $basepath = $dir;
}
$in_name = $filename;

my $centers = load_data($basepath.'/'.$in_name);

my $m= new Matrix<Rational>(@{$centers->[0]->{'vertices'}->[0]});
my $p=new Polytope(POINTS=>$m);
my $size = new Int($p->N_LATTICE_POINTS+1);

my ($out_name) = $in_name =~ m/([0-9]+.pdata)/i;
$out_name = "partial_merged_list.".$size."v.".$out_name;

my $total_centers = scalar(@$centers);

my $list = [];
my $list_hash = {};
foreach my $c (@$centers) {
    polytopes_from_center($c,$size,$list,$list_hash);
}
save_data($list,$basepath.'/'.$out_name);