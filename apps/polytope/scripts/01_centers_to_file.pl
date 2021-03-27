
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

script("polytope_centers.pl");

my ($in_name, $targetpath, $basepath) = @ARGV;

$Verbose::credits=0;
my($filename, $dir, $suffix) = fileparse($in_name);
if ( !defined($basepath) ) {
    $basepath = $dir;
}
if ( !defined($targetpath) ) {
    $targetpath = $dir;
}
$in_name = $filename;
my $out_name = "partial_centers.".$in_name;

my $list = load_data($dir.'/'.$in_name);
my $centers = centers_list_of_polytopes($list);
save_data($centers,$targetpath.'/'.$out_name);

1;