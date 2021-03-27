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

my ($basepath,$prefix, $n) = @ARGV;
my @files = glob( $basepath . '/'.$prefix.'*' );
    
$Verbose::credits=0;
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

            $centers->[$index]->{'center_string'} = $center_string;
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

my $total_number = scalar(@$centers);
print "found ".$total_number." centers\n";
print "now splitting into files\n";
my $out_name = $prefix =~ s/partial_centers/centers/r;
$out_name =~ s/.pdata//;

my $t = $total_number;
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

my $partial_list = [];
my $total_count = 0;
$count = 0;
my $file = 0;
foreach (@$centers) {
    push @$partial_list, $_;
    $count++;
    $total_count++;
    if ( $count % $n == 0 ) {
        my $num = sprintf "%0*d", $digits, $file;
        my $name = $basepath.'/'.$out_name.".".$num.".pdata";
        save_data($partial_list,$name);
        $file++;
        $partial_list = [];
        $count = 0;
    }
}
if ( scalar(@$partial_list) > 0 ) {
    my $num = sprintf "%0*d", $digits, $file;
    my $name = $basepath.'/'.$out_name.".".$num.".pdata";
    save_data($partial_list,$name);
    $file++;
}

print "wrote ".$total_count." centers into ".$file." files\n";
