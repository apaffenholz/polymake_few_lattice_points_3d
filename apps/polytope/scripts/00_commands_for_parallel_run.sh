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


#!/bin/bash -x

basepath=$1
scriptpath=$2
nv=$3;
split=$4
new_split=$5

nv_prev=$((nv-1));

mkdir -p ${basepath}/${nv}v

find ${basepath}/${nv_prev}v -regextype sed -regex  ".*/${nv_prev}v.[0-9].*pdata" -exec sem  --id centers -j 50 perl/polymake --script ${scriptpath}/01_centers_to_file.pl {} ${basepath}/${nv}v \;
sem --id centers --wait
echo "partial centers computed"

perl/polymake --script ${scriptpath}/02_merge_centers_from_file.pl ${basepath}/${nv}v partial_centers.${nv_prev}v ${split}
find ${basepath}/${nv}v -regextype sed -regex  ".*/centers.${nv_prev}v.[0-9].*pdata" -exec sem --id merge -j 50 perl/polymake --script ${scriptpath}/03_merged_polytopes_for_center_list.pl {} \;
sem --id merge --wait
echo "partial merged polytopes computed"

perl/polymake --script ${scriptpath}/04_read_matrices_from_files.pl ${basepath}/${nv}v partial_merged_list.${nv}v ${split} ${new_split}
find ${basepath}/${nv}v -regextype sed -regex  ".*/string.*pdata" -exec sem --id convert -j 50 perl/polymake --script ${scriptpath}/05_mat_string_to_poly.pl {} \;
sem --id convert --wait
echo "vertex list converted to polytopes"

find ${basepath}/${nv}v -name "partial_*" -exec rm {} \;
find ${basepath}/${nv}v -name "centers.*" -exec rm {} \;
find ${basepath}/${nv}v -name "string.*"  -exec rm {} \;