# Copyright (c) 2003-2017 Broad Institute, Inc., Massachusetts Institute of Technology, and Regents of the University of California.  All rights reserved.
#!/bin/sh
execDir=`dirname $0`

zip1=$1
zip2=$2

base1=`basename $1`
base2=`basename $2`
bare1=${base1%%.zip}
bare2=${base2%%.zip}
diffDir1=`mktemp -d $bare1.XXXXXX`
diffDir2=`mktemp -d $bare2.XXXXXX`

unzip -q $zip1 -d $diffDir1
unzip -q $zip2 -d $diffDir2

# Diff only selected files out of the ZIP
diff -i --strip-trailing-cr -q $diffDir1/edb/$bare1.rnk $diffDir2/edb/$bare2.rnk
status=$?
diff -i --strip-trailing-cr -q $diffDir1/edb/results.edb $diffDir2/edb/results.edb
status=$(( $? + status ))

# The generated GMT files sometimes have a different order when there are multiple GMT inputs.
# This should probably be changed in the long run, but for now we just compare the sorted outputs
# This does *NOT* happen with "regular" GSEA so this script differs at this point.
sort $diffDir1/edb/gene_sets.gmt > $diffDir1/edb/sorted.matrix
sort $diffDir2/edb/gene_sets.gmt > $diffDir2/edb/sorted.matrix

diff -i --strip-trailing-cr -q $diffDir1/edb/sorted.matrix $diffDir2/edb/sorted.matrix
status=$(( $? + status ))

rm -rf $diffDir1 $diffDir2
exit $status
