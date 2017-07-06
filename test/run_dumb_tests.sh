#!/usr/bin/zsh
#
# this script runs a lot of tests and produces a lot of output which needs to
# be filtered later on to extract the timing information
#

num_trials=$1

DIMENSIONS=2
FILE_NAME=dummesh

n=1
while [[ $n -le $num_trials ]]
do
    o=2
    while [[ $o -le 5 ]]
    do
        i=1
        while [[ $i -le 9 ]]
        do
            rd=$((RANDOM*RANDOM % 1e6)) # random number up to a million
            rd=$(( (rd/1e6 - 0.5) * 0.2)) # turn it into random number between -0.2 and 0.2
            num_samples=$(((1+$rd)*$i*1e$o))
            num_samples=${num_samples%.*}
            echo "run number "$n,$o,$i" with "$num_samples"points"
            rbox D$DIMENSIONS $num_samples > $FILE_NAME
            cat $FILE_NAME | qhull d i >> $FILE_NAME
            validmesh=1
            ./qhull2ply dummesh || validmesh=0
            if (( $validmesh ))
                # qhull2ply fails if qhull gave us a mesh that we cannot process,
                # so we try again
                then echo "test run!"
                ./yamltest ./dumbexample.yaml >> dumb_testlog
                echo "test done" >> dumb_testlog
                i=$((i+1))
            else echo "qhull tricked us! run again!"
            fi
        done
        o=$((o+1))
    done
    n=$((n+1))
done

#
# here we filter the output produced just above
#

echo "faces , edges , vertices , Eigen , coeff_flow" > results.csv

#
# TEST END
#
###################################################################
#
# PROCESSING THE OUTPUT
#

# filter the lines that contain the data we want to keep
# "test done" --> signals the end of a test
# number of --> vertices, edges, faces
# bounding chain --> results from eigen and coeff flow
cat dumb_testlog | \
    grep '\(^test done$\|number of\|bounding chain\)' >> results.csv

# remove the parts of the lines that aren't the data (i.e. the useless text)
sed 's/\(.*number of [a-z]\+: \|test done\|\(calculated.\+\)\@=[0-9]\+\.\|calculated.*cycles \| seconds\)//g' -i results.csv

# put the data from each test into a single comma separated line
sed ':a;/[0-9]$/{N;s/\n/ , /;ba}' -i results.csv

# delete the line endings
sed 's/ *, *$//' -i results.csv

# deleting blank lines and lines where tests didn't complete due to a gsimp::out_of_context exception
# (meaning that there was probably some error while making the complex)
sed '/\(^[0-9\. ]\+\(, [0-9\. ]\+\)\{3\}$\|^$\)/d' -i results.csv

