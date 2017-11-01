#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=2,walltime=2:00:00,mem=15GB
#PBS -j oe

module load compilers/gcc/6.2.0
module load openmpi/1.10.1-gcc

cd ${PBS_O_WORKDIR}
timings="/home/armen.abnousi/shingle/cluster_timings32"

operation="grappolo_to_mrmpi_formatted_file"

nprocs=2

test_num_hash=$1
test_metis_filename=$2
base_metis_filename=$3
fvalue_filename=$4

test_graph_dictionary=$test_metis_filename"_dictionary"
test_grappolo_communities=$test_metis_filename"_clustInfo"
test_output_communities=$test_grappolo_communities"_mrmpi_formated"
base_graph_dictionary=$base_metis_filename"_dictionary"
base_grappolo_communities=$base_metis_filename"_clustInfo"
base_output_communities=$base_grappolo_communities"_mrmpi_formated"
rename=0

min_cluster_size=10
min_print_stats_cluster_size=50 #not used
largest_cluster_count_to_compare=-1

pagesize=4196
seq_count=1000
base_type="shingling"
setname="setname"

echo "calling grappolo for test set $test_metis_filename" >> $timings
SECONDS=0
/home/armen.abnousi/shingle/graph_formater/grappolo-05-2014/driverForGraphClustering -f 5 -o -c $test_metis_filename
/home/armen.abnousi/shingle/graph_formater/graph_formater_driver $operation $test_graph_dictionary $test_grappolo_communities $test_output_communities
echo "armtime: test grappolo time: $SECONDS seconds" >> $timings

echo "calling grappolo for base set $base_metis_filename" >> $timings
SECONDS=0
if [ ! -f $base_output_communities ]; then
	echo "File not found!"
	/home/armen.abnousi/shingle/graph_formater/grappolo-05-2014/driverForGraphClustering -f 5 -o -c $base_metis_filename
	/home/armen.abnousi/shingle/graph_formater/graph_formater_driver $operation $base_graph_dictionary $base_grappolo_communities $base_output_communities
fi
echo "armtime: base grappolo time: $SECONDS seconds" >> $timings

echo "calling nmi compute\n" >> $timings
SECONDS=0
nmap=$nprocs
mpiexec -n $nprocs /home/armen.abnousi/shingle/exec_nmi_shingling_vs_pfam $nmap $pagesize $test_output_communities $base_output_communities $seq_count $min_cluster_size $min_print_stats_cluster_size $largest_cluster_count_to_compare $test_num_hash $base_type $fvalue_filename $setname $rename
echo "armtime: nmi computation time: $SECONDS seconds" >> $timings
