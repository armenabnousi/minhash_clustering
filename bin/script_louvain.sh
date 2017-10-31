#!/bin/bash
#PBS -k o
#PBS -l nodes=4:intel:ppn=16,walltime=6:00:00,mem=550GB
##PBS -l nodes=3:smp:ppn=32+1:ppn=30,mem=360GB
#PBS -m abe
#PBS -N 64_zc_shlouv_dbld
#PBS -M armen.abnousi@wsu.edu
#PBS -j oe

module load compilers/gcc/6.2.0
module load openmpi/1.10.1-gcc

cd ${PBS_O_WORKDIR}
#rm outputs/new_l200_k6sh2m2cc-1h*

nproc=64
cat $PBS_NODEFILE >> "node_list_for_this_run"$nproc
setnumber_array=$(seq 70 1 70)
setname="large"
groundset="bact11"
set_initial=${setname:0:1}

score_threshold=0.3
length_threshold=0.8
similarity_function="similar"
#similarity_function="homologous"
variable="_num_hash_fs"
#variable="_shingling_iters"
#variable="_cc_iters"
difference=40
num_hash_fs_array=$(seq 140 1 140)
start_hash_number=41
shingling_iters_array=$(seq 2 1 2)
cc_iters_array=$(seq -1 1 -1)
min_cluster_size=10
collect_clusters=0
rename=0
print_cluster_stats_limit=-1
final_cluster_size_limit=100 #this is extra for now
#cc_algorithm="union_find"
#cc_algorithm="label_propagation"
#cc_algorithm="hash_to_min"
cc_algorithm="nocc"
com_detection="usc_louvain"
#com_detection="grappolo"
#com_detection="no_com_detection"
metis_prefix="/data/hpcbio/users/armen/shingling/proc"$nproc"/test_"
fvalue_filename="/data/hpcbio/users/armen/shingling/proc"$nproc"/fvalue_temp"
fvalue_threshold=1.1

num_process_array=($nproc)
page_size_array=(10000)
#page_size_array=(4096)
num_maps_array=($nproc)
separation_delta_array=(30000)
kmer_length_array=$(seq 6 1 6)
shingle_size_array=$(seq 2 1 2)

declare test_variable="test"$variable
declare base_variable="base"$variable
#cat $PBS_NODEFILE #> /home/armen.abnousi/boostclust/nodelist.log

get_seq_count () {
        if [ "$1" == "large70" ]; then count=94599
	elif [ "$1" == "large200" ]; then count=217936
        elif [ "$1" == "small2" ]; then count=1938
	elif [ "$1" == "small3" ]; then count=2056
	elif [ "$1" == "small4" ]; then count=2139
	elif [ "$1" == "small5" ]; then count=2565
	elif [ "$1" == "small6" ]; then count=1479
	elif [ "$1" == "small7" ]; then count=1494
	elif [ "$1" == "small8" ]; then count=2037
	elif [ "$1" == "small9" ]; then count=808
	elif [ "$1" == "small10" ]; then count=1424
	elif [ "$1" == "small11" ]; then count=1542
	elif [ "$1" == "NC_g11" ]; then count=13571
	elif [ "$1" == "test_small2" ]; then count=200
	elif [ "$1" == "sample5" ]; then count=5186
        fi
        echo $count
}

#printf "base\tprocs\tmaps\tpage\tdelta\tseqs\tk\tsh\tset\ttest_h\tbase_h\tt_m\tb_m\tt_cc\tb_cc\tdiff\tcsl\tt_lrg\tnmi_lfk\tnmi_max\n" >> $nmi_output
for num_process in ${num_process_array[@]}
do
	for page_size in ${page_size_array[@]}
	do
		for num_maps in ${num_maps_array[@]}
		do
			for test_cc_iters in ${cc_iters_array[@]}
			do
				for test_shingling_iters in ${shingling_iters_array[@]}
				do
					for sep_delta in ${separation_delta_array[@]}
					do
						for kmer_length in ${kmer_length_array[@]}
						do
							for setnumber in ${setnumber_array[@]}
							do
								for test_num_hash_fs in ${num_hash_fs_array[@]}
								do
									for shingle_size in ${shingle_size_array[@]}
									do
										sth=`echo "$score_threshold * 10"|bc`
										lth=`echo "$length_threshold * 10"|bc`
										sth_rounded=`echo "$sth" | awk '{printf "%d\n", $1}'`
										lth_rounded=`echo "$lth" | awk '{printf "%d\n", $1}'`
										finalized_dir="finalized_clusters/finalized_clusters_"$similarity_function"_"$set_initial$setnumber"_sth"$sth_rounded"_lth"$lth_rounded"_proc"$nproc"_nocc_wf1score"
										mkdir $finalized_dir
										nmi_output="/home/armen.abnousi/shingle/zc_shingle_looped_"$set_initial$setnumber"m"$test_shingling_iters"cc"$test_cc_iters"h1to500_diff"$difference"_cslim"$min_cluster_size"_sth"$sth"_lth"$lth"_"$similarity_function"_proc"$nproc"_nocc_wf1score"
										seq_count=$(get_seq_count $setname$setnumber)
										base_num_hash_fs=$test_num_hash_fs
										base_shingling_iters=$test_shingling_iters
										base_cc_iters=$test_cc_iters
										declare $base_variable=$(expr ${!test_variable} - $difference)
										#printf "shin\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t"  $num_process $num_maps $page_size $sep_delta $seq_count $kmer_length $shingle_size $setnumber $test_num_hash_fs $base_num_hash_fs $test_shingling_iters $base_shingling_iters $test_cc_iters $base_cc_iters $difference $min_cluster_size >> $nmi_output
#										printf "t_h\tb_h\ttkst\ttvst\ttcct\ttct\tbkst\tbvst\tbcct\tbct\tnmit\tmert\ttotlt\tlrgst\tcrem\tsrem\tsstay\ts_th\tl_th\tnscomp\tnmi_lfk\tpre_nmi\tnmi_max\tkmer_h\tkmercollate\tgraph\tcustomhash\tredun\n" >> $nmi_output
										#file_name="/home/armen.abnousi/shingle/datasets/s"$setnumber"m6regions.fasta"
										#file_name="/home/armen.abnousi/shingle/datasets/s"$setnumber"_pfam_regions"
										#if [ "$setname" == "small" ];
										#	then file_name="datasets/"$set_initial$setnumber"m6regions_by_pfam29_smoothw10.fasta_numbered"
										#else
											#file_name="datasets/"$setname$setnumber"/"$setname$setnumber"doms_by_nadda_trained_by_"$groundset"pfam29_numbered"
											#file_name="datasets/sample5/sample5doms_extracted_from_large70doms_nadda_numbered"
											#file_name="datasets/"$setname$setnumber"/"$setname$setnumber"doms_by_nadda_6mer_trained_by_bact11th50cdhit04_minlen50_numbered"
											#file_name="datasets/NC_g11/NC_g11doms_by_nadda_6mer_trained_by_bact11th50cdhit04_minlen5"
											file_name="datasets/"$setname$setnumber"/"$setname$setnumber"_domains_extracted_from_pfam_numbered"
											#file_name="datasets/test_small2/test_small4_domains_extracted_from_pfam_numbered"
											#file_name="datasets/tigerfam61/tigerfam_all_seeds_numbered.fasta"
										#fi
										clusters_output_filename="outputs/not_exist_new_"$set_initial$setnumber
										#echo "starting"
										setid=$set_initial$setnumber
										#setid="t61"
										#setid="sample5"
										echo $setid
										#echo $num_maps $page_size $sep_delta $file_name $clusters_output_filename $seq_count $kmer_length $shingle_size $test_num_hash_fs $test_shingling_iters $test_cc_iters $nmi_output $min_cluster_size $rename $collect_clusters $print_cluster_stats_limit $cc_algorithm $score_threshold $length_threshold $similarity_function $finalized_dir $com_detection $setid $start_hash_number $metis_prefix $fvalue_filename $fvalue_threshold $difference
										mpiexec -n $num_process ./execshingle_louvain_doubled $num_maps $page_size $sep_delta $file_name $clusters_output_filename $seq_count $kmer_length $shingle_size $test_num_hash_fs $test_shingling_iters $test_cc_iters $nmi_output $min_cluster_size $rename $collect_clusters $print_cluster_stats_limit $cc_algorithm $score_threshold $length_threshold $similarity_function $finalized_dir $com_detection $setid $start_hash_number $metis_prefix $fvalue_filename $fvalue_threshold $difference #$final_cluster_size_limit
										#echo "done"
#										echo $(which mpiexec)
#										mpirun -n valgrind --tool=memcheck --log-file=log_valgrind --leak-check=full --show-reachable=yes ./execshingle_louvain 1000 2048 30000 /home/armen.abnousi/shingle/datasets/small2/small2doms_by_nadda_trained_by_bact11pfam29_numbered outputs/new_s2 1938 6 2 3 2 2 2 -1 -1 zc_shingle_s2m2cc-1h1to50_diff1_cslim 10 1 0 -1
#										num_process=1
#										mpiexec -n $num_process valgrind --tool=memcheck --leak-check=full --show-reachable=yes --log-file="log_valgrind_looped_small2" ./execshingle_louvain $num_maps $page_size $sep_delta $file_name $clusters_output_filename $seq_count $kmer_length $shingle_size $test_num_hash_fs $test_shingling_iters $test_cc_iters $nmi_output $min_cluster_size $rename $collect_clusters $print_cluster_stats_limit $cc_algorithm $score_threshold $length_threshold $similarity_function $finalized_dir $com_detection #$final_cluster_size_limit
#										mpiexec -n $num_process valgrind --tool=massif --stacks=yes --time-unit=ms --log-file="log_valgrind_looped_small2" ./execshingle_louvain $num_maps $page_size $sep_delta $file_name $clusters_output_filename $seq_count $kmer_length $shingle_size $test_num_hash_fs $test_shingling_iters $test_cc_iters $nmi_output $min_cluster_size $rename $collect_clusters $print_cluster_stats_limit $cc_algorithm $score_threshold $length_threshold $similarity_function $finalized_dir $com_detection $setid #$final_cluster_size_limit
									done
								done
							done
						done
					done
				done
			done
		done
	done
done
rm "zc_polo*"
#TO DEBUG RUN "gdb --args ./execshingle $num_maps $page_size $kmer_length $sep_delta $file_name $num_hash_fs $shingle_size
