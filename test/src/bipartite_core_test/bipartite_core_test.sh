#!/bin/bash
PATH='/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/bin'
export PATH

PREFIX=$HOME

#thread_number1=16

#small_data1=$PREFIX'/wbai/bipartite_graph/small_data/non_temporal_graph/'
#./bipartite_core_shuffle_non_temporal_bipartite_edges $small_data1 $thread_number1
##shellcheck disable=SC2045
#for file in $(ls $small_data1)
#do
#    if [ -f $small_data1$file ]
#    then
#        echo $file
#  	    ./bipartite_core_decomposition_test $small_data1 $file $thread_number1
#  	    ./bipartite_core_insertion_test $small_data1 $file $thread_number1
#        ./bipartite_core_removal_test $small_data1 $file $thread_number1
#        ./bipartite_core_insertion_removal_test $small_data1 $file $thread_number1
#    fi
#done
#
#small_data2=$PREFIX'/wbai/bipartite_graph/small_data/unique_temporal_graph/'
##shellcheck disable=SC2045
#for file in $(ls $small_data2)
#do
#    if [ -f $small_data2$file ]
#    then
#        echo $file
#  	    ./bipartite_core_decomposition_test $small_data2 $file $thread_number1
#  	    ./bipartite_core_insertion_test $small_data2 $file $thread_number1
#        ./bipartite_core_removal_test $small_data2 $file $thread_number1
#        ./bipartite_core_insertion_removal_test $small_data2 $file $thread_number1
#    fi
#done

thread_number2=16

#large_data1=$PREFIX'/wbai/bipartite_graph/large_data/non_temporal_graph/'
#./bipartite_core_shuffle_non_temporal_bipartite_edges $large_data1 $thread_number2
## shellcheck disable=SC2045
#for file in $(ls $large_data1)
#do
#   if [ -f $large_data1$file ]
#   then
#       echo $file
#       ./bipartite_core_decomposition_test $large_data1 $file $thread_number2
#       ./bipartite_core_insertion_compare_test $large_data1 $file $thread_number2
#       ./bipartite_core_removal_compare_test $large_data1 $file $thread_number2
#        ./bipartite_core_insertion_size_test $large_data1 $file $thread_number2
#        ./bipartite_core_removal_size_test $large_data1 $file $thread_number2
#	     ./bipartite_core_insertion_thread_test $large_data1 $file $thread_number2
#       ./bipartite_core_removal_thread_test $large_data1 $file $thread_number2
#   fi
#done


large_data2=$PREFIX'/wbai/bipartite_graph/large_data/unique_temporal_graph/'
# shellcheck disable=SC2045
for file in $(ls $large_data2)
do
   if [ -f $large_data2$file ]
   then
       echo $file
#       ./bipartite_core_decomposition_test $large_data2 $file $thread_number2
#       ./bipartite_core_insertion_compare_test $large_data2 $file $thread_number2
#       ./bipartite_core_removal_compare_test $large_data2 $file $thread_number2
#        ./bipartite_core_insertion_size_test $large_data2 $file $thread_number2
#        ./bipartite_core_removal_size_test $large_data2 $file $thread_number2
	     ./bipartite_core_insertion_thread_test $large_data2 $file $thread_number2
       ./bipartite_core_removal_thread_test $large_data2 $file $thread_number2
   fi
done


