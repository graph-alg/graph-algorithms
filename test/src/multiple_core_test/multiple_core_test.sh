#!/bin/bash
PATH=/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/bin
export PATH

thread_number=14

#small_data1=$HOME'/wbai/multiple_graph/processed_data/small_data/unipartite_graph/'
#rm -rf $small_data1/log
##shellcheck disable=SC2045
#for file in $(ls $small_data1)
#do
#     if [ -f $small_data1$file ]
#     then
#         ./multiple_core_index_performance_test $small_data1 $file $thread_number
#  	     ./multiple_core_performance_test $small_data1 $file $thread_number
#     fi
#done
#
#small_data2=$HOME'/wbai/multiple_graph/processed_data/small_data/bipartite_graph/'
#rm -rf $small_data2/log
##shellcheck disable=SC2045
#for file in $(ls $small_data2)
#do
#     if [ -f $small_data2$file ]
#     then
#         ./multiple_core_index_performance_test $small_data2 $file $thread_number
#  	     ./multiple_core_performance_test $small_data2 $file $thread_number
#     fi
#done


large_data2=$HOME'/wbai/multiple_graph/processed_data/large_data/bipartite_graph/'
#rm -rf $large_data2/log
# shellcheck disable=SC2045
for file in $(ls $large_data2)
do
   if [ -f $large_data2$file ]
   then
     #./multiple_core_index_scalability_test $large_data2 $file $thread_number
	   ./multiple_core_compare_scalability_test $large_data2 $file $thread_number
	   ./multiple_core_size_scalability_test $large_data2 $file $thread_number
	   ./multiple_core_thread_scalability_test $large_data2 $file $thread_number
   fi
done


