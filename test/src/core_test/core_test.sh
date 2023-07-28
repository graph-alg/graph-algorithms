#!/bin/bash
PATH='/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/bin'
export PATH

thread_number=14

#small_data1=$HOME'/wbai/unipartite_graph/small_data/unique_temporal_graph/'
#rm -rf $small_data1/'log/'
##shellcheck disable=SC2045
#for file in $(ls $small_data1); do
#  if [ -f $small_data1$file ]; then
#    echo $file
#    ./core_decomposition_test $small_data1 $file $thread_number
#    ./traversal_insertion_memory_test $small_data1 $file $thread_number
#    ./order_insertion_memory_test $small_data1 $file $thread_number
#    ./jes_insertion_memory_test $small_data1 $file $thread_number
#    ./quasi_insertion_memory_test $small_data1 $file $thread_number
#    ./parallel_insertion_memory_test $small_data1 $file $thread_number
#    ###
#    ./traversal_removal_memory_test $small_data1 $file $thread_number
#    ./order_removal_memory_test $small_data1 $file $thread_number
#    ./jes_removal_memory_test $small_data1 $file $thread_number
#    ./quasi_removal_memory_test $small_data1 $file $thread_number
#    ./parallel_removal_memory_test $small_data1 $file $thread_number
#    ###
#    ./traversal_insertion_removal_memory_test $small_data1 $file $thread_number
#    ./order_insertion_removal_memory_test $small_data1 $file $thread_number
#    ./jes_insertion_removal_memory_test $small_data1 $file $thread_number
#    ./quasi_insertion_removal_memory_test $small_data1 $file $thread_number
#    ./parallel_insertion_removal_memory_test $small_data1 $file $thread_number
#  fi
#done
##
#small_data2=$HOME'/wbai/unipartite_graph/small_data/non_temporal_graph/'
#rm -rf $small_data2'log/'
#./shuffle_non_temporal_edges $small_data2 $thread_number
##shellcheck disable=SC2045
#for file in $(ls $small_data2); do
#  if [ -f $small_data2$file ]; then
#    echo $file
#    ./core_decomposition_test $small_data2 $file $thread_number
#    ./traversal_insertion_memory_test $small_data2 $file $thread_number
#    ./order_insertion_memory_test $small_data2 $file $thread_number
#    ./jes_insertion_memory_test $small_data2 $file $thread_number
#    ./quasi_insertion_memory_test $small_data2 $file $thread_number
#    ./parallel_insertion_memory_test $small_data2 $file $thread_number
#    ###
#    ./traversal_removal_memory_test $small_data2 $file $thread_number
#    ./order_removal_memory_test $small_data2 $file $thread_number
#    ./jes_removal_memory_test $small_data2 $file $thread_number
#    ./quasi_removal_memory_test $small_data2 $file $thread_number
#    ./parallel_removal_memory_test $small_data2 $file $thread_number
#    ###
#    ./traversal_insertion_removal_memory_test $small_data2 $file $thread_number
#    ./order_insertion_removal_memory_test $small_data2 $file $thread_number
#    ./jes_insertion_removal_memory_test $small_data2 $file $thread_number
#    ./quasi_insertion_removal_memory_test $small_data2 $file $thread_number
#    ./parallel_insertion_removal_memory_test $small_data2 $file $thread_number
#  fi
#done
#
#large_data1=$HOME'/wbai/unipartite_graph/large_data/unique_temporal_graph/'
#rm -rf $large_data1'log/'
## shellcheck disable=SC2045
#for file in $(ls $large_data1); do
#  if [ -f $large_data1$file ]; then
#    echo $file
#    ./core_decomposition_test $large_data1 $file $thread_number
#    ./core_insertion_compare_test $large_data1 $file $thread_number
#    ./core_removal_compare_test $large_data1 $file $thread_number
#    ./core_insertion_removal_compare_test $large_data1 $file $thread_number
#  fi
#done
#
#large_data2=$HOME'/wbai/unipartite_graph/large_data/non_temporal_graph/'
#rm -rf $large_data2'log/'
#./shuffle_non_temporal_edges $large_data2 $thread_number
## shellcheck disable=SC2045
#for file in $(ls $large_data2); do
#  if [ -f $large_data2$file ]; then
#    echo $file
#    ./core_decomposition_test $large_data2 $file $thread_number
#    ./core_insertion_compare_test $large_data2 $file $thread_number
#    ./core_removal_compare_test $large_data2 $file $thread_number
#    ./core_insertion_removal_compare_test $large_data2 $file $thread_number
#  fi
#done
##
#large_data3=$HOME'/wbai/unipartite_graph/large_data/random_graph/'
#./shuffle_non_temporal_edges $large_data3 $thread_number
#rm -rf $large_data3'log/'
## shellcheck disable=SC2045
#for file in $(ls $large_data3); do
#  if [ -f $large_data3$file ]; then
#    echo $file
#    ./core_decomposition_test $large_data3 $file $thread_number
#    ./core_insertion_size_test $large_data3 $file $thread_number
#    ./core_removal_size_test $large_data3 $file $thread_number
#    ./core_insertion_removal_size_test $large_data3 $file $thread_number
#  fi
#done
##
#large_data4=$HOME'/wbai/unipartite_graph/large_data/huge_graph/'
#./shuffle_non_temporal_edges $large_data4 $thread_number
##rm -rf $large_data4'log/'
## shellcheck disable=SC2045
#for file in $(ls $large_data4); do
#  if [ -f $large_data4$file ]; then
#    echo $file
#    ./core_decomposition_test $large_data4 $file $thread_number
#    ./core_insertion_thread_test $large_data4 $file $thread_number
#    ./core_removal_thread_test $large_data4 $file $thread_number
#    ./core_insertion_removal_thread_test $large_data4 $file $thread_number
#  fi
#done

#large_data5=$HOME'/wbai/unipartite_graph/large_data/huge_graph/'
##./shuffle_non_temporal_edges $large_data5 $thread_number
##rm -rf $large_data5'log/'
## shellcheck disable=SC2045
#for file in $(ls $large_data5); do
#  if [ -f $large_data5$file ]; then
#    echo $file
#    ./core_target_insertion_test $large_data5 $file $thread_number
#    ./core_target_removal_test $large_data5 $file $thread_number
#    ./core_target_insertion_removal_test $large_data5 $file $thread_number
#  fi
#done

large_data6=$HOME'/wbai/unipartite_graph/large_data/huge_graph/'
#./shuffle_non_temporal_edges $large_data6 $thread_number
#rm -rf $large_data6'log/'
# shellcheck disable=SC2045
for file in $(ls $large_data6); do
  if [ -f $large_data6$file ]; then
    echo $file
    ./core_insertion_performance_test $large_data6 $file $thread_number
    ./core_removal_performance_test $large_data6 $file $thread_number
  fi
done



