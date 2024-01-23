#!/bin/bash
PATH='/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/bin'
export PATH

thread_number1=12

small_data=$HOME'/baiwen/bipartite_graph/small_data/'
./wing_shuffle_non_temporal_bipartite_edges $small_data $thread_number1
 #shellcheck disable=SC2045
for file in $(ls $small_data)
do
     if [ -f $small_data$file ]
     then
         echo $file
         ./wing_decomposition_test $small_data $file $thread_number1
  	     ./wing_insertion_test $small_data $file $thread_number1
         ./wing_removal_test $small_data $file $thread_number1
         ./wing_insertion_removal_test $small_data $file $thread_number1
     fi
 done

thread_number2=24
large_data=$HOME'/baiwen/bipartite_graph/large_data/'
./wing_shuffle_non_temporal_bipartite_edges $large_data $thread_number2
# shellcheck disable=SC2045
for file in $(ls $large_data)
do
   if [ -f $large_data$file ]
   then
       echo $file
      ./wing_decomposition_test $large_data $file $thread_number2
      ./wing_insertion_compare_test $large_data $file $thread_number2
      ./wing_removal_compare_test $large_data $file $thread_number2
      ./wing_insertion_size_test $large_data $file $thread_number2
      ./wing_removal_size_test $large_data $file $thread_number2
      ./wing_insertion_thread_test $large_data $file $thread_number2
      ./wing_removal_thread_test $large_data $file $thread_number2
   fi
done


