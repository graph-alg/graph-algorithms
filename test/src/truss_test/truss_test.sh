#!/bin/bash
PATH='/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/bin'
export PATH

thread_number=16

small_data1=$HOME'/wbai/unipartite_graph/processed_data/non_temporal_data/small_data/'
./truss_shuffle_non_temporal_edges $small_data1 $thread_number
 #shellcheck disable=SC2045
for file in $(ls $small_data1)
do
     if [ -f $small_data1$file ]
     then
         echo $file
         ./truss_decomposition_test $small_data1 $file $thread_number
  	     ./truss_insertion_test $small_data1 $file $thread_number
         ./truss_removal_test $small_data1 $file $thread_number
         ./truss_insertion_removal_test $small_data1 $file $thread_number
     fi
 done

small_data2=$HOME'/wbai/unipartite_graph/processed_data/unique_temporal_data/small_data/'
  #shellcheck disable=SC2045
for file in $(ls $small_data2)
do
      if [ -f $small_data2$file ]
      then
          echo $file
          ./truss_decomposition_test $small_data2 $file $thread_number
   	     ./truss_insertion_test $small_data2 $file $thread_number
          ./truss_removal_test $small_data2 $file $thread_number
          ./truss_insertion_removal_test $small_data2 $file $thread_number
      fi
  done

 random_data=$HOME'/wbai/unipartite_graph/processed_data/non_temporal_data/random_data/'
./truss_shuffle_non_temporal_edges $random_data $thread_number
# shellcheck disable=SC2045
for file in $(ls $random_data)
  do
     if [ -f $random_data$file ]
     then
        echo $file
       ./truss_decomposition_test $random_data $file $thread_number
       ./truss_insertion_test $random_data $file $thread_number
       ./truss_removal_test $random_data $file $thread_number
        ./truss_insertion_removal_test $random_data $file $thread_number
        ./truss_insertion_compare_test $random_data $file $thread_number
        ./truss_removal_compare_test $random_data $file $thread_number
        ./truss_insertion_size_test $random_data $file $thread_number
        ./truss_removal_size_test $random_data $file $thread_number
        ./truss_insertion_thread_test $random_data $file $thread_number
        ./truss_removal_thread_test $random_data $file $thread_number
     fi
  done

large_data=$HOME'/wbai/unipartite_graph/processed_data/non_temporal_data/large_data/'
./truss_shuffle_non_temporal_edges $large_data $thread_number
# shellcheck disable=SC2045
for file in $(ls $large_data)
do
   if [ -f $large_data$file ]
   then
       echo $file
     ./truss_decomposition_test $large_data $file $thread_number
     ./truss_insertion_test $large_data $file $thread_number
     ./truss_removal_test $large_data $file $thread_number
      ./truss_insertion_removal_test $large_data $file $thread_number
      ./truss_insertion_compare_test $large_data $file $thread_number
      ./truss_removal_compare_test $large_data $file $thread_number
      ./truss_insertion_size_test $large_data $file $thread_number
      ./truss_removal_size_test $large_data $file $thread_number
      ./truss_insertion_thread_test $large_data $file $thread_number
      ./truss_removal_thread_test $large_data $file $thread_number
   fi
done


