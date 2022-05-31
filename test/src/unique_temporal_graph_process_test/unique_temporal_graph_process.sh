#!/bin/bash
PATH='/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/bin'
export PATH

input_data_path=$HOME'/wbai/unipartite_graph/original_data/unique_temporal_graph/'
output_data_path=$HOME'/wbai/unipartite_graph/processed_data/unique_temporal_graph/'
thread_number=14

 #shellcheck disable=SC2045
./unique_temporal_graph_process_test $input_data_path $output_data_path $thread_number



