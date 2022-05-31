#!/bin/bash
PATH='/bin:/sbin:/usr/bin:/usr/sbin:/usr/local/bin:/usr/local/bin'
export PATH

input_data_path=$HOME'/baiwen/original_data/'
output_data_path=$HOME'/baiwen/processed_data/'
thread_number=6

 #shellcheck disable=SC2045
./temporal_graph_process $input_data_path $output_data_path $thread_number



