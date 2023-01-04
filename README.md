# scnu

## Description
This is a commonly used  library of experiments for published papers, 
including several libraries, such as graph, io, random, string, etc.

## Software Architecture
contents in each category

1. include  
   common head files of each component
  
2. lib  
   shared library of each component
   
3. src  
   source code of each component   
   
4. test  
   test project for each component   

## Components
1. container  
contain some useful containers, such thread safe container and various list  

2. bipartite_core  
implement decomposition and maintenance algorithms for (<i>i,j</i>)-core

3. graph  
common graph structure for various types of graph, e.g. general graph, bipartite graph

4. io  
implement input/output functions for various types of graphs

5. logger  
implement a simple_logger system to record the information of projects

6. random  
implement some functions to get random number with the given distributions

7. string  
implement some extra functions for `std::string`

8. thread  
implement some thread classes for supporting multi-threads running

9. time
implement time related classes, e.g. timer for recording running time

...

## Compile
1. `cmake` CMakeList.txt in the root fold of the project and then `make`
2. Executable and library files are located in the lib and test fold respectively
 

 

   


