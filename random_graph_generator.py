# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings
import networkx as nx
from networkx.utils.random_sequence import powerlaw_sequence


#preference graph
def generate_ba_graph(e):
    G = nx.barabasi_albert_graph(pow(2, e), 128)
    print_graph(G, 'ba_graph')

#preference graph
def generate_dba_graph(e):
    G = nx.dual_barabasi_albert_graph(pow(2, e), 192, 64, 0.5)
    print_graph(G, 'dba_graph')

def generate_power_law_cluster_graph(e):
    G = nx.powerlaw_cluster_graph(pow(2, e), 128, 0.5)
    print_graph(G, 'pl_graph')


# note that the vertex id start at 0
def print_graph(G, file_name):
    f = open(file_name,'w')

    for u,v in nx.edges(G):
        f.writelines(str(u+1) + "," + str(v + 1) + '\n')
    f.close()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':

    #parameter setting
    e = 19
    file_path = './'

    generate_ba_graph(e)
    generate_dba_graph(e)
    generate_power_law_cluster_graph(e)
    



# See PyCharm help at https://www.jetbrains.com/help/pycharm/
