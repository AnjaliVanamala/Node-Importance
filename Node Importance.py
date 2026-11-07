"""
Created on Tue Jun 11 10:39:11 2024

@author: Anjali Vanamala
"""
# IMPORTANT!! Things you should change can be found with ctrl-F -> SET THIS. 

import math
import os
import pandas as pd
import numpy as np
import re
from scipy.stats import spearmanr
"""
Node Class definition, but honestly it should probably be called an edge class. 
x is the numerical representation of the first taxa for the edge and y is the 
second. Weight is the absolute value of the correlation coefficient between
the two taxa, with sign being -1 if the correlation is negative and 1 if
positive. 
"""
class Node():
    def __init__(self, num, num2, weight):
        self.x = num 
        self.y = num2
        self.weight = abs(weight)
        # figures out sign.
        if (weight != 0):
            self.sign = (weight)/abs(weight)
        else:
            self.sign = 1
    # str outputs relevant info, mostly used for debugging purposes. 
    def __str__(self):
        return(str(self.x) + " " + str(self.y) + " " + str(self.weight))

"""
Function to build the correlation matrix. Assumes relative abundance inputs.
"""
def make_cor(file): 
    # load in relative abundance file
    taxa = pd.read_csv(file, delimiter="\t")
    print(taxa)
    # s is the number of metadata columns SET THIS. 
    s = 1
    l = len(taxa.columns) - s
    # The start (inclusive) and end (exclusive) of your different treatment/
    # time pairings, if you've already separated before inputting the abundance
    # you can just set it to 0 and then len(taxa.rows). SET THIS
    start = [0, 7, 14, 22]
    end = [7, 14, 22, 46]
    # debug
    #print(taxa.iloc[0:len(taxa.columns), 1:2])
    # Set the range value to the number of time - treatment pairs SET THIS
    for k in range(4):
        # Initialize matrix
        mat = np.zeros((l, l))
        weight = np.zeros((l, l))

        # Extract column names
        names = taxa.columns[s:l+s]

        # Extract and format full names <- this code is iffy depends on what you 
        # want/how qiime outputs stuff can fail, the order is also just the order
        # in the abundance file you inputted if you want to do it manually. 
        #order = [re.split(r"[__//.]+", col) for col in taxa.columns[3:l+3]]
        #print(order)
        #fullname = [f"c_{o[3]}.{o[4]}.{o[5]}.{o[6]}" for o in order]
        #legend = pd.DataFrame({'labels': range(1, 65), 'fullname': fullname})

        # initialize correlations
        cors = np.zeros(l)
        color = []
        # debug/check that your manual starts and ends were correct. Will print
        # to screen.
        print(taxa.iloc[start[k]:end[k], 0:s])
        for i in range(s, l+s):
            # Grabs first non-metadata column. 
            cors[i - s] = 0
            # Does spearman with each column and all other columns. Technically this is redundant 
            # and I could just to columns above, but I like to check it to make sure nothing was
            # selected wrong. 
            for j in range(s, l+s):
                try:
                    cor_stat, p_value = spearmanr(taxa.iloc[start[k]:end[k], i], taxa.iloc[start[k]:end[k], j])
                    if not np.isnan(p_value) and p_value < 0.01 and i != j:
                        cors[i - s] += abs(cor_stat)
                        # Colors are specific to if you're putting this into a graphing software 
                        # that needs edge colors. 
                        color.append("red" if cor_stat > 0 else "blue")
                        mat[j - s, i - s] = 1
                        weight[j - s, i - s] = cor_stat
                except:
                    continue

        # Save results
        np.fill_diagonal(mat, 0)
        np.fill_diagonal(weight, 0)
        print(weight)
        # File to save to, will just go from 0-n where n is the number of different 
        # treatment time pairs you have SET THIS 
        pd.DataFrame(weight).to_csv(f"Corr_matrix{k}.txt", sep="\t", index=False)
    
    # To save labels if you've uncommented that. 
    #pd.DataFrame(legend).to_csv("labels-SpringMLabelAnna.txt", sep="\t", index=False)

"""
Takes the saved correlation matrix and loads it back in. Yes this is kind of redundant,
but I like to be able to use these functions for other things! 
"""
def load_adj(p):
    lines = []
    # SET THIS to the same as the make_cor function.
    file = open(f"Corr_matrix{p}.txt", "r").readlines()
    # Goes through and adds each edge. 
    for line in file[1:len(file)]:
        lin = line.rstrip().split('\t')
        lines.append([float(x) for x in lin])
    adj = [[] for i in range(len(lines))]
    for i in range(len(lines)):
       for j in range(len(lines[i])):
           if (lines[i][j] != 0):
               adj[i].append(Node(i, j, lines[i][j]))
    return(adj)
    
"""
Graph class definition. This started as a prim's algorithm project, this is a remnant 
you can delete it. 
"""
class Graph():
    def __init__(self, taxa):
        self.graph = [[0 for column in range(taxa)]
                      for row in range(taxa)]
    
# Sort adjacencies near zero 
    def max_edge(self, nodelist, adj):
        maxe = Node(-1, -1, 0)
        for i in nodelist:
            for j in range(len(adj[i])):
                if (adj[i][j].y not in nodelist and adj[i][j].weight > maxe.weight):
                    maxe = adj[i][j]
        return(maxe)
            
# Add to graph 
    def add(self, maxe):
        self.graph[maxe.x][maxe.y] = maxe.weight * maxe.sign
        return()

# If find_max returns 0, add next node that's not in list 
    def next_n(self, taxanum, nodelist):
        for i in range(taxanum):
            if (i not in nodelist):
                nodelist.append(i)
                break
        return()
# When all original nodes are in the nodes list stop
    def Prims(self, adj):
        taxa = len(adj)
        nodelist = [0]
        while (len(nodelist) < taxa):
            maxe = self.max_edge(nodelist, adj)
            if maxe.weight != 0:
                self.add(maxe)
                nodelist.append(maxe.y)
            else:
                self.next_n(taxa, nodelist)
        return(self.graph)

"""
Gets average abundance for each taxa. We ended up using node importance for 
both size and color of the network but you can use abundance for one 
if you'd like. 
"""
def get_abun(file):
    # Same thing as the previous s. SET THIS
    s = 1
    lines = []
    file = open(file, "r").readlines()
    for line in file[1:len(file)]:
        lis = line.rstrip().split('\t')
        lines.append([float(x) for x in lis[s:len(lis)]])
    lines = list(map(list, zip(*lines)))
    lines = [sum(x)/len(x) for x in lines]
    print(lines)
    return(lines)
# Main 
if __name__ == '__main__':
    # The abundance file SET THIS, file format expectation is that columns 
	# are taxa and rows are treatments. 
    file = "RelativeAbundanceFile.txt"
    make_cor(file)

    # debug
    #print("correlations done")
    # Range value is the number of treatment/time pairs you have. SET THIS
    for k in range(4):
        adj = load_adj(k)
        adj0 = load_adj(0)
        adj1 = load_adj(1)
        adj2 = load_adj(2)
        adj3 = load_adj(3)
        # debug
        #print("load adj done")
        # Prims algorithm remnants 
        """
        adj = [[], [], [], [], []]
            matrix = [[0, 2, 0, 6, 0],
            [2, 0, 3, 8, 5],
                [0, 3, 0, 0, 7],
                [6, 8, 0, 0, 9],
                [0, 5, 7, 9, 0]]
        
        for i in range(len(matrix)):
            for j in range(len(matrix[i])):
                adj[i].append(Node(i, j, matrix[i][j]))
                adj[j].append(Node(j, i, matrix[i][j]))
        
        graph = Graph(len(adj))
        print("graph done")
        prims = graph.Prims(adj)
        print("prims done")
        
        print("abundance done")
        f = open(f"prims_graph{k}.txt", "w") 
        for i in prims:
            for j in i:
                f.write(str(j) + "\t")
            f.write("\n")
        f.close()
        """
        abun = get_abun(file)
        # Weights for nodes 
        nodeImp = []
        for i in range(len(adj)):
            summed = 0
            # Literally just the sum of the absolute value of connected node 
            # correlation coefficients. 
            for j in adj[i]:
                summed += j.weight
                    # debug
                    #print("CHECK", i, j.y, len(adj[i]), len(adj[j.y]), (j.weight), ((len(adj[i]) * len(adj[j.y]))/(((1-j.weight)*100)**2)))
                    #print("correlation", i, j)
                #print(i, math.log(summed/len(adj[i])))
            nodeImp.append(summed)
            
        # Can normalize importances (I suggest against this as correlations are already essentially 
        # normalized to 0-1) 
        #nodeImp = [imp/max(nodeImp) for imp in nodeImp]

        # Prims stuff again.     
        """
        f = open(f"prims_graph_flat{k}.txt", "w") 
        f.write(f"Taxa1.{k}\tTaxa1Imp.{k}\tTaxa2.{k}\tCorrelation\n")
        for i in range(len(prims)):
            for j in range(len(prims[i])):
                if (prims[i][j] != 0):
                    #print(i)
                    f.write(str(i) + "\t" + str(nodeImp[i]) + "\t" + str(j) + "\t" + str(prims[i][j]) + "\n")
                    #f.write(str(j) + "\t" + str(j) + "\t" + str(i) + "\t" + str(prims[i][j]) + "\n")
                    if i not in list1:
                        list1.append(i)
                    if j not in list1:
                        list1.append(j)
            f.write(str(i) + "\t" + str(nodeImp[i]) + "\t" + str(i) + "\t" + str(prims[i][i]) + "\n")
        f.close()
        """
        # Flattens and prints the correlation matrices into a form that can be directly used in Cytoscape
        # Also outputs a list of node importances. 
        list1 = []
        list2 = []
        # Cytoscape file SET THIS
        f = open(f"adjacency_matrix_flat{k}.txt", "w") 
        # Node importance file SET THIS
        n = open(f"node_importance{k}.txt", "w")
        f.write(f"Taxa1.{k}\tTaxa1Importance.{k}\tTaxa2.{k}\tTaxa1AvgAbundance.{k}\tCorrelation.{k}\n")
        print(len(adj0), len(adj2))
        for i in range(len(adj)):
            for j in range(len(adj[i])):
                if (adj[i][j] != 0):
                    if (adj[i][j].x, adj[i][j].y) not in list2 and (adj[i][j].y, adj[i][j].x) not in list2:
                        check = 1
                        for p in range(len(adj0[i])):
                            if (adj0[i][p].y == adj[i][j].y and k != 0):
                                check = 0
                        for p in range(len(adj1[i])):
                            if (adj1[i][p].y == adj[i][j].y and k != 1):
                                check = 0
                        for p in range(len(adj2[i])):
                            if (adj2[i][p].y == adj[i][j].y and k != 2):
                                check = 0
                        for p in range(len(adj3[i])):
                            if (adj3[i][p].y == adj[i][j].y and k != 3):
                                check = 0
                        if check:
                            f.write(str(adj[i][j].x) + "\t" + str(nodeImp[adj[i][j].x]) + "\t" + str(adj[i][j].y) + "\t" + str(abun[adj[i][j].x]) + "\t" + str(adj[i][j].weight*adj[i][j].sign) + "\n")
                            #print(str(adj[i][j].x) + "\t" + str(adj[i][j].y) + "\t" + str(nodeImp[adj[i][j].x]) + "\n")
                            list2.append((adj[i][j].x, adj[i][j].y))
            n.write(str(nodeImp[i]) + "\n")
            # Explicitly setting self correlations to 0
            f.write(str(i) + "\t" + str(nodeImp[i]) + "\t" + str(i) + "\t" + str(abun[i]) + "\t" + str(0) + "\n")        
        f.close()
                
