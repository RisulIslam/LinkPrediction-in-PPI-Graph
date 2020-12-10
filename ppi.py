import sys, codecs
from unidecode import unidecode
import re
import matplotlib.pyplot as plt

import networkx as nx
from statistics import median
from collections import Counter
import numpy as np
import matplotlib.patches as mpatches


G = nx.Graph()
GM = nx.Graph()
#threshold=0.0

fh=codecs.open("modified human 500.txt","r",encoding='utf-8')
linesh=fh.readlines()

fm=codecs.open("modified mouse 500.txt","r",encoding='utf-8')
linesm=fm.readlines()
def find_cc(x,y):
    p=[]
    c=[]
    cc=[]
    cum=0
    for i in range(0,len(x)):
        prob=y[i]/np.sum(y)
        p.insert(len(p),prob)
        cum=cum+prob
        c.insert(len(c),cum)
        cc.insert(len(cc),1-cum)
        #print("P: ",prob,"   C: ",cum)
    return c

green_patch = mpatches.Patch(color='green', label='Human')
red_patch = mpatches.Patch(color='red', label='Mouse')

def draw_distribution(x,y,p,q,xlabel,ylabel):
    
    ch=find_cc(x,y)
    #print(ch)
    plt.errorbar(x,ch,fmt='g.')
    cm=find_cc(p,q)
    plt.errorbar(p,cm,fmt='r.')
    #plt.title("# of 6-mers with distinct frequency vs 6-mer frequency for Random DNA sequence")
    
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend(handles=[green_patch,red_patch],loc=4)
    #plt.show()
    directory="fig500cent/cdf "+xlabel+" VS "+ylabel
    plt.savefig(directory)
    plt.close()

def find_distribution_from_counter(L):
    
    c=Counter(L)
    items=[]
    freq=[]
    for k, v in sorted(c.items(), key=lambda x: x[0], reverse=False):
        items.append(k)
        freq.append(v)
    
    #items=list(c.keys())
    #freq=list(c.values())
    #print(items)
    #print(freq)
    return items,freq

def find_basic_statistics(G,GM):
    # number of nodes and edges
    nonode=G.number_of_nodes()
    noedge=G.number_of_edges()
    nonodem=GM.number_of_nodes()
    noedgem=GM.number_of_edges()
    print("\n\nHuman basic info colection\n")
    print("Number of nodes h: ",nonode," Number of edges h: ",noedge)
    
    # human basic info collection
    nodelist=list(G.nodes)
    maxd=0
    maxdnodes=[]
    dlist=[]
    for node in nodelist:
        d=G.degree[node]
        dlist.append(d)
        if d>maxd:
            maxd=d
    for node in nodelist:
        d=G.degree[node]
        if d==maxd:
            maxdnodes.append(node)
    print("Max degree: ",maxd)
    print("Max degree nodes list:")
    print(maxdnodes)
    mediand=median(dlist)
    print("Median degree: ",mediand)
    avgd=sum(dlist)/len(dlist)
    print("Avg degree: ",avgd)
    mind=min(dlist)
    print("Min degree: ",mind)
    no_connected_component=nx.algorithms.components.number_connected_components(G)
    print("Number of connected component: ",no_connected_component)

    
    print("\n\nmouse basic info colection\n")
    print("Number of nodes m: ",nonodem," Number of edges m: ",noedgem)
    nodelist=list(GM.nodes)
    maxd=0
    maxdnodes=[]
    dlist=[]
    for node in nodelist:
        d=GM.degree[node]
        dlist.append(d)
        if d>maxd:
            maxd=d
    for node in nodelist:
        d=GM.degree[node]
        if d==maxd:
            maxdnodes.append(node)
    print("Max degree m: ",maxd)
    print("Max degree nodes list m:")
    print(maxdnodes)
    mediand=median(dlist)
    print("Median degree m: ",mediand)
    avgd=sum(dlist)/len(dlist)
    print("Avg degree m: ",avgd)
    mind=min(dlist)
    print("Min degree m: ",mind)
    no_connected_component=nx.algorithms.components.number_connected_components(GM)
    print("Number of connected component m: ",no_connected_component)

def find_list_from_map(map_d_cent):
    
    d_cent_list=[]
    for key in map_d_cent.keys():
        d_cent_list.append(map_d_cent[key])
    return d_cent_list

def find_centralities(G,GM):
    
    print("DC>>>")
    map_d_cent=nx.algorithms.centrality.degree_centrality(G)
    d_cent_list=find_list_from_map(map_d_cent)
    map_d_centm=nx.algorithms.centrality.degree_centrality(GM)
    d_cent_listm=find_list_from_map(map_d_centm)
    x,y=find_distribution_from_counter(d_cent_list)
    p,q=find_distribution_from_counter(d_cent_listm)
    draw_distribution(x,y,p,q,"Degree Centrality","Frequency")
    
    print("EVC>>>")    
    map_ev_cent=nx.eigenvector_centrality(G)
    ev_cent_list=find_list_from_map(map_ev_cent)
    x,y=find_distribution_from_counter(ev_cent_list)
    map_ev_cent=nx.eigenvector_centrality(GM)
    ev_cent_list=find_list_from_map(map_ev_cent)
    p,q=find_distribution_from_counter(ev_cent_list)
    draw_distribution(x,y,p,q,"Eigenvector Centrality","Frequency")
    
    print("CLC>>>")
    map_cl_cent=nx.algorithms.centrality.closeness_centrality(G)
    cl_cent_list=find_list_from_map(map_cl_cent)
    x,y=find_distribution_from_counter(cl_cent_list)
    print("CLC>>>>>>>")
    map_cl_cent=nx.algorithms.centrality.closeness_centrality(GM)
    cl_cent_list=find_list_from_map(map_cl_cent)
    p,q=find_distribution_from_counter(cl_cent_list)
    draw_distribution(x,y,p,q,"Closeness Centrality","Frequency")

    print("BC>>>")    
    map_b_cent=nx.algorithms.centrality.betweenness_centrality(G)
    b_cent_list=find_list_from_map(map_b_cent)
    x,y=find_distribution_from_counter(b_cent_list)
    print("BC>>>>>>>")
    map_b_cent=nx.algorithms.centrality.betweenness_centrality(GM)
    b_cent_list=find_list_from_map(map_b_cent)
    p,q=find_distribution_from_counter(b_cent_list)
    draw_distribution(x,y,p,q,"Betweenness Centrality","Frequency")

    print("SGC>>>")
    map_sg_cent=nx.algorithms.centrality.subgraph_centrality(G)
    sg_cent_list=find_list_from_map(map_sg_cent)
    x,y=find_distribution_from_counter(sg_cent_list)
    print("SGC>>>>>>>")
    map_sg_cent=nx.algorithms.centrality.subgraph_centrality(GM)
    sg_cent_list=find_list_from_map(map_sg_cent)
    p,q=find_distribution_from_counter(sg_cent_list)
    draw_distribution(x,y,p,q,"Subgraph Centrality","Frequency")
    

def find_radius_center_diameter(G):
    print("Connected components >>>>")
    sgs=list(nx.connected_component_subgraphs(G))
    
    #print(sgs[2].nodes)
    r_list=[]
    d_list=[]
    c_lists=[]
    for sg in sgs:
        print("Number of nodes in this connected components: ",sg.number_of_nodes())
        r=nx.algorithms.distance_measures.radius(sg)
        r_list.append(r)
        d=nx.algorithms.distance_measures.diameter(sg)
        d_list.append(d)
        center_list=nx.algorithms.distance_measures.center(sg)
        c_lists.append(center_list)
    print("Center node list: ")
    print(c_lists)
    print("Rasius: ",r_list)
    print("Diameter: ",d_list)

def find_intersection_lower_pred(elower,pred_edges):
    count=0
    common_list=[]
    for u,v in pred_edges:
        if (u,v) in elower or (v,u) in elower:
            common_list.append((u,v))
            count+=1
    return count,common_list

def find_union(l1,l2):
    res=[]
    for u,v in l1:
        if (u,v) not in l2 or (v,u) not in l2:
            res.append((u,v))
    for u,v in l2:
        res.append((u,v))
    return res

def make_new_G_on_predicted_links(G,threshold):
    print("Predicting links>>>")
    n=list(G.nodes)
    edge_pairs_list=G.edges
    #print("Old Edges:",edge_pairs_list)
    
    print("pairs making")
    
    n_pairs=[]
    pred_edgesh=[]
    pred_another=[]
    for i in range(0,len(n)):
        
        
       
        for j in range(i+1,len(n)):
            n1=n[i]
            n2=n[j]
            pair=(n1,n2)
            if pair not in edge_pairs_list:
                n_pairs.append(pair)
            if (j+1)%5000==0:
                break
        if (i+1)%10==0:    
            print("probability calculating: ",i)    
            #predj = nx.resource_allocation_index(G,n_pairs) #Resource allocation similarity
            #predr =nx.jaccard_coefficient(G,n_pairs) # Jaccard Similarity
            predj = nx.adamic_adar_index(G,n_pairs)
            print("Adding new edges")
            """
            map_pair_pr={}
            for u,v,pr in predr:
                if pr>0.5:
                    map_pair_pr[(u,v)]=pr
            print("map len: ",len(map_pair_pr))
            """
            
            for u,v,pj in predj:
                if pj>threshold:
                    pred_edgesh.append((u,v))
                    #G.add_edge(u,v,weight=1)
                    #print("Finally added")
            """
            for u,v,pj in predr:
                if pj>threshold:
                    pred_another.append((u,v))
                    #G.add_edge(u,v,weight=1)
            count,union=find_intersection_lower_pred(pred_edgesh,pred_another)# UNION nad INTERSECTION
            print(len(pred_edgesh),len(pred_another),len(union))
            pred_edgesh=union
            """
            
            n_pairs=[]
        if (i+1)%10==0:
            break
        
        
    return pred_edgesh
                    

def find_maximalclique(G,Gm):
    print("MC G")
    I = nx.algorithms.clique.find_cliques(G)
    hclq_nodes=[]
    hclq_number=0
    m=0
    hc_list=[]
    for clique_nodes in I:
        #print("Clique vertices : ",clique_nodes)
        m+=1
        hc_list.append(len(clique_nodes))
        if len(clique_nodes)>hclq_number:
            hclq_number=len(clique_nodes)
            #hclq_nodes.append(clique_nodes)
            hclq_nodes=clique_nodes
        if m%5000==0:
            print(m)
        if m>50000:
            break
    x,y=find_distribution_from_counter(hc_list)
    print("Human largest maximal clique number:",hclq_number)
    print("Human clique vertices:")
    print(hclq_nodes)

    print("MC GM")
    I = nx.algorithms.clique.find_cliques(GM)
    mclq_nodes=[]
    mclq_number=0
    m=0
    mc_list=[]
    for clique_nodes in I:
        #print("Clique vertices : ",clique_nodes)
        m+=1
        mc_list.append(len(clique_nodes))
        if len(clique_nodes)>mclq_number:
            mclq_number=len(clique_nodes)
            #mclq_nodes.append(clique_nodes)
            mclq_nodes=clique_nodes
        if m%5000==0:
            print(m)
        if m>50000:
            break
    p,q=find_distribution_from_counter(mc_list)
    draw_distribution(x,y,p,q,"Maximal Clique Number","Frequency Probability")
    print("Mouse largest maximal clique number:",mclq_number)
    print("Mouse clique vertices:")
    print(mclq_nodes)

    return hclq_nodes,mclq_nodes

def create_low_edgelist_for_handm():
    
    max_new=5000000

    
    print("Human>>")
    fh=codecs.open("modified human lower500.txt","r",encoding='utf-8')
    linesh=fh.readlines()
    lower_edge_listh=[]
    k=0
    g=nx.Graph()
    for line in linesh:
        line=line.replace('\n','')
        line=line.replace('\r','')
        #print(len(line))
        t=line.split(' ')
        #print(len(t[1]),"and",len(t[0]))
        g.add_edge(t[0],t[1],weight=1)
        k+=1
        if k>max_new:
            break
    print("Mouse>>")    
    #print("Lower 500 edges number:",len(lower_edge_listh))
    fm=codecs.open("modified mouse lower500.txt","r",encoding='utf-8')
    linesm=fm.readlines()
    lower_edge_listm=[]
    k=0
    gm=nx.Graph()
    for line in linesm:
        line=line.replace('\n','')
        line=line.replace('\r','')
        #print(len(line))
        t=line.split(' ')
        #print(len(t[1]),"and",len(t[0]))
        gm.add_edge(t[0],t[1],weight=1)
        k+=1
        if k>max_new:
            break
    #print("Lower 500 edges number:",len(lower_edge_listm))

    return g.edges, gm.edges


            
        

def make_graph(linesh,linesm):
    elementnumber=1000000
    print("Reading human data")
    j=1
    for line in linesh:
        line=line.replace('\n','')
        line=line.replace('\r','')
        #print(len(line))
        t=line.split(' ')
        #print(len(t[1]),"and",len(t[0]))
        G.add_edge(t[0],t[1],weight=1)
        
        j+=1
        if j>elementnumber:
            break
    print("Reading mouse data")
    j=1
    for line in linesm:
        line=line.replace('\n','')
        line=line.replace('\r','')
        #print(len(line))
        t=line.split(' ')
        #print(len(t[1]),"and",len(t[0]))
        GM.add_edge(t[0],t[1],weight=1)
        
        j+=1
        if j>elementnumber:
            break
    find_basic_statistics(G,GM)
    
    #find_centralities(G,GM)
    print("\n\nFinding radius centers for human\n")
    #find_radius_center_diameter(G)
    print("\n\nFinding radius centers for Mouse\n")
    #find_radius_center_diameter(GM)
    #hclq_nodes,mclq_nodes=find_maximalclique(G,GM)

    print("Creating lower edge lists")
    elowerh,elowerm=create_low_edgelist_for_handm()
    print("Lower edge list human: ",len(elowerh))
    print("Lower edge list mouse: ",len(elowerm))
    threshold=0.0
    while threshold<0.9:
        pred_edgesh=make_new_G_on_predicted_links(G,threshold)
        pred_edgesm=make_new_G_on_predicted_links(GM,threshold)
        
        common_counth,common_listh=find_intersection_lower_pred(elowerh,pred_edgesh)
        common_countm,common_listm=find_intersection_lower_pred(elowerm,pred_edgesm)

        totalpred_h=len(pred_edgesh)+1
        common_percent=(common_counth*100)/totalpred_h
        
        print("For human threshold > ",threshold)
        print("Predicted new links human: ",totalpred_h)
        print("Common count percentage:",common_percent)
        print("False predicted percentage:",100-common_percent)

        totalpred_m=len(pred_edgesm)+1
        common_percent=(common_countm*100)/totalpred_m
        print("For mouse threshold > ",threshold)
        print("Predicted new links mouse: ",totalpred_m)
        print("Common count percentage:",common_percent)
        print("False predicted percentage:",100-common_percent)

        threshold +=0.1
    
    
    """
    newhclq_nodes,newmclq_nodes=find_maximalclique(Gmod,GMmod)
    diffhclq_nodes=list(set(newhclq_nodes)-set(hclq_nodes))
    diffmclq_nodes=list(set(newmclq_nodes)-set(mclq_nodes))
    print("Difference h nodes: ",len(diffhclq_nodes))
    print(diffhclq_nodes)
    print("Difference m nodes: ",len(diffmclq_nodes))
    print(diffmclq_nodes)
    """
    
    
    
    
    

make_graph(linesh,linesm)
        
