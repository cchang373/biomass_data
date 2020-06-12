#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 23 14:20:34 2018

@author: cchang373
"""

import imolecule
from rdkit import Chem
from imolecule import generate
import networkx as nx
import matplotlib.pyplot as plt
import networkx.algorithms.isomorphism as iso
import pickle


def generate_mol_dict(data,data_format):
    """generate molecule dictionary including bond,charge,location,element information"""
    mol_dict=eval(generate(data,data_format))
    return mol_dict

def to_dict(Graph):
    mol={}
    atoms=[]
    bonds=[]
    nodes=[]
    
    for node in Graph.nodes():
        atoms.append(Graph.node[node].copy())
        nodes.append(node)
    
    for edge in Graph.edges():
        #print(Graph.edges())
        node_1=edge[0]
        node_2=edge[1]
        bond=[node_1,node_2]
        #print(Graph.edges[node_1,node_2])
        bonds.append({'atoms':bond,'order':Graph.edges[node_1,node_2]['order']}) #treat double/triple bonds as single bond
    
    mol['atoms']=atoms
    mol['bonds']=bonds    
    return mol

def mol_graph(mol_dict):
    """build molecule graph"""
    G_mol=nx.Graph()
    atoms=mol_dict['atoms']
    bonds=mol_dict['bonds']
    
    node_num=[]
    for i,atom in enumerate(atoms):
        s_a={}
        a=atom['element']+str(atom['location'])
        s_a['element']=atom['element']
        s_a['location']=str(atom['location'])
        #s_a['charge']=atom['charge']
        node_num.append((a,s_a))
    G_mol.add_nodes_from(node_num)
    
    for bond in bonds:
        idxs=bond['atoms']
        node_1=atoms[idxs[0]]['element']+str(atoms[idxs[0]]['location'])
        node_2=atoms[idxs[1]]['element']+str(atoms[idxs[1]]['location'])
        #node_1=idxs[0]
        #node_2=idxs[1]
        G_mol.add_edge(node_1,node_2)
        G_mol.edges[node_1,node_2]['order']=bond['order']
        G_mol.edges[node_1,node_2]['atoms']=[node_1,node_2]
    
    return G_mol
                    

def unique_atom_number(mol_dict):
    atoms=mol_dict['atoms']
    unique_atom=[]
    for atom in atoms:
        element=atom['element']
        if element not in unique_atom:
            unique_atom.append(element)
    return len(unique_atom)

"""
all_groups=[]#to save all the recursive bond breaking groups
x=[]#to save the number of bond breakings
group_num=[]#to save the number of accumulate groups after each recursion
groups_recur=[]#to save accumulate groups after each recursion
count_id=[]# to save groups with (element,location)
"""
def subgroup(Graph_list,all_groups=[],x=[],group_num=[],groups_recur=[],count_id=[]):
    
    species=[]
    
    em=iso.categorical_edge_match('order',0)
    #nm_0=iso.categorical_node_match(['element','location'],['X','[0,0,0]'])
    nm_1=iso.categorical_node_match(['element'],['X'])
    
    for G in Graph_list:
        if len(G.edges()) !=0: 
            new_graph=[]
            for edge in G.edges():
                new_G=G.copy()
                new_G.remove_edge(*edge)
                new_graph.append(new_G)
        
            all_frags=[]
            for graph in new_graph:
                if len(graph.edges()) != 0:
                    k=nx.k_components(graph)
                    #print k
                    #for nk in k:
                    frags=k[1]
                    #print frags
                    if len(frags) == 1:
                        for sets in frags:
                            subg=graph.subgraph(sets)
                            #nx.draw(subg,with_labels=True)
                            #plt.show()
                            if subg not in all_frags:
                                all_frags.append(subg)
                            frag_0=[node for node in graph.nodes() if node not in sets]
                            if len(frag_0) != 0:
                                subg_0=graph.subgraph(frag_0)
                                if subg_0 not in all_frags:
                                    all_frags.append(subg_0)
                                    
                    else:
                        for sets in frags:
                            subg=graph.subgraph(sets)
                            if subg not in all_frags:
                                all_frags.append(subg)
                else:
                    for node in graph.nodes():
                        subg=graph.subgraph(node)
                        if subg not in all_frags:
                            all_frags.append(subg)
               
     #       for subg in all_frags:
      #          if True not in [nx.is_isomorphic(subg,G,nm_0,em) for G in species]:
       #             species.append(subg)
                    
            for specy in all_frags:
                if True not in [nx.is_isomorphic(specy,group,nm_1,em) for group in all_groups]:
                    all_groups.append(specy)
                    species.append(specy)
            
            #for subg in all_frags:
             #   if True not in [nx.is_isomorphic(subg,group_l,nm_0,em) for group_l in count_id]:
              #      count_id.append(subg)
                    
        else:
            continue
            
    x.append(1)
    group_num.append((len(x),len(all_groups)))
    groups_recur.append(all_groups)
    print group_num
    print len(groups_recur)
    
    #group_dict=count_num(count_id,all_groups)#save accumulate graph and number for each recursion    
    
    all_group_dict=all_groups_dict(all_groups)
    
    with open('groups_e_'+str(len(x))+'.pickle','wb') as f:
        pickle.dump((species,all_groups,all_group_dict,x,group_num,groups_recur,count_id),f)
    
    if len(species) == 0:
        return all_groups#count_id
    
    return subgroup(species,all_groups=all_groups,x=x,group_num=group_num,groups_recur=groups_recur,count_id=count_id)
    #return all_groups   

def count_num(count_id,all_groups):
    
    group_dict={}
    nm=iso.categorical_node_match('element','X')
    em=iso.categorical_edge_match('order',0)
    
    for id_group in all_groups:
        i=0
        for group in count_id:
            if nx.is_isomorphic(id_group,group,nm,em):
                i +=1
        group_dict[id_group]=i
    return group_dict
        
        
def all_groups_dict(all_groups):
    all_group_dict=[]
    for group in all_groups:
        mol=to_dict(group)
        all_group_dict.append(mol)
        
    return all_group_dict   
                          
#name=input("molecule smiles notation: ")
#molecules=['O=C=O','[H][H]','C(=O)O','[C-]#[O+]','C(=O)OC']
mol_dict=generate_mol_dict('C(C(C(CO)O)O)O','smi')
G_mol=mol_graph(mol_dict)
all_groups=subgroup([G_mol])
#print(all_groups)
all_group_dict=all_groups_dict(all_groups)
#print all_group_dict
#name='sorbic'
with open('groups_'+'C(C(C(CO)O)O)O'+'_2.pickle','wb') as f:
    pickle.dump((all_groups,all_group_dict),f)
#for graph in all_groups:
 #   nx.draw(graph,with_labels=True)
  #  plt.show()
#group_dict=count_num(count_id,all_groups)

"""
with open('groups_'+name+'.pickle','wb') as f:
    #save all the unique groups graph and group dictionary to a pickle file
    pickle.dump((all_groups,all_group_dict),f)



import pandas as pd

V_mol_0=pd.read_excel('my_data.xlsx')
V_mol=V_mol_0.iloc[:,1]
mol_list_0=[]
mol_list=[]
for i in V_mol:
    mol_list_0.append(str(i))
for i in mol_list_0:
    mol_dict=generate_mol_dict(i,'smi')
    G_mol=mol_graph(mol_dict)
    mol_list.append([G_mol])

recursion(mol_list,mol_list_0)

mol_list_0=['C=C','CC','C','O']
mol_list=[]
for i in mol_list_0:
    mol_dict=generate_mol_dict(i,'smi')
    G_mol=mol_graph(mol_dict)
    mol_list.append([G_mol])

recursion(mol_list,mol_list_0)
"""   
        
