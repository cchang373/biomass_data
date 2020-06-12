#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  8 14:41:38 2018

@author: cchang373
"""

import imolecule
from imolecule import generate
import networkx as nx
import matplotlib.pyplot as plt
import networkx.algorithms.isomorphism as iso
import pickle
from pickle import load
import json
from json import load

def to_dict(smi):
    mol_dict=eval(generate(smi,'smi'))
    return mol_dict

def to_graph(mol_dict):
    G_mol=nx.Graph()
    atoms=mol_dict['atoms']
    bonds=mol_dict['bonds']
    
    node_num=[]
    for i,atom in enumerate(atoms):
        s_a={}
        element=atoms[i]['element']
        location=atoms[i]['location']
        #charge=atoms[i]['charge']
        a=atom['element']+str(location)
        s_a['element']=element
        s_a['location']=location
        #s_a['charge']=charge
        node_num.append((a,s_a))
    G_mol.add_nodes_from(node_num)
    
    for bond in bonds:
        idxs=bond['atoms']
        order_n=bond['order']
        node_1=node_num[idxs[0]][0]
        node_2=node_num[idxs[1]][0]
        #node_1=atoms.index({'location':idxs[0][1:],'element':idxs[0][0]})
        #node_2=atoms.index({'location':idxs[1][1:],'element':idxs[1][0]})
        G_mol.add_edge(node_1,node_2,order=order_n)
        #G_mol.add_edge(node_1,node_2)
    return G_mol


def one_bond_breaking(mol_name,G_mol_list,bond_breaking={'CH':[],'OH':[],'CO':[],'CC':[],'CCd':[],'COd':[]},
                      original={'CH':[],'OH':[],'CO':[],'CC':[],'CCd':[],'COd':[]},bonds=[],nums=0):
    """
    G_mol_list should be one original molecule graph
    this function return one bond breaking results of this molecule
    """
    
    for G_mol in G_mol_list:
        #break one chemical bond
        #nx.draw(G_mol,with_labels=True)
        #plt.show()
        new_mol=[]
        elements=[]
        #print G_mol.edges(data=True)
        if len(G_mol.edges()) > 0:
            for edge in G_mol.edges():
                element_1=edge[0]
                element_2=edge[1]
                new_G_mol=G_mol.copy()
                new_G_mol.remove_edge(*edge)
        
                elements.append([element_1,element_2])
                new_mol.append(new_G_mol)
    
        all_products=[]
        element_products=[]
    
        em=iso.categorical_edge_match('order',0)
        nm=iso.categorical_node_match(['element'],['X'])
        
        for num,i in enumerate(new_mol):
            if len(i.nodes()) > 2:
                k=nx.k_components(i)
                products=[]
                for nk in k:
                    frags=k[nk]
                    if len(frags) == 1:
                        for sets in frags:
                            subg=i.subgraph(sets)
                            frag_0=[node for node in i.nodes() if node not in sets]
                            subg_0=i.subgraph(frag_0)
                            products.append(subg)
                            products.append(subg_0)
                    else:
                        for sets in frags:
                            subg=i.subgraph(sets)
                            products.append(subg)
                
                if len(products) != 0:                
                    element_products.append(elements[num])        
                    all_products.append(products)
                    
            elif len(i.nodes()) == 2:
                element_products.append(elements[num])
                products=[]
                for node in i.nodes():
                    subg=i.subgraph(node)
                    products.append(subg)
                all_products.append(products)
                
            #print element_products
        
        for num,element in enumerate(element_products):
            if 'C' in element[0] and 'H' in element[1] :
                bond_breaking['CH'].append(all_products[num])
                original['CH'].append(G_mol)
            elif 'C' in element[1] and 'H' in element[0]:
                bond_breaking['CH'].append(all_products[num])
                original['CH'].append(G_mol)
            elif 'C' in element[0] and 'O' in element[1]:
                if G_mol[element[0]][element[1]]['order'] == 2:
                    bond_breaking['COd'].append(all_products[num])
                    original['COd'].append(G_mol)
                elif G_mol[element[0]][element[1]]['order'] == 1:
                    bond_breaking['CO'].append(all_products[num])
                    original['CO'].append(G_mol)
            elif 'C' in element[1] and 'O' in element[0]:
                if G_mol[element[0]][element[1]]['order'] == 2:
                     bond_breaking['COd'].append(all_products[num])
                     original['COd'].append(G_mol)
                elif G_mol[element[0]][element[1]]['order'] == 1:
                    bond_breaking['CO'].append(all_products[num])
                    original['CO'].append(G_mol)
            elif 'C' in element[0] and 'C' in element[1]:
                #print element
                if G_mol[element[0]][element[1]]['order'] == 2:
                    bond_breaking['CCd'].append(all_products[num])
                    original['CCd'].append(G_mol)
                elif G_mol[element[0]][element[1]]['order'] == 1:
                    bond_breaking['CC'].append(all_products[num])
                    original['CC'].append(G_mol)
            else:
                bond_breaking['OH'].append(all_products[num])
                original['OH'].append(G_mol)
        #print bond_breaking
    
        bond={}
        bond['CC']=bond_breaking['CC'][:]
        bond['CH']=bond_breaking['CH'][:]
        bond['CO']=bond_breaking['CO'][:]
        bond['OH']=bond_breaking['OH'][:]
        bond['CCd']=bond_breaking['CCd'][:]
        bond['COd']=bond_breaking['COd'][:]
        #print bond
        
        original_m={}
        original_m['CC']=original['CC'][:]
        original_m['CH']=original['CH'][:]
        original_m['CO']=original['CO'][:]
        original_m['OH']=original['OH'][:]
        original_m['CCd']=original['CCd'][:]
        original_m['COd']=original['COd'][:]
        
        origins={'CH':[],'OH':[],'CO':[],'CC':[],'CCd':[],'COd':[]}
        for key in bond_breaking.keys():
            for num,products in enumerate(bond_breaking[key]):
                bond[key].remove(products)
                origin=original_m[key][num]
                if len(bond[key]) != 0:
                    if True in [nx.is_isomorphic(products[0],ps[1],nm,em) and nx.is_isomorphic(products[1],ps[0],nm,em) for ps in bond[key]]:
                        origins[key].append(origin)
                    elif True in [nx.is_isomorphic(products[0],ps[0],nm,em) and nx.is_isomorphic(products[1],ps[1],nm,em) for ps in bond[key]]:
                        origins[key].append(origin)
                    
                    else:
                        bond[key].append(products)
                else:
                    bond[key].append(products)
        
        for key in origins.keys():
            for i in origins[key]:
                original_m[key].remove(i)
    
    new_mol_list=[]
    for key in bond.keys():
        for products in bond[key]:
            for graph in products:
                new_mol_list.append(graph)
    
    bond_l=[]           
    for key in bond.keys():
        #print key
        l=len(bond[key])
        bond_l.append(l)
    bonds.append(bond_l)
    #print bonds
    if nums > 1:
        if bonds[nums] == bonds[nums-1]:
            return bond,original_m
    nums += 1
    
    with open(mol_name+'/rxn/%s.pickle' % (str(nums)),'w+') as f:
        pickle.dump([new_mol_list,bond,original_m,bonds,nums],f)
    
    return one_bond_breaking(mol_name,G_mol_list=new_mol_list,bond_breaking=bond,original=original_m,bonds=bonds,nums=nums)

def unique_reactions(mol_name_list,json_file_list):
    """
    this function is used to find unique reactions for a list of molecules
    input json_file_list should be a list of json files [vectors,reactants,products]
    input mol_name_list should be a list of molecule names corresponds to json_file_list
    """
    reactions_p={'CH':[],'OH':[],'CO':[],'CC':[],'CCd':[],'COd':[]}
    reactions_r={'CH':[],'OH':[],'CO':[],'CC':[],'CCd':[],'COd':[]}
    reactions_v={'CH':[],'OH':[],'CO':[],'CC':[],'CCd':[],'COd':[]}
    reactions_n={'CH':[],'OH':[],'CO':[],'CC':[],'CCd':[],'COd':[]}
    for num,json_file in enumerate(json_file_list):
        mol_name=mol_name_list[num]
        p_file=json.load(open(json_file))[2]
        r_file=json.load(open(json_file))[1]
        v_file=json.load(open(json_file))[0]
        for key in p_file.keys():
            for idxs,products in enumerate(p_file[key]):
                if True in [products[0] == product[0] and products[1] == product[1] for product in reactions_p[key]]:
                    for product in reactions_p[key]:
                        if products[0] == product[0] and products[1] == product[1]:
                            #print product
                            products_idxs=reactions_p[key].index(products)
                            reactions_n[key][products_idxs].append(mol_name)
                            #print reactions_n[key][products_idxs]
                elif True in [products[0] == product[1] and products[1] == product[0] for product in reactions_p[key]]:
                    for product in reactions_p[key]:
                        if products[0] == product[1] and products[1] == product[0]:
                            #print product
                            products_idxs=reactions_p[key].index(product)
                            reactions_n[key][products_idxs].append(mol_name)
                            #print reactions_n[key][products_idxs]
                else:
                    reactions_p[key].append(products)
                    reactions_r[key].append(r_file[key][idxs])
                    reactions_v[key].append(v_file[key][idxs])
                    reactions_n[key].append([mol_name])
    return [reactions_v,reactions_r,reactions_p,reactions_n]

#reactions_v,reactions_r,reactions_p,reactions_n=unique_reactions(['fumaric','propionic'],
 #                                                               ['groups_rxn/fumaric/rxn/fumaric_10_0.json',
  #                                                               'groups_rxn/propionic/rxn/propionic_10_0.json'])
#with open('groups_rxn/unique_reaction/fumaric_propionic_10_0.json','w+') as f:
 #   json.dump([reactions_v,reactions_r,reactions_p,reactions_n],f)

def bond_to_xyz(react_dict,product_dict,mol_name,path):
    """
    convert bond_breaking results to xyz file
    the input pickle_file should be a pickle file contains
    {'C-O':[],'C-H':[],...} dictionary
    """
    #product_dict=load(open(pickle_file))[1]
    #react_dict=load(open(pickle_file))[2]
    element_index={'O':8,'C':6,'H':1}
    
    for key in react_dict.keys():
        for num,product in enumerate(react_dict[key]):
            num_atoms=len(product.nodes())
            fh=open(path+'/'+key+str(num)+'.xyz','w+')
            fh.write("%i\n \n" % num_atoms)
            fh.close()
            for i,node in enumerate(product.nodes()):
                atom_name=node[0]
                location_0=node[2:-1]
                location_1=location_0.split(',')
                x=float(location_1[0])
                y=float(location_1[1])
                z=float(location_1[2])
                index=element_index[atom_name]
                fh=open(path+'/'+key+str(num)+'.xyz','a+')
                fh.write("%s       %f       %f       %f        %i\n" % (atom_name,x,y,z,index))
                fh.close()
    
    for key in product_dict.keys():
        for num,products in enumerate(product_dict[key]):
            for idx,product in enumerate(products):
                num_atoms=len(product.nodes())
                f=open(path+'/'+key+'_'+str(num)+str(idx)+'.xyz','w+')
                f.write("%i\n \n" % num_atoms)
                f.close()
                for i,node in enumerate(product.nodes()):
                    atom_name=node[0]
                    location_0=node[2:-1]
                    location_1=location_0.split(',')
                    x=float(location_1[0])
                    y=float(location_1[1])
                    z=float(location_1[2])
                    index=element_index[atom_name]
                    f=open(path+'/'+key+'_'+str(num)+str(idx)+'.xyz','a+')
                    f.write("%s       %f       %f       %f        %i\n" % (atom_name,x,y,z,index))
                    f.close()
    return




    




            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
