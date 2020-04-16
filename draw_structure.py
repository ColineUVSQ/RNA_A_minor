'''
Created on 10 oct. 2018

@author: Coline Gi

Dessiner les graphes globaux version matplotlib
'''
import pickle
import networkx as nx
import matplotlib.pyplot as plt
from networkx.classes.function import get_node_attributes
import os
from matplotlib.patches import FancyArrowPatch, Circle
import numpy as np
from recup_data.constantes import EXTENSION_PATH_TAILLE

'''pour dessiner des aretes courbes'''
def draw_network(G,edges, edges_sans_label, pos,ax,sg=None):
    
#     print(edges)
#     print(edges_sans_label)
    nodes = []
    for elt in edges :
        if elt[0] not in nodes :
            nodes.append(elt[0])
        if elt[1] not in nodes :
            nodes.append(elt[1])
    
    for n in nodes:
        c=Circle(pos[n],radius=0.02,alpha=0.5)
        ax.add_patch(c)
        G.node[n]['patch']=c
        x,y=pos[n]
    seen={}
    lab = 0
    for u,v,data in G.edges(data=True) :
#         print(u)
#         print(v)
        
        try :
            lab = edges_sans_label.index((u,v))
        except ValueError :
            lab = None
#         print(lab)
#         print(data["label"])
        if lab != None :

            n1=G.node[u]['patch']
            n2=G.node[v]['patch']
            rad=0.1
            if (u,v) in seen:
                rad=seen.get((u,v))
                rad=(rad+np.sign(rad)*0.1)*-1
            alpha=1.0
            #color='k'
            if (u,v) in [(liste_pos_motif[0],liste_pos_motif[1]), (liste_pos_motif[0],liste_pos_motif[4]), (liste_pos_motif[1],liste_pos_motif[0]), (liste_pos_motif[1],liste_pos_motif[3]), (liste_pos_motif[1],liste_pos_motif[4]), (liste_pos_motif[1],liste_pos_motif[0]), (liste_pos_motif[2],liste_pos_motif[3]), (liste_pos_motif[3],liste_pos_motif[2]), (liste_pos_motif[4],liste_pos_motif[0]), (liste_pos_motif[4],liste_pos_motif[1])] :
                color = 'red'
            elif data["label"] == 'CWW' :
                color = 'blue'
            else :
                color = 'black'
            
            e = FancyArrowPatch(n1.center,n2.center,patchA=n1,patchB=n2,
                                arrowstyle='<|-|>',
                                connectionstyle='arc3,rad=%s'%rad,
                                mutation_scale=10.0,
                                lw=2,
                                alpha=alpha,
                                edgecolor=color, 
                                label=data["label"])
            if G.nodes[u]["coordonnees"][0] == G.nodes[v]["coordonnees"][0] :
                if G.nodes[u]["coordonnees"][1] > G.nodes[v]["coordonnees"][1] :
                    plt.text(G.nodes[u]["coordonnees"][0]-0.15, G.nodes[v]["coordonnees"][1] + 0.25, data["label"], fontsize=8)
                else :
                    plt.text(G.nodes[u]["coordonnees"][0]-0.15, G.nodes[u]["coordonnees"][1] + 0.25, data["label"], fontsize=8)
            elif G.nodes[u]["coordonnees"][1] == G.nodes[v]["coordonnees"][1] :
                if G.nodes[u]["coordonnees"][0] > G.nodes[v]["coordonnees"][0] :
                    plt.text(G.nodes[v]["coordonnees"][0]+0.25, G.nodes[u]["coordonnees"][1]-0.15, data["label"], fontsize=8)
                else :    
                    plt.text(G.nodes[u]["coordonnees"][0]+0.25, G.nodes[u]["coordonnees"][1]-0.15, data["label"], fontsize=8)
            seen[(u,v)]=rad
            ax.add_patch(e)
    return e

def draw_structure(G, liste_pos_motif, num, rep, nom_fichier_coord):
    print(G.nodes.data())
    nx.set_node_attributes(G, (33,33), "coordonnees")
    #fichier.write(str(element)+"\n") 
    #fichier.write(str(G.number_of_nodes())+"\n") 
    #print(G.edges.data())
    #print(occ["a_minor"])
    G.nodes[liste_pos_motif[0]]["coordonnees"] = (0.0,0.5)
    G.nodes[liste_pos_motif[1]]["coordonnees"] = (2.0,0.5)
    G.nodes[liste_pos_motif[2]]["coordonnees"] = (0.0,0.0)
    G.nodes[liste_pos_motif[3]]["coordonnees"] = (2.0,0.0)
    G.nodes[liste_pos_motif[4]]["coordonnees"] = (2.75,0.5)
    
#                 liste_nodes = []
#                 
#                 for elt in occ["a_minor"] :
#                     liste_nodes.append(elt)
#                 for noeud in G.nodes() :
#                     if noeud not in liste_nodes :
#                         liste_nodes.append(noeud)
    
    for noeud in G.nodes() :
        #voisins = G[noeud] 
                if len(G[noeud]) == 0 :
                    print(noeud)
                
                coordonnees_noeud = G.nodes[noeud]["coordonnees"]
                for pred in G.predecessors(noeud) :
                    if G.nodes[pred]["coordonnees"] == (33,33) :
                        coordonnees = []
                        for node in G.nodes() :
                            coordonnees.append(G.nodes[node]["coordonnees"])
                        for edge in G[pred][noeud] :
                            if G[pred][noeud][edge]["label"] == "B53" :
                                if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                    G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                else :
                                    #fichier.write("probleme\n")
                                    print("probleme")
                            else :
#                                                         print(G.nodes[noeud])
                                if 1 in G.nodes[pred]["chaine"]:
                                    if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme") 
                                    
                                elif 2 in G.nodes[pred]["chaine"]  :
                                    if (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme") 
                                elif 3 in G.nodes[pred]["chaine"] :
                                    if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme") 
                                else : # chaine = 4
                                    if (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme")                           
                for succ in G.successors(noeud) :
                    if G.nodes[succ]["coordonnees"] == (33,33) :
                        coordonnees = []
                        for node in G.nodes() :
                            coordonnees.append(G.nodes[node]["coordonnees"])
                        for edge in G[noeud][succ] :
                            if G[noeud][succ][edge]["label"] == "B53" :
                                if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                    G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                else :
                                    #fichier.write("probleme\n")
                                    print("probleme")
                            else :
                                if 1 in G.nodes[succ]["chaine"]:
                                    if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme") 
                                    
                                elif 2 in G.nodes[succ]["chaine"]  :
                                    if (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme") 
                                elif 3 in G.nodes[succ]["chaine"] :
                                    if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme") 
                                else : # chaine = 4
                                    if (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                        G.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                    else : 
                                        #fichier.write("probleme\n")
                                        print("probleme")             
#                             print(G.nodes.data())
#                 for noeud in G.nodes() :
#                     fichier.write(str(noeud) + " " + str(G.nodes[noeud])+"\n")
#                 for u,v,edata in G.edges(data=True) :
#                     fichier.write(str(u)+ "'" +str(v) + " " + str(edata["label"])+"\n")
                  
    pos = get_node_attributes(G, 'coordonnees')
    #print(pos)
    
    fig = plt.figure(figsize =(5,12))
    fig.suptitle(num, fontsize=16)
    
    
    '''si plusieurs aretes se superposent, on utilise des courbes '''
    courbes = []  
    courbes_avec_label = []       
    for noeud in G.nodes() :
        courbes_x = []
        courbes_y = []  
        coordonnees_x = []
        coordonnees_y = []
        for voisin in G[noeud] :
            if (noeud, voisin) in G.edges() and (voisin, noeud) in G.edges() :
#                                         print(voisin)
#                                         print(G.nodes[voisin]["coordonnees"][1])
#                                         print(noeud)
#                                         print(G.nodes[noeud]["coordonnees"][1])
                if G.nodes[voisin]["coordonnees"][0] == G.nodes[noeud]["coordonnees"][0] :
                    
                    for elt in coordonnees_x :
                        if G.nodes[voisin]["coordonnees"][0] == elt[1] :
                            courbes_x.append((elt[0], voisin))
                if G.nodes[voisin]["coordonnees"][1] == G.nodes[noeud]["coordonnees"][1] :
#                                             print("ramou")
                    for elt in coordonnees_y :
                        if G.nodes[voisin]["coordonnees"][1] == elt[1] :
                            courbes_y.append((elt[0], voisin))
            coordonnees_x.append((voisin, G.nodes[voisin]["coordonnees"][0]))
            coordonnees_y.append((voisin, G.nodes[voisin]["coordonnees"][1]))
#                                 print(courbes_y)
        if len(courbes_x) > 0 :
            for elt in courbes_x :
                if (G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] < G.nodes[elt[1]]["coordonnees"][1]) or (G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[0]]["coordonnees"][1] and G.nodes[noeud]["coordonnees"][1] > G.nodes[elt[1]]["coordonnees"][1]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) > abs(G.nodes[elt[1]]["coordonnees"][1] - G.nodes[noeud]["coordonnees"][1]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                            for edge in G[noeud][elt[0]] :
                                courbes_avec_label.append((noeud, elt[0], G[noeud][elt[0]][edge]["label"]))
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
                            for edge in G[noeud][elt[1]] :
                                courbes_avec_label.append((noeud, elt[1], G[noeud][elt[1]][edge]["label"]))
        if len(courbes_y) > 0 :
            for elt in courbes_y :
                if (G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] < G.nodes[elt[1]]["coordonnees"][0]) or (G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[0]]["coordonnees"][0] and G.nodes[noeud]["coordonnees"][0] > G.nodes[elt[1]]["coordonnees"][0]) :
                    if abs(G.nodes[elt[0]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) > abs(G.nodes[elt[1]]["coordonnees"][0] - G.nodes[noeud]["coordonnees"][0]) :
                        if (noeud, elt[0]) not in courbes and (elt[0], noeud) not in courbes  :
                            courbes.append((noeud, elt[0]))
                            for edge in G[noeud][elt[0]] :
                                courbes_avec_label.append((noeud, elt[0], G[noeud][elt[0]][edge]["label"]))
                    else :
                        if (noeud, elt[1]) not in courbes and (elt[1], noeud) not in courbes :
                            courbes.append((noeud, elt[1]))
                            for edge in G[noeud][elt[1]] :
                                courbes_avec_label.append((noeud, elt[1], G[noeud][elt[1]][edge]["label"]))
#                             print(courbes)
    
    with open(rep+nom_fichier_coord, "wb") as fichier_sortie :
        mon_pickler = pickle.Pickler(fichier_sortie)
        mon_pickler.dump(G)
    
    print(G.nodes.data())
    
    green_edges = []
    blue_edges = []
    black_edges = []
    red_edges = []
    #black_edges = [edge for edge in G.edges() if edge not in red_edges]
    for u,v,edata, in G.edges(data=True) :
            if (u,v) in [(liste_pos_motif[0],liste_pos_motif[1]), (liste_pos_motif[0],liste_pos_motif[4]), (liste_pos_motif[1],liste_pos_motif[0]), (liste_pos_motif[1],liste_pos_motif[3]), (liste_pos_motif[1],liste_pos_motif[4]), (liste_pos_motif[1],liste_pos_motif[0]), (liste_pos_motif[2],liste_pos_motif[3]), (liste_pos_motif[2], liste_pos_motif[0]), (liste_pos_motif[3],liste_pos_motif[2]), (liste_pos_motif[4],liste_pos_motif[0]), (liste_pos_motif[4],liste_pos_motif[1])] :
                red_edges.append((u,v))
            elif edata["label"] == "B53" :
                green_edges.append((u,v))
            elif edata["label"] == "CWW" :
                blue_edges.append((u,v)) 
                print(u,v,edata)
            else :
                black_edges.append((u,v))

    edge_labels=dict([((u,v,),d["motif"])for u,v,d in G.edges(data=True) if d["label"] != 'B53' and d["label"] != 'CWW' and ((u,v) not in courbes and (v,u) not in courbes)])
    #print(edge_labels)
    node_labels=dict([(u,(u, d["chaine"]))for u,d in G.nodes(data=True)])## if d["type"] != None])
    #node_labels=dict([(u, (u,d["type"], d["poids"])) if d["type"] != None else (u, (u)) for u,d in G.nodes(data=True) ])
    #node_labels=dict([(u, (u)) for u,d in G.nodes(data=True) ])
    #print(node_labels)
    node_colors = ['red' if node in liste_pos_motif else 'pink' for node in G.nodes()]
    nx.draw_networkx_nodes(G, pos, node_size=150, node_color=node_colors)
                #nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                #           node_color = values, node_size = 500)
    nx.draw_networkx_labels(G, pos, labels = node_labels, font_size = 8)
    #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=blue_edges, edge_color='b', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=green_edges, edge_color='g', edge_labels = edge_labels)
    #nx.draw_networkx_edges(G, pos, edgelist=black_edges, edge_labels = edge_labels)
    
    edge_colors = ['black' if (u,v) in black_edges else 'blue' if (u,v) in blue_edges else 'green' if (u,v) in green_edges else 'red' for u,v in G.edges() if (u,v) not in courbes and (v,u) not in courbes]
#         print(black_edges)
#         print(red_edges)
#         print(edge_colors)
    #nx.draw_networkx_edge_labels(G,pos)
    edges_list = [(u,v,)for u,v in G.edges() if ((u,v) not in courbes and (v,u) not in courbes)] 
#                             print(edges_list)
    nx.draw_networkx_edge_labels(G,pos, edge_labels = edge_labels, font_size=8)
    nx.draw_networkx_edges(G,pos, edgelist = edges_list, edge_color=edge_colors)
    
    ax=plt.gca()
    draw_network(G,courbes_avec_label, courbes, pos,ax)
    print("petit rat")
    plt.axis('off')
    
    for u,v,data in G.edges(data=True) :
        print(u,v,data)
    
    #plt.savefig(rep+str(num)+".png") # save as png
    #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
    plt.show()
    plt.clf()
    plt.close()
    
    
if __name__ == '__main__':
    for i in range(4, 5) :
        with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            tab_aminor = mon_depickler.load()
            with open("grands_graphes_new_data_taille_%d_test.pickle"%i, 'rb') as fichier :
                mon_depickler = pickle.Unpickler(fichier)
                dico_graphes = mon_depickler.load()
                #     
                
                types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]
                liste_tout = []
                with open("resolutions.pickle", 'rb') as fichier_pickle :
                            mon_depickler = pickle.Unpickler(fichier_pickle)
                            resolutions = mon_depickler.load()
                for elt in types_arn :
                           
                                   
                    with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                                    mon_depickler = pickle.Unpickler(fichier_sortie)
                                    liste_a_garder_noms = mon_depickler.load()    
                                             
                                    for element in liste_a_garder_noms :
                                        if element == ('4w2f', 16) :
                                            print(element)
                                            print(elt)
                                        if resolutions[element[0]] <= 3.0 :
                                            if element in liste_tout :
                                                print(element)
                                                print(elt)
                                            if element not in liste_tout :
                                                liste_tout.append(element)
                         
                print(liste_tout)
                
                os.makedirs("Graphes_globaux/graphes_sequence/taille_%d_new_data"%i, exist_ok=True)

                for elt in liste_tout :
                    if elt in dico_graphes.keys() :
                        if elt == ('4y4o', 2) :
#                         ok = False
#                         for noeud in dico_graphes[elt][0].nodes() :
#                             liste_voisins = []
#                             for voisin in dico_graphes[elt][0][noeud] :
#                                 if voisin in liste_voisins :
#                                     ok = True
#                                 liste_voisins.append(voisin)
#                         
#                                     
#                         if ok :
                            G = dico_graphes[elt]
                            liste_pos_motif = []
                            for noeud, data in G[1].nodes(data=True) :
                                print(noeud)
                                liste_pos_motif.append((data["num_ch"], data["position"][0]))
    #                                                 
                            draw_structure(G[0], liste_pos_motif, elt[0]+"_"+str(elt[1]), "Graphes_globaux/graphes_sequence/taille_%d_new_data/"%i, "graphe_global_%s_%s_avec_coord_test.pickle"%(elt[0], elt[1]))
                        #exit()                    
                                    
                                
                                                    
        
            
                