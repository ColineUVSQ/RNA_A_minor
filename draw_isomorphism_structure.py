'''
Created on 4 févr. 2019

@author: coline

Dessiner les superpositions pour les graphes globaux version matplotlib
'''

from recup_data.sous_graphe_commun_max_version_grands_graphes import nombre_aretes_arcs

import pickle
import networkx as nx
import matplotlib.pyplot as plt
from networkx.classes.function import get_node_attributes
import os
from matplotlib.patches import FancyArrowPatch, Circle
import numpy as np
from recup_data.calcul_sim import calcul_sim_non_cov_sans_motif

''' renvoie vrai si le graphe_commun passe en argument contient le sommet passe en argument
renvoie faux sinon 
le parametre num permet de savoir quel index de noeud il faut regarder (premier index pour le premier graphe ou deuxième index pour le deuxième graphe)'''

def contient_sommet(sommet, graphe_commun, num ):
    for noeud in graphe_commun.nodes() :
        if noeud[num] == sommet :
            return True
    return False

''' renvoie vrai si le graphe_commun passe en argument contient l arete passee en argument
renvoie faux sinon 
le parametre num permet de savoir quel index d'arete il faut regarder (premier index pour le premier graphe ou deuxième index pour le deuxième graphe)'''

def contient_arete(arete, graphe_commun, num):
    for edge in graphe_commun.edges() :
        if edge[0][num] == arete[0] and edge[1][num] == arete[1] or edge[1][num] == arete[0] and edge[0][num] == arete[1]   :
            return True
    return False

'''dessine des aretes courbes'''
def draw_network(G,edges, edges_sans_label, GC, occ_aminor, pos,ax,sg=None):
    
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
#         
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
            if u in occ_aminor and v in occ_aminor :
                color = 'red'
            elif elt[2] == 'CWW' :
                color = 'blue'
            else :
                color = 'black'
            
            for arete_1, arete_2, data_GC in GC.edges(data=True) :
                if (arete_1 == u and arete_2 == v) or (arete_1 == v and arete_2 == u) :
                    if data_GC["type"] == 'COV' and data["label"] == 'B53' :
                        epaisseur = 3
                    elif data_GC["type"] == "CAN" and data["label"] == 'CWW' :
                        epaisseur = 3
                    elif data_GC["type"] == "NON_CAN" and data["label"] != 'CWW' and data["label"] != 'B53':
                        epaisseur = 3
                    else :
                        epaisseur = 1
                else :
                    epaisseur = 1
            
            
            e = FancyArrowPatch(n1.center,n2.center,patchA=n1,patchB=n2,
                                arrowstyle='<|-|>',
                                connectionstyle='arc3,rad=%s'%rad,
                                mutation_scale=10.0,
                                lw=2,
                                alpha=alpha,
                                edgecolor=color,
                                linewidth=epaisseur)
#             if G.nodes[u]["coordonnees"][0] == G.nodes[v]["coordonnees"][0] :
#                 if G.nodes[u]["coordonnees"][1] > G.nodes[v]["coordonnees"][1] :
#                     plt.text(G.nodes[u]["coordonnees"][0]-0.15, G.nodes[v]["coordonnees"][1] + 0.25, data["label"], fontsize=8)
#                 else :
#                     plt.text(G.nodes[u]["coordonnees"][0]-0.15, G.nodes[u]["coordonnees"][1] + 0.25, data["label"], fontsize=8)
#             elif G.nodes[u]["coordonnees"][1] == G.nodes[v]["coordonnees"][1] :
#                 if G.nodes[u]["coordonnees"][0] > G.nodes[v]["coordonnees"][0] :
#                     plt.text(G.nodes[v]["coordonnees"][0]+0.25, G.nodes[u]["coordonnees"][1]-0.15, data["label"], fontsize=8)
#                 else :    
#                     plt.text(G.nodes[u]["coordonnees"][0]+0.25, G.nodes[u]["coordonnees"][1]-0.15, data["label"], fontsize=8)
            seen[(u,v)]=rad
            ax.add_patch(e)
    return e

if __name__ == '__main__':
    
    with open("fichier_comp_grands_graphes_V2.pickle", 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        dico_graphe = mon_depickler.load()
        with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            tab_aminor = mon_depickler.load()
            with open("grands_graphes.pickle", 'rb') as fichier :
                mon_depickler = pickle.Unpickler(fichier)
                dico_graphes = mon_depickler.load()
            
                for comp in dico_graphe.keys() :
                    #comp = (('1FJG', 'A', 48,8), ('5J5B', 'BA', 48,23))
                    fig, axs=plt.subplots(figsize=(10, 12), nrows=1, ncols=2)
                    #fig=plt.figure()
                    compteur = 0
                    columns = 2
                    rows = 1
                    
                    nom = ""
                    GC = dico_graphe[comp].copy()
                    print(comp[0])
                    print(comp[1])
                    compteur_arc_arete_1, compteur_arc_arete_2 = nombre_aretes_arcs(dico_graphes[comp[0]], dico_graphes[comp[1]])
                    recalcul_sim = (((GC.number_of_nodes() + GC.number_of_edges())*(GC.number_of_nodes() + GC.number_of_edges()))/((dico_graphes[comp[0]].number_of_nodes() + compteur_arc_arete_1)*(dico_graphes[comp[1]].number_of_nodes() + compteur_arc_arete_2)))
                    sim_non_cov_sans_motif = calcul_sim_non_cov_sans_motif(dico_graphes[comp[0]], dico_graphes[comp[1]], GC)
                    
                    fig.suptitle('sim = '+str(round(sim_non_cov_sans_motif[0],2)), fontsize=16)
                    for element in comp :
                        
                        #fig.suptitle('sim = '+str(round(dico_sim[elt],2)), fontsize=16)
                        '''on cree un subplot par graphe de la comparaison'''
                        fig.add_subplot(rows, columns, compteur%2+1)
                        #fig.axis('off')
                        #axs[compteur%2].xaxis.set_visible(False)
                        axs[compteur%2].axis("off")
                        axs[compteur%2].set_title(element)
                    
                        print(element)
               
                        #     
                        for occ in tab_aminor :
                            num = (occ["num_PDB"], occ["num_ch"], occ["num_motif"],occ["num_occ"])
                            if num == element :
#                                 print(element)
    #                             for occ in tab_aminor :
    #                                 if occ["num_PDB"] == '5J7L' and occ["num_ch"] == 'DA' and occ["num_motif"] == 197 and occ["num_occ"] == 4 :
    #                                     num = ('5J7L', 'DA', 197, 4)
                                G = dico_graphes[num]
                                nx.set_node_attributes(G, (33,33), "coordonnees")
                                #fichier.write(str(element)+"\n") 
                                #fichier.write(str(G.number_of_nodes())+"\n") 
                                #print(G.edges.data())
                                #print(occ["a_minor"])
                                G.nodes[occ["a_minor"][0]]["coordonnees"] = (0.0,0.5)
                                G.nodes[occ["a_minor"][1]]["coordonnees"] = (2.0,0.5)
                                G.nodes[occ["a_minor"][2]]["coordonnees"] = (0.0,0.0)
                                G.nodes[occ["a_minor"][3]]["coordonnees"] = (2.0,0.0)
                                G.nodes[occ["a_minor"][4]]["coordonnees"] = (3.0,0.5)
                                
                #                 liste_nodes = []
                #                 
                #                 for elt in occ["a_minor"] :
                #                     liste_nodes.append(elt)
                #                 for noeud in G.nodes() :
                #                     if noeud not in liste_nodes :
                #                         liste_nodes.append(noeud)
                                
                                for noeud in G.nodes() :
                                    #voisins = G[noeud] 
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
                                #print(G.nodes.data())
                        #                 for noeud in G.nodes() :
                        #                     fichier.write(str(noeud) + " " + str(G.nodes[noeud])+"\n")
                        #                 for u,v,edata in G.edges(data=True) :
                        #                     fichier.write(str(u)+ "'" +str(v) + " " + str(edata["label"])+"\n")
                                              
                                pos = get_node_attributes(G, 'coordonnees')
                                #print(pos)
                                
#                                 fig = plt.figure(figsize =(5,12))
#                                 fig.suptitle(num, fontsize=16)
                                
                                ''' on regarde si des aretes risquent de se superposer, si c'est le cas, il faudra mettre des aretes courbes '''
                                courbes = []  
                                courbes_avec_label = []       
                                for noeud in G.nodes() :
                                    courbes_x = []
                                    courbes_y = []  
                                    coordonnees_x = []
                                    coordonnees_y = []
                                    for voisin in G[noeud] :
                                        if (noeud, voisin) in G.edges() and (voisin, noeud) in G.edges() :
#                                             print(voisin)
#                                             print(G.nodes[voisin]["coordonnees"][1])
#                                             print(noeud)
#                                             print(G.nodes[noeud]["coordonnees"][1])
                                            if G.nodes[voisin]["coordonnees"][0] == G.nodes[noeud]["coordonnees"][0] :
                                                
                                                for elt in coordonnees_x :
                                                    if G.nodes[voisin]["coordonnees"][0] == elt[1] :
                                                        courbes_x.append((elt[0], voisin))
                                            if G.nodes[voisin]["coordonnees"][1] == G.nodes[noeud]["coordonnees"][1] :
                                                print("ramou")
                                                for elt in coordonnees_y :
                                                    if G.nodes[voisin]["coordonnees"][1] == elt[1] :
                                                        courbes_y.append((elt[0], voisin))
                                        coordonnees_x.append((voisin, G.nodes[voisin]["coordonnees"][0]))
                                        coordonnees_y.append((voisin, G.nodes[voisin]["coordonnees"][1]))
                                    #print(courbes_y)
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
                                #print(courbes)
                                green_edges = []
                                blue_edges = []
                                black_edges = []
                                red_edges = []
                                weights = []
                                #black_edges = [edge for edge in G.edges() if edge not in red_edges]
                                for u,v,edata, in G.edges(data=True) :
                                        if u in occ["a_minor"] and v in occ["a_minor"] :
                                            red_edges.append((u,v))
                                        elif edata["label"] == "B53" :
                                            green_edges.append((u,v))
                                        elif edata["label"] == "CWW" :
                                            blue_edges.append((u,v)) 
                                        else :
                                            black_edges.append((u,v))
                                        
                                        print(str((u,v)))
                                        
                                        ''' on mettra en plus gras les aretes en commun '''
                                        if (u,v) not in courbes and (v,u) not in courbes :
                                            if contient_arete((u,v), GC, compteur%2) :
                                                weights.append(3)
                                                print(True)
                                            else :
                                                weights.append(1)
                                                print(False)
                        
                                edge_labels=dict([((u,v,),d["label"])for u,v,d in G.edges(data=True) if d["label"] != 'B53' and d["label"] != 'CWW' and ((u,v) not in courbes and (v,u) not in courbes)])
                                #print(edge_labels)
                                node_labels=dict([(u,(d["nt"]))for u,d in G.nodes(data=True)])## if d["type"] != None])
                                #node_labels=dict([(u, (u,d["type"], d["poids"])) if d["type"] != None else (u, (u)) for u,d in G.nodes(data=True) ])
                                #node_labels=dict([(u, (u)) for u,d in G.nodes(data=True) ])
                                #print(node_labels)
                                orange_nodes = []
                                pink_nodes = []
                                red_nodes = []
                                '''on mettra en orange les sommets en commun'''
                                for noeud in G.nodes() :
                                    if noeud in occ["a_minor"] :
                                        red_nodes.append(noeud)
                                    elif contient_sommet(noeud, GC, compteur%2) : 
                                        orange_nodes.append(noeud)
                                    else :
                                        pink_nodes.append(noeud)
                                node_colors = ['pink' if node in pink_nodes else 'orange' if node in orange_nodes else 'red' for node in G.nodes()]
                                print(node_colors)
                                #node_colors = ['red' if node in occ["a_minor"] else 'pink' for node in G.nodes()]
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
                                #print(edges_list)
                                #nx.draw_networkx_edge_labels(G,pos, edge_labels = edge_labels, font_size=8)
                                nx.draw_networkx_edges(G,pos, edgelist = edges_list, edge_color=edge_colors,  width=weights)
                                
                                ax=plt.gca()
                                draw_network(G,courbes_avec_label, courbes, GC, occ["a_minor"], pos,ax)
                                print("petit rat")
                                plt.axis('off')
                                nom = nom + "_" +str(element)
                                compteur = compteur+1
    
                                #plt.savefig("graphes_sequence/"+str(num)+".png") # save as png
                                #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
    #                                 plt.show()
    #                                 plt.clf()
    #                                 plt.close()
                    if nom+".png" not in os.listdir("graphes_sequence_V2") :
                        plt.savefig("graphes_sequence_V2/"+nom+".png") # save as png
         
               
                    #plt.show()
                    plt.clf()
                    plt.close()

                