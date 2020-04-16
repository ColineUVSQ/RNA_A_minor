'''
Created on 12 déc. 2018

@author: coline
'''
import pickle
import networkx as nx
import os
import matplotlib.pyplot as plt
with open("graphs_2.92.pickle", 'rb') as fichier_tot :
    mon_depickler = pickle.Unpickler(fichier_tot)
    graphe_tot = mon_depickler.load()
    
#     for fic in os.listdir("graphes_extension/") :
#         if "pickle" in fic :
    fic = "fichier_1FJG_A_48_8.pickle"    
    with open("graphes_extension/"+fic, 'rb') as fichier_graphe :
                mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
                graphe = mon_depickler_graphe.load()
                occ = fic.split("_")[1]
                chaine = fic.split("_")[2]
                num = (occ, chaine)
            
                pos_1 = graphe.nodes[1]["position"][0]
                pos_2 = graphe.nodes[2]["position"][0]
                pos_3 = graphe.nodes[3]["position"][0]
                pos_4 = graphe.nodes[4]["position"][0]
                pos_5 = graphe.nodes[5]["position"][0]
                
                tab_sommets = []
                print(tab_sommets)
                for j in [2,4,1,3] :
                    pos = graphe.nodes[j]["position"][0]
                    tab_sommets.append(pos)
                    for i in range(10) :
                        if j == 1 or j == 4 :
                            tab_sommets.append(pos+i)
                        else :
                            tab_sommets.append(pos-i)
                tab_sommets.append(pos_5)
                    
                print(tab_sommets)
                graphe_mol_subgraph = graphe_tot[num].subgraph(tab_sommets)
                print(graphe_mol_subgraph.nodes.data())
                print(graphe_mol_subgraph.edges.data())
                
                graphe_mol = nx.Graph(graphe_mol_subgraph)
                nx.set_node_attributes(graphe_mol, (33,33), "coordonnees")
                graphe_mol.nodes[pos_1]["coordonnees"] = (0.0,0.5)
                graphe_mol.nodes[pos_2]["coordonnees"] = (2.0,0.5)
                graphe_mol.nodes[pos_3]["coordonnees"] = (0.0,0.0)
                graphe_mol.nodes[pos_4]["coordonnees"] = (2.0,0.0)
                graphe_mol.nodes[pos_5]["coordonnees"] = (3.0,0.5)
                 
                #print(graphe_molraphe_mol.nodes())
                for i in range(len(tab_sommets)) :
                    noeud = tab_sommets[i]
                    #voisins = G[noeud]
                    print(noeud)
                    print("pred")
                    coordonnees_noeud = graphe_mol.nodes[noeud]["coordonnees"]
                    
                    if i < 12 or i > 33 : 
                        print(i)
                        for pred in graphe_tot[num].predecessors(noeud) :
                            print(pred)
#                             if pred not in graphe_mol.nodes() :
#                                 dict_pred = graphe_tot[num].nodes[pred]
#                                 dict_pred.update({'coordonnees' : (33,33)})
#                                 graphe_mol.add_node(pred, **dict_pred)
#                                 dict_edge = graphe_tot[num].edges[pred, noeud]
#                                 graphe_mol.add_edge(pred, noeud, **dict_edge) 
                            if pred in graphe_mol.nodes():
                                if graphe_mol.nodes[pred]["coordonnees"] == (33,33) :
                                    coordonnees = []
                                    for node in graphe_mol.nodes() :
                                        coordonnees.append(graphe_mol.nodes[node]["coordonnees"])
                                    if graphe_mol[pred][noeud]["label"] == "B53" :
                                        if i < 12 : 
                                            if (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                            elif (coordonnees_noeud[0], coordonnees_noeud[1]+1.0) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+1.0)
                                            elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                            elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                            elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                            else :
                                                print("probleme")
                                        else :
                                            if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                            elif (coordonnees_noeud[0], coordonnees_noeud[1]-1.0) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-1.0)
                                            elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                            elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                            elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                            else :
                                                print("probleme")
                                    else :
                                            if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                            elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                            elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                            elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                            elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                            elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                graphe_mol.nodes[pred]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                            else : 
                                                print("probleme") 
                    else :
                        print(i)
                        for succ in graphe_tot[num].successors(noeud) :
#                             if succ not in graphe_mol.nodes() :
#                                 dict_succ = graphe_tot[num].nodes[succ]
#                                 dict_succ.update({'coordonnees' : (33,33)})
#                                 graphe_mol.add_node(succ, **dict_succ) ## etoiles pour desempacter le dictionnaire 
#                                 dict_edge = graphe_tot[num].edges[noeud, succ]
#                                 print(dict_edge)
#                                 graphe_mol.add_edge(noeud, succ, **dict_edge)
                            if succ in graphe_mol.nodes() :
                                if graphe_mol.nodes[succ]["coordonnees"] == (33,33) :
                                    coordonnees = []
                                    for node in graphe_mol.nodes() :
                                        coordonnees.append(graphe_mol.nodes[node]["coordonnees"])
                                        if graphe_mol[noeud][succ]["label"] == "B53" :
                                                if i < 23 :
                                                    if (coordonnees_noeud[0], coordonnees_noeud[1]+0.5) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.5)
                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]+1.0) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+1.0)
                                                    elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                                    elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                    else :
                                                        print("probleme")
                                                else : 
                                                    if (coordonnees_noeud[0], coordonnees_noeud[1]-0.5) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.5)
                                                    elif (coordonnees_noeud[0], coordonnees_noeud[1]-1.0) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-1.0)
                                                    elif (coordonnees_noeud[0]-0.5, coordonnees_noeud[1]) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.5, coordonnees_noeud[1])
                                                    elif (coordonnees_noeud[0]+0.5, coordonnees_noeud[1]) not in coordonnees :
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.5, coordonnees_noeud[1])
                                                    elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                    elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                        graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                    else :
                                                        print("probleme")
                                        else :
                                                if (coordonnees_noeud[0]-0.75, coordonnees_noeud[1]) not in coordonnees :
                                                    graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.75, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0]+0.75, coordonnees_noeud[1]) not in coordonnees :
                                                    graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.75, coordonnees_noeud[1])
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]-0.75) not in coordonnees :
                                                    graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]-0.75)
                                                elif (coordonnees_noeud[0], coordonnees_noeud[1]+0.75) not in coordonnees :
                                                    graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0], coordonnees_noeud[1]+0.75)
                                                elif (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25) not in coordonnees:
                                                    graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]+0.25, coordonnees_noeud[1]+0.25)
                                                elif (coordonnees_noeud[0]-0.25, coordonnees_noeud[1]) not in coordonnees: 
                                                    graphe_mol.nodes[succ]["coordonnees"] = (coordonnees_noeud[0]-0.25, coordonnees_noeud[1])
                                                else : 
                                                    print("probleme") 
                                    
                
                #fichier.write(str(G.nodes.data())+"\n")
                #fichier.write(str(G.edges.data())+"\n")                
                pos = nx.get_node_attributes(graphe_mol, 'coordonnees')
                
                #plt.figure(figsize =(5,12))
                           
                red_edges = [(pos_1,pos_2),(pos_2,pos_1),(pos_1,pos_3),(pos_3,pos_1),(pos_1,pos_5),(pos_5,pos_1),(pos_2,pos_4),(pos_4,pos_2),(pos_3,pos_4),(pos_4,pos_3),(pos_2,pos_5),(pos_5,pos_2)]
                green_edges = []
                blue_edges = []
                black_edges = []
                weights = []
                #black_edges = [edge for edge in G.edges() if edge not in red_edges]
                for u,v,edata, in graphe_mol.edges(data=True) :
                    if (u,v) not in red_edges :
#                             if contient_arete((u,v), GC, compteur%2) :
#                                 orange_edges.append((u,v))
                        if edata["label"] == "B53" :
                            green_edges.append((u,v))
                        elif edata["label"] == "CWW" :
                            blue_edges.append((u,v)) 
                        else :
                            black_edges.append((u,v))
            
                #edge_labels=dict([((u,v,),d["label"])for u,v,d in G.edges(data=True)])
                #print(edge_labels)
                #node_labels=dict([(u,(d["nt"], d["type"]))for u,d in G.nodes(data=True)])## if d["type"] != None])
                node_labels=dict([(u, (u, d["nt"]))for u,d in graphe_mol.nodes(data=True)])
                print(node_labels)
                
                
                
                nx.draw_networkx_nodes(graphe_mol, pos, node_size=150, node_color='pink')
                            #nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                            #           node_color = values, node_size = 500)
                nx.draw_networkx_labels(graphe_mol, pos, labels = node_labels, font_size=5)
                #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', edge_labels = edge_labels)
                #nx.draw_networkx_edges(G, pos, edgelist=blue_edges, edge_color='b', edge_labels = edge_labels)
                #nx.draw_networkx_edges(G, pos, edgelist=green_edges, edge_color='g', edge_labels = edge_labels)
                #nx.draw_networkx_edges(G, pos, edgelist=black_edges, edge_labels = edge_labels)
                
                edge_labels = dict([((u,v), (d["label"])) for u,v,d in graphe_mol.edges(data=True)])
                edge_colors = ['black' if edge in black_edges else 'red' if edge in red_edges else 'blue' if edge in blue_edges else 'green' for edge in graphe_mol.edges()]
                print(black_edges)
                print(red_edges)
                print(edge_colors)
                #nx.draw_networkx_edge_labels(G,pos)
                nx.draw_networkx_edges(graphe_mol,pos,edge_color=edge_colors, edge_label=edge_labels)
                #axs[compteur%2].set_axis_off()
                #plt.savefig("graphes_extension/"+element[:len(element)-7]+".png") # save as png
                #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
                
                
                plt.axis('off')              
                plt.show()
                #plt.savefig("graphes_comparaison/"+nom+".png") # save as png
                plt.clf()
                plt.close()