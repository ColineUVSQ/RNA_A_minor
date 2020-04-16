'''
Created on 18 oct. 2018

@author: Coline Gi
'''
import os
import pickle
import networkx as nx
from networkx.algorithms.coloring.tests.test_coloring import empty_graph
import itertools
from itertools import combinations

#for i in range(len(os.listdir("graphes_extension/fichiers_pickle/"))) :
#    for j in range(i+1, len(os.listdir("graphes_extension/fichiers_pickle/"))) :
#        element1 = os.listdir("graphes_extension/fichiers_pickle/")[i]
#        element2 = os.listdir("graphes_extension/fichiers_pickle/")[j]
element1 = "fichier_1FJG_A_48_8.pickle"
element2 = "fichier_1FJG_A_48_11.pickle"         
with open("graphes_extension/fichiers_pickle/"+element1, 'rb') as fichier1 :
    mon_depickler1 = pickle.Unpickler(fichier1)
    graphe1 = mon_depickler1.load()
    with open("graphes_extension/fichiers_pickle/"+element2, 'rb') as fichier2 :
        mon_depickler2 = pickle.Unpickler(fichier2)
        graphe2 = mon_depickler2.load()
        print(element1)
        print(element2)
        print(graphe1.nodes.data())
        print(graphe1.edges.data())
        print(graphe2.nodes.data())
        print(graphe2.edges.data())
        limite1 = -1
        
        if graphe1.nodes[4]["position"] - graphe1.nodes[2]["position"] < 0 :
            i = 6
            while i < graphe1.number_of_nodes() and( graphe1.nodes[i]["position"] - graphe1.nodes[4]["position"] > 0  or graphe1.nodes[i]["position"] - graphe1.nodes[4]["position"] < -10 ):
                i = i+1
            limite1 = i
        else :
            i = 6
            while i < graphe1.number_of_nodes() and (graphe1.nodes[i]["position"] - graphe1.nodes[4]["position"] < 0  or graphe1.nodes[i]["position"] - graphe1.nodes[4]["position"] > 10 ):
                i = i+1
            limite1 = i
        
        if graphe2.nodes[4]["position"] - graphe2.nodes[2]["position"] < 0 :
            i = 6
            while i < graphe2.number_of_nodes() and (graphe2.nodes[i]["position"] - graphe2.nodes[4]["position"] > 0  or graphe2.nodes[i]["position"] - graphe2.nodes[4]["position"] < -10) :
                i = i+1
            limite2 = i
        else :
            i = 6
            while i < graphe2.number_of_nodes() and (graphe2.nodes[i]["position"] - graphe2.nodes[4]["position"] < 0  or graphe2.nodes[i]["position"] - graphe2.nodes[4]["position"] > 10) :
                i = i+1
            limite2 = i
        
        print(limite1)
        print(limite2)
        result = nx.Graph()
        for i in range(1,6) :
            result.add_node((i,i), type=graphe1.nodes[i]["type"])
                  
        print(result.nodes.data())
        
        #for nbr, datadict in graphe1.adj[2].items():     ##nbr donne le numero du sommet voisin et datadict les infos de l'arete entre les deux   
        #    if nbr != 1 and nbr != 5 and nbr != 4 :
        #        voisin21 = nbr
        #for nbr, datadict in graphe2.adj[2].items():
        #    if nbr != 1 and nbr != 5 and nbr != 4 :
        #        voisin22 = nbr 
        
        
#         graphe1_nodes_vus = []
#         graphe2_nodes_vus = []
#         i = 6
#         j = 6
#         vieux_i = 6
#         vieux_j = 6
#         nouveau_sommet1 = []
#         nouveau_sommet2 = []
#         
#         while i < limite1 :
#             vieux_i = i
#             print(i)
#             print(j)
#             j = i
#             while j < limite2 and i < limite1  :
#                 vieux_j = j
#                 print("i="+str(i))
#                 print("j="+str(j))
#                 #ancien_sommet_principal1 = -1
#                 #ancien_sommet_principal2 = -1
#                 if graphe1.nodes[i]["type"] == graphe2.nodes[j]["type"] and graphe1.nodes[i]["type"] != "None" : ## chercher une bijection
#                     sommet_principal1 = i
#                     sommet_principal2 = j
#                     result.add_node((sommet_principal1,sommet_principal2), type=graphe1.nodes[sommet_principal1]["type"])
#                     graphe1_nodes_vus.append(sommet_principal1)
#                     graphe2_nodes_vus.append(sommet_principal2)
#                     premier = True
#                     while premier == True or (len(nouveau_sommet1) > 0 and len(nouveau_sommet2) > 0) :
#                         if premier == True :
#                             premier = False
#                     #while sommet_principal1 != ancien_sommet_principal1 and sommet_principal2 != ancien_sommet_principal2 :
#                         #ancien_sommet_principal1 = sommet_principal1
#                         #ancien_sommet_principal2 = sommet_principal2
#                         print("ramousnif")
#                         print(sommet_principal1)
#                         print(sommet_principal2)
#                         for nbr1, datadict1 in graphe1.adj[sommet_principal1].items():
#                             for nbr2, datadict2 in graphe2.adj[sommet_principal2].items() :
#                                 print(nbr1)
#                                 print(nbr2)
#                                 print(sommet_principal1)
#                                 print(sommet_principal2)
#                                 print(graphe1.edges[sommet_principal1,nbr1]["weight"])
#                                 print(graphe2.edges[sommet_principal2,nbr2]["weight"])
#                                 if nbr1 not in graphe1_nodes_vus and nbr2 not in graphe2_nodes_vus :
#                                     if datadict1["label"] == datadict2["label"] and abs(graphe1.edges[sommet_principal1,nbr1]["weight"] - graphe2.edges[sommet_principal2,nbr2]["weight"]) <= 2 :
#                                         if (nbr1, nbr2) in result.nodes() :
#                                             #result.add_node((sommet_principal1,sommet_principal2), type=graphe1.nodes[sommet_principal1]["type"])
#                                             result.add_edge((sommet_principal1,sommet_principal2),(nbr1,nbr2), label= datadict1["label"])
#                                             #graphe1_nodes_vus.append(sommet_principal1)
#                                             #graphe2_nodes_vus.append(sommet_principal2)
#                                         else :
#                                             if graphe1.nodes[nbr1]["type"] == graphe2.nodes[nbr2]["type"] :
#                                                 result.add_node((nbr1,nbr2), type=graphe1.nodes[nbr1]["type"])
#                                                 #result.add_node((sommet_principal1,sommet_principal2), type=graphe1.nodes[i]["type"])
#                                                 result.add_edge((sommet_principal1,sommet_principal2),(nbr1,nbr2), label= datadict1["label"])
#                                                 nouveau_sommet1.append(nbr1)
#                                                 nouveau_sommet2.append(nbr2) 
#                                                 
#                                                 graphe1_nodes_vus.append(sommet_principal1)
#                                                 graphe2_nodes_vus.append(sommet_principal2)
#                                                 
#                                                 graphe1_nodes_vus.append(nbr1)
#                                                 graphe2_nodes_vus.append(nbr2)
#                         if len(nouveau_sommet1) != 0 and len(nouveau_sommet2) != 0 :
#                             sommet_principal1 = nouveau_sommet1[0]
#                             sommet_principal2 = nouveau_sommet2[0]
#                             del(nouveau_sommet1[0])
#                             del(nouveau_sommet2[0])
#                     
#                     print(graphe1_nodes_vus)
#                     print(graphe2_nodes_vus)
#                 if vieux_j == j :
#                     j = j+1
#                         
#                 while i in graphe1_nodes_vus :
#                     i = i+1
#                 while j in graphe2_nodes_vus :
#                     j = j+1  
#                 print(i)
#                 print(j)                      
#                 print(result.nodes.data())
#                 print(result.edges.data())
#             if vieux_i == i :
#                 i = i+1
#                     #result.remove_nodes_from(result.nodes())
#                     #for i in range(1,6) :
#                     #    result.add_node((i,i), type=graphe1.nodes[i]["type"])
#         i = limite1
#         j = limite2                       
#         while i <= graphe1.number_of_nodes() :
#             vieux_i = i
#             print(i)
#             print(j)
#             j = limite2
#             while j <= graphe2.number_of_nodes() and i <= graphe1.number_of_nodes() :
#                 vieux_j = j
#                 print("i="+str(i))
#                 print("j="+str(j))
#                 #ancien_sommet_principal1 = -1
#                 #ancien_sommet_principal2 = -1
#                 if graphe1.nodes[i]["type"] == graphe2.nodes[j]["type"] and graphe1.nodes[i]["type"] != "None" : ## chercher une bijection
#                     sommet_principal1 = i
#                     sommet_principal2 = j
#                     result.add_node((sommet_principal1,sommet_principal2), type=graphe1.nodes[sommet_principal1]["type"])
#                     graphe1_nodes_vus.append(sommet_principal1)
#                     graphe2_nodes_vus.append(sommet_principal2)                        
#                     premier = True
#                     while premier == True or (len(nouveau_sommet1) > 0 and len(nouveau_sommet2) > 0) :
#                         if premier == True :
#                             premier = False
#                     #while sommet_principal1 != ancien_sommet_principal1 and sommet_principal2 != ancien_sommet_principal2 :
#                         #ancien_sommet_principal1 = sommet_principal1
#                         #ancien_sommet_principal2 = sommet_principal2
#                         print("ramousnif")
#                         print(sommet_principal1)
#                         print(sommet_principal2)
#                         for nbr1, datadict1 in graphe1.adj[sommet_principal1].items():
#                             for nbr2, datadict2 in graphe2.adj[sommet_principal2].items() :
#                                 print(nbr1)
#                                 print(nbr2)
#                                 print(sommet_principal1)
#                                 print(sommet_principal2)
#                                 print(graphe1.edges[sommet_principal1,nbr1]["weight"])
#                                 print(graphe2.edges[sommet_principal2,nbr2]["weight"])
#                                 if nbr1 not in graphe1_nodes_vus and nbr2 not in graphe2_nodes_vus :
#                                     if datadict1["label"] == datadict2["label"] and abs(graphe1.edges[sommet_principal1,nbr1]["weight"] - graphe2.edges[sommet_principal2,nbr2]["weight"]) <= 2 :
#                                         if (nbr1, nbr2) in result.nodes() :
#                                             result.add_edge((sommet_principal1,sommet_principal2),(nbr1,nbr2), label= datadict1["label"])
#                                             graphe1_nodes_vus.append(sommet_principal1)
#                                             graphe2_nodes_vus.append(sommet_principal2)
#                                         else :
#                                             if graphe1.nodes[nbr1]["type"] == graphe2.nodes[nbr2]["type"] :
#                                                 result.add_node((nbr1,nbr2), type=graphe1.nodes[nbr1]["type"])
#                                                 ##result.add_node((sommet_principal1,sommet_principal2), type=graphe1.nodes[i]["type"])
#                                                 result.add_edge((sommet_principal1,sommet_principal2),(nbr1,nbr2), label= datadict1["label"])
#                                                 nouveau_sommet1.append(nbr1)
#                                                 nouveau_sommet2.append(nbr2) 
#                                                 
#                                                 graphe1_nodes_vus.append(sommet_principal1)
#                                                 graphe2_nodes_vus.append(sommet_principal2)
#                                                 
#                                                 graphe1_nodes_vus.append(nbr1)
#                                                 graphe2_nodes_vus.append(nbr2)
#                         if len(nouveau_sommet1) != 0 and len(nouveau_sommet2) != 0 :
#                             sommet_principal1 = nouveau_sommet1[0]
#                             sommet_principal2 = nouveau_sommet2[0]
#                             del(nouveau_sommet1[0])
#                             del(nouveau_sommet2[0])
#                     
#                     print(graphe1_nodes_vus)
#                     print(graphe2_nodes_vus)
#                 if vieux_j == j :
#                     j = j+1
#                         
#                 while i in graphe1_nodes_vus :
#                     i = i+1
#                 while j in graphe2_nodes_vus :
#                     j = j+1  
#                 print(i)
#                 print(j)                      
#                 print(result.nodes.data())
#                 print(result.edges.data())
#             if vieux_i == i :
#                 i = i+1            
        
        
        with open("fichier_comparaison.txt", 'w') as fichier :
            
            paires = []
            for element in itertools.product(range(6,limite1),range(6,limite2)):
                paires.append(element)
            print(paires)  
            #print(paires)
    
            #print(list(combinations(paires, 4)))
            
            liste = []
            combinaisons = list(combinations(paires,4))
            for elt in combinaisons :
                #print(elt)
                pas_bon = False
                i = 0
                while i < len(elt) :
                    j = i + 1
                    while j < len(elt) :
                        if elt[i][0] == elt[j][0] or elt[i][1] == elt[j][1] :
                            pas_bon = True
                        j = j+1
                    i = i+1
                if pas_bon == False :
                    liste.append(elt) 
            
            print("ramousnif")       
            #print(len(liste))        
            #print(liste)
            
            print("ramousnif")
            deja_vu1 = []
            deja_vu2 = []
            paires_faites = []
            i = 6
            j = 6
            compteur = 0
            premiere_fois = True
            resultat = nx.Graph()
             
            for elt in liste :
                print("compteur: "+str(compteur))
                fichier.write("compteur: "+str(compteur)+"\n")
                if premiere_fois :
                    premiere_fois = False
                else :
                    resultat = empty_graph()
                #print(resultat.nodes.data())
                #print(paires_faites)
                for paire in elt :
                    if graphe1.nodes[paire[0]]["type"] == graphe2.nodes[paire[1]]["type"] : ##and graphe1.nodes[i]["type"] != None :
                            print("ramous")
                            resultat.add_node((paire[0],paire[1]), type=graphe1.nodes[i]["type"])
                compteur = compteur + 1
                print("ramou")
                 
                for noeud1 in resultat.nodes() :
                    for noeud2 in resultat.nodes() :
                        if noeud1 != noeud2 :
                            arete_ok = [0,0]
                            type_arete = [-1,-1]
                            poids = [-1,-1]
                            for nbr, datadict in graphe1.adj[noeud1[0]].items() :
                                #print(datadict)
                                if nbr == noeud2[0] :
                                    arete_ok[0] = 1
                                    if datadict["label"] == "B53" :
                                        type_arete[0] = datadict["label"]
                                    elif datadict["label"] == "CWW" :
                                        type_arete[0] = datadict["label"]
                                    else : 
                                        type_arete[0] = "NON_CAN"
                                    poids[0] = datadict["weight"]
                            for nbr, datadict in graphe2.adj[noeud1[1]].items() :
                                if nbr == noeud2[1] :
                                    arete_ok[1] = 1
                                    if datadict["label"] == "B53" :
                                        type_arete[1] = datadict["label"]
                                    elif datadict["label"] == "CWW" :
                                        type_arete[1] = datadict["label"]
                                    else : 
                                        type_arete[1] = "NON_CAN"
                                    poids[1] = datadict["weight"]
                            if arete_ok[0] == 1 and arete_ok[1] == 1 and type_arete[0] == type_arete[1] and abs(poids[0] - poids[1]) <= 2 :
                                resultat.add_edge(noeud1,noeud2,type = type_arete[1])
                 
                 
                for k in range(1,6) :
                    resultat.add_node((k,k), type=graphe1.nodes[k]["type"])
                     
                resultat.add_edge((1,1),(2,2), label="CSS", weight=0)
                resultat.add_edge((1,1),(3,3), label="B53", weight=0)
                resultat.add_edge((1,1),(5,5), label="TSS", weight=0)
                resultat.add_edge((2,2),(4,4), label="B53", weight=0)
                resultat.add_edge((2,2),(5,5), label="CWW", weight=0)
                resultat.add_edge((3,3),(4,4), label="CSS", weight=0)
                 
                for k in range(1,6) :
                    voisins_1 = []
                    type_arete_1 = []
                    poids_1 = []
                    for nbr, datadict in graphe1.adj[k].items() :
                        for l in range(6,resultat.number_of_nodes()+1) :
                            if nbr == l :
                                voisins_1.append(l)
                                type_arete_1.append(datadict["label"])
                                poids_1.append(datadict["weight"])
                    voisins_2 = []
                    type_arete_2 = []
                    poids_2 = []
                    for nbr, datadict in graphe2.adj[k].items() :
                        for l in range(6,resultat.number_of_nodes()+1) :
                            if nbr == l :
                                voisins_2.append(l)
                                type_arete_2.append(datadict["label"])
                                poids_2.append(datadict["weight"])
                    compteur_1 = 0
                    compteur_2 = 0
                    for elt1 in voisins_1 :
                        for elt2 in voisins_2 :
                            if elt1 == elt2 and type_arete_1[compteur_1] == type_arete_2[compteur_2] and abs(poids_1[compteur_1] - poids_2[compteur_2]) <= 2 :
                                resultat.add_edge((k,k), (elt1,elt2), type=type_arete_2[compteur_2])
                            compteur_2 = compteur_2+1
                        compteur_1 = compteur_1+1
                 
                a_enlever = []
                for noeud in resultat.nodes() :
                    if len(resultat.adj[noeud]) == 0 and (resultat.nodes[noeud]["type"] == 3 or resultat.nodes[noeud]["type"] == 2 or resultat.nodes[noeud]["type"] == None) :
                        a_enlever.append(noeud)
                resultat.remove_nodes_from(a_enlever)
                              
                fichier.write(str(resultat.nodes.data())+"\n")
                fichier.write(str(resultat.edges.data())+"\n")
                    
                         
                               
            
            
                
                
                    
        
        
        
                     