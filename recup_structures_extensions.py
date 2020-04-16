'''
Created on 29 janv. 2019

@author: coline
'''
import pickle
import networkx as nx
import copy
from recup_data.extension_new import stockage_motif, type_sommet

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

''' renvoie le type du sommet new_position
voisins : l'ensemble de ses voisins
G : le graphe qui est en train d'etre cree
graphe_base : le graphe de structure dont on part pour creer le graphe d'extension
taille_motif : le nombre de nucleotides contenant le motif qu'on veut etendre (5 pour A-minor) '''
def type_sommet(voisins, new_position, G, graphe_base, taille_motif):
    
    for i in range(1,taille_motif) :
        for noeud, data in G.nodes(data=True) :
            if data["motif"] == i and new_position == noeud :
                return i+10
    
    type_sommet_actuel = -1
    if len(voisins) == 0 :
        type_sommet_actuel = 0
    if len(voisins) == 1 : ##Pas de voisin en dehors de la sequence (ou alors il manque la liaison B53)
        for voisin in voisins :
            if graphe_base[new_position][voisin]["label"] == "B53": 
                type_sommet_actuel = 0
            elif graphe_base[new_position][voisin]["label"] == 'CWW' and ((graphe_base.nodes[new_position]["nt"] == 'A' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'A') or (graphe_base.nodes[new_position]["nt"] == 'C' and graphe_base.nodes[voisin]["nt"] == 'G') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'C') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) ::
                type_sommet_actuel = 1
            else :
                type_sommet_actuel = 3
    elif len(voisins) == 2 : ## Un autre voisin en dehors de la sequence
        for voisin in voisins :
            label_voisin = graphe_base[new_position][voisin]["label"]
            if label_voisin != "B53" :
                if label_voisin == 'CWW' and ((graphe_base.nodes[new_position]["nt"] == 'A' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'A') or (graphe_base.nodes[new_position]["nt"] == 'C' and graphe_base.nodes[voisin]["nt"] == 'G') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'C') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                #if label_voisin == 'CWW' :
                ##sommet de type 0
                    type_sommet_actuel = 1       
                    #compteur_tige2 = voisin
                else : ##sommet de type 3

                    type_sommet_actuel = 3
    else : ##Plus d'un autre voisin en dehors de la sequence
        for voisin in voisins : #Recherche d'une liaison can de tige
            label_voisin = graphe_base[new_position][voisin]["label"]
            if label_voisin != "B53" :
                if label_voisin == 'CWW' and ((graphe_base.nodes[new_position]["nt"] == 'A' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'A') or (graphe_base.nodes[new_position]["nt"] == 'C' and graphe_base.nodes[voisin]["nt"] == 'G') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'C') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                #if label_voisin == 'CWW' :
                    type_sommet_actuel = 2
                    #compteur_tige2 = voisin
        if type_sommet_actuel == -1 :
            type_sommet_actuel = 3
            
    return type_sommet_actuel

def recup_structure(taille_ext) :
    with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            tab_aminor = mon_depickler.load()
        
            with open("graphs_2.92.pickle", 'rb') as fichier_tout :
                mon_depickler_graphes = pickle.Unpickler(fichier_tout)
                graphes = mon_depickler_graphes.load()
            
                dico_graphes = {}
                for occ in tab_aminor :
                    est_la = False
                    for elt in liste :
                        if occ["num_PDB"]+"_"+ occ["num_ch"]+"_"+ str(occ["num_motif"])+"_"+ str(occ["num_occ"]) == elt :
                            est_la = True
                            
                    if est_la == False :#and occ["num_PDB"] == '4L81' and occ["num_ch"] == 'A' and occ["num_motif"] == 25 and occ["num_occ"] == 77 :
                        num = (occ["num_PDB"], occ["num_ch"])
                        
                        tab_id = []
                        graphe_a_minor = nx.MultiDiGraph()
                        
                        for j in range(5) :
                            pos = occ["a_minor"][j]
                        
                            if pos not in graphe_a_minor.nodes() :
                                #tab_id.append(pos)
                                attrs = graphes[num].nodes[pos]
                                attrs.update({'chaine' : [j+1], 'motif' : j+1})

                                graphe_a_minor.add_node(pos, **attrs)
#                                 for voisin in graphes[num][pos] :
#                                     
#                                     deja_vu = False
#                                     for voisin_deja_vu in graphe_a_minor.successors(pos) :
#                                         for edge in graphe_a_minor[pos][voisin_deja_vu] :
#                                             if voisin == voisin_deja_vu and (voisin in graphes[num].successors(pos) and graphe_a_minor[pos][voisin_deja_vu][edge]["label"] == graphes[num].edges[pos, voisin]["label"]) or (voisin in graphes[num].predecessors(pos) and graphe_a_minor[pos][voisin_deja_vu][edge]["label"] == graphes[num].edges[voisin, pos]["label"]) :
#                                                 deja_vu = True
#                                     for voisin_deja_vu in graphe_a_minor.predecessors(pos) :
#                                         for edge in graphe_a_minor[voisin_deja_vu][pos] :
#                                             if voisin == voisin_deja_vu and (voisin in graphes[num].successors(pos) and graphe_a_minor[voisin_deja_vu][pos][edge]["label"] == graphes[num].edges[pos, voisin]["label"]) or (voisin in graphes[num].predecessors(pos) and graphe_a_minor[voisin_deja_vu][pos][edge]["label"] == graphes[num].edges[voisin, pos]["label"]) :
#                                                 deja_vu = True
#                                                 
#                                     if deja_vu == False :
#                                         if voisin not in graphe_a_minor.nodes() :
#                                             attrs = graphes[num].nodes[voisin]
#                                             attrs.update({'chaine' : [j+1], 'motif' : False})
#                                             graphe_a_minor.add_node(voisin, **attrs)
#                                             
#                                         if voisin in graphes[num].successors(pos) :
#                                             graphe_a_minor.add_edge(pos, voisin, **graphes[num].edges[pos, voisin])
#                                             if graphes[num].edges[pos, voisin]["label"] != 'B53' :
#                                                 graphe_a_minor.add_edge(voisin, pos, **graphes[num].edges[voisin, pos])
#                                         else :
#                                             graphe_a_minor.add_edge(voisin, pos, **graphes[num].edges[voisin, pos])
#                                             if graphes[num].edges[voisin, pos]["label"] != 'B53' :
#                                                 graphe_a_minor.add_edge(pos, voisin **graphes[num].edges[pos, voisin])
#                                                 
#                                         if j+1 not in graphe_a_minor.nodes[voisin]["chaine"] :
#                                             graphe_a_minor.nodes[voisin]["chaine"].append(j+1)
#                             
                            
                        for j in range(4) :
                                pos = occ["a_minor"][j]
                                print(pos)
                                print(graphes[num].number_of_nodes())
                                for i in range(0,taille_ext) :
                                    print(i)
                                    if i == 7 :
                                            print("petit rat")
                                            print(pos-i)
                                    if pos - i > 0 :
                                        
                                        if pos - i not in graphe_a_minor.nodes() :
                                            #tab_id.append(pos - i)
                                            attrs = graphes[num].nodes[pos-i]
                                            attrs.update({'chaine' : [j+1], 'motif' : False})
                                            graphe_a_minor.add_node(pos-i, **attrs)
                                        
#                                         
#                                         if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
#                                             print(pos-i)
#                                             print(graphes[num][pos-i])
#                                             print(pos+i)
#                                             print(graphes[num][pos+i])
                                        for voisin in graphes[num][pos-i] :
                                            
                                            
                                            deja_vu = False
                                            for voisin_deja_vu in graphe_a_minor.successors(pos-i) :
                                                for edge in graphe_a_minor[pos-i][voisin_deja_vu] :
                                                    if voisin == voisin_deja_vu and ((voisin in graphes[num].successors(pos-i) and graphe_a_minor[pos-i][voisin_deja_vu][edge]["label"] == graphes[num].edges[pos-i, voisin]["label"]) or (voisin in graphes[num].predecessors(pos-i) and graphe_a_minor[pos-i][voisin_deja_vu][edge]["label"] == graphes[num].edges[voisin, pos-i]["label"])) :
                                                        deja_vu = True
                                            for voisin_deja_vu in graphe_a_minor.predecessors(pos-i) :
                                                for edge in graphe_a_minor[voisin_deja_vu][pos-i] :
                                                    if voisin == voisin_deja_vu and ((voisin in graphes[num].successors(pos-i) and graphe_a_minor[voisin_deja_vu][pos-i][edge]["label"] == graphes[num].edges[pos-i, voisin]["label"]) or (voisin in graphes[num].predecessors(pos-i) and graphe_a_minor[voisin_deja_vu][pos-i][edge]["label"] == graphes[num].edges[voisin, pos-i]["label"])) :
                                                        deja_vu = True
                                                        
                                            if deja_vu == False :                                                 
                                                if voisin not in graphe_a_minor.nodes() :
                                                    attrs = graphes[num].nodes[voisin]
                                                    attrs.update({'chaine' : [j+1], 'motif' : False})
                                                    graphe_a_minor.add_node(voisin, **attrs)
                                                
                                                if voisin in graphes[num].successors(pos-i) :
                                                        graphe_a_minor.add_edge(pos-i, voisin, **graphes[num].edges[pos-i, voisin])
                                                        if graphes[num].edges[pos-i, voisin]["label"] != 'B53' :
                                                            graphe_a_minor.add_edge(voisin, pos-i, **graphes[num].edges[voisin, pos-i])
                                                else :
                                                        graphe_a_minor.add_edge(voisin, pos-i, **graphes[num].edges[voisin, pos-i])
                                                        if graphes[num].edges[voisin, pos-i]["label"] != 'B53' :
                                                            graphe_a_minor.add_edge(pos-i, voisin **graphes[num].edges[pos-i, voisin])
                                                            
                                        if (pos-i, pos-i+1) in graphe_a_minor.edges() : 
                                            vu = False
                                            for edge in graphe_a_minor[pos-i][pos-i+1] :
                                                if graphe_a_minor[pos-i][pos-i+1][edge]["label"] == 'B53'  :
                                                    vu = True
                                            if vu == False :
                                                graphe_a_minor.add_edge(pos-i, pos-i+1, label= 'B53', long_range=False)                
#                                                 if j+1 not in graphe_a_minor.nodes[voisin]["chaine"] :
#                                                     graphe_a_minor.nodes[voisin]["chaine"].append(j+1)      
#                                             if ((pos-i, voisin) in graphes[num].edges() and graphes[num].edges[pos-i,voisin]["label"] != 'B53') or ((voisin, pos-i) in graphes[num].edges() and graphes[num].edges[voisin, pos-i]["label"] != 'B53') :
#                                                 if voisin not in tab_id :
#                                                     tab_id.append(voisin)
                                        
                                        if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
                                            print("pos + i : " + str(pos+i))
                                            
                                        
                                            
                                    if pos + i < graphes[num].number_of_nodes()   :  
                                        
                                        if pos + i not in graphe_a_minor.nodes() :
                                            #tab_id.append(pos + i)
                                            attrs = graphes[num].nodes[pos+i]
                                            attrs.update({'chaine' : [j+1], 'motif' : False})
                                            graphe_a_minor.add_node(pos+i, **attrs)
                                         
                                        for voisin in graphes[num][pos+i] :
                                            
                                            if i < 9 or (voisin in graphes[num].successors(pos+i) and graphes[num].edges[pos+i, voisin]["label"] != 'B53') or (voisin in graphes[num].predecessors(pos+i) and graphes[num].edges[voisin, pos+i]["label"] != 'B53') : 
                                                
                                                deja_vu = False
                                                for voisin_deja_vu in graphe_a_minor.successors(pos+i) :
                                                    for edge in graphe_a_minor[pos+i][voisin_deja_vu] :
                                                        if voisin == voisin_deja_vu and ((voisin in graphes[num].successors(pos+i) and graphe_a_minor[pos+i][voisin_deja_vu][edge]["label"] == graphes[num].edges[pos+i, voisin]["label"]) or (voisin in graphes[num].predecessors(pos+i) and graphe_a_minor[pos+i][voisin_deja_vu][edge]["label"] == graphes[num].edges[voisin, pos+i]["label"])) :
                                                            deja_vu = True
                                                            
                                                            if pos+i == 867 : 
                                                                print("petit rat")
                                                                print(voisin)
                                                                print(voisin_deja_vu)
                                                            
                                                for voisin_deja_vu in graphe_a_minor.predecessors(pos+i) :
                                                    for edge in graphe_a_minor[voisin_deja_vu][pos+i] :
                                                        if voisin == voisin_deja_vu and ((voisin in graphes[num].successors(pos+i) and graphe_a_minor[voisin_deja_vu][pos+i][edge]["label"] == graphes[num].edges[pos+i, voisin]["label"]) or (voisin in graphes[num].predecessors(pos+i) and graphe_a_minor[voisin_deja_vu][pos+i][edge]["label"] == graphes[num].edges[voisin, pos+i]["label"])) :
                                                            deja_vu = True
                                                            if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
                                                                print("petit rat")
                                                
                                                if pos+i == 867 : 
                                                        print("ramou")
                                                        print(voisin)
                                                        print(deja_vu)
                                                
                                                if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
                                                    print(deja_vu)            
                                                if deja_vu == False and (i < taille_ext -1 or \
                                                                          ((pos+i,voisin) in graphes[num].edges() and graphes[num].edges[pos+i, voisin]["label"] != 'B53') or ((voisin, pos+i) in graphes[num].edges() and graphes[num].edges[voisin, pos+i]["label"] != 'B53')) :
                                                   
                                                    
                                                    if voisin not in graphe_a_minor.nodes() :
                                                        attrs = graphes[num].nodes[voisin]
                                                        attrs.update({'chaine' : [j+1], 'motif' : False})
                                                        graphe_a_minor.add_node(voisin, **attrs)
                                                        
                                                    if voisin in graphes[num].successors(pos+i) :
                                                            graphe_a_minor.add_edge(pos+i, voisin, **graphes[num].edges[pos+i, voisin])
                                                            if graphes[num].edges[pos+i, voisin]["label"] != 'B53' :
                                                                graphe_a_minor.add_edge(voisin, pos+i, **graphes[num].edges[voisin, pos+i])
                                                    else :
                                                            graphe_a_minor.add_edge(voisin, pos+i, **graphes[num].edges[voisin, pos+i])
                                                            if graphes[num].edges[voisin, pos+i]["label"] != 'B53' :
                                                                graphe_a_minor.add_edge(pos+i, voisin **graphes[num].edges[pos+i, voisin])
#                                                     if j+1 not in graphe_a_minor.nodes[voisin]["chaine"] :
#                                                         graphe_a_minor.nodes[voisin]["chaine"].append(j+1)
                                                
#                                         if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
#                                             print("petit rat")
                                        
                                                
                                        if (pos+i-1, pos+i) in graphe_a_minor.edges() :
                                            vu = False
                                            for edge in graphe_a_minor[pos+i-1][pos+i] :
                                                if graphe_a_minor[pos+i-1][pos+i][edge]["label"] == 'B53' :
                                                    vu = True
                                            if vu == False :
                                                graphe_a_minor.add_edge(pos+i-1, pos+i, label= 'B53', long_range=False)
#                                                 if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
#                                                     print("gros rat")
#                                                     print(pos+i-1)
                                            
#                                             if ((pos+i, voisin) in graphes[num].edges() and graphes[num].edges[pos+i,voisin]["label"] != 'B53') or ((voisin, pos+i) in graphes[num].edges() and graphes[num].edges[voisin, pos+i]["label"] != 'B53') :
#                                                 if voisin not in tab_id :
#                                                     tab_id.append(voisin)
#                         if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
#                             print(occ["a_minor"])
#                             print(tab_id)
#                             print(len(tab_id))
                            
                        #print(tab_id)
                        
                        #graphe_a_minor = graphes[num].subgraph(tab_id)
                        #print(graphe_a_minor.nodes.data())
                        
                        if not graphe_a_minor.has_edge(occ["a_minor"][0], occ["a_minor"][1]) :
                            graphe_a_minor.add_edge(occ["a_minor"][0], occ["a_minor"][1], label="CSS", long_range=True)
                            graphe_a_minor.add_edge(occ["a_minor"][1], occ["a_minor"][0], label="CSS", long_range=True)
                        if not graphe_a_minor.has_edge(occ["a_minor"][2], occ["a_minor"][3]) :
                            graphe_a_minor.add_edge(occ["a_minor"][2], occ["a_minor"][3], label="CSS", long_range=True)
                            graphe_a_minor.add_edge(occ["a_minor"][3], occ["a_minor"][2], label="CSS", long_range=True)
                        if not graphe_a_minor.has_edge(occ["a_minor"][0], occ["a_minor"][4]) :
                            graphe_a_minor.add_edge(occ["a_minor"][0], occ["a_minor"][4], label="TSS", long_range=True)
                            graphe_a_minor.add_edge(occ["a_minor"][4], occ["a_minor"][0], label="TSS", long_range=True)
                            
                        dico_graphes.update({(occ["num_PDB"], occ["num_ch"], occ["num_motif"], occ["num_occ"]) : graphe_a_minor.copy()})
                        
                #print(dico_graphes)     
                #print(graphes)   
                with open("grands_graphes_taille_%s.pickle"%taille_ext, 'wb') as fichier_ecriture :
                    mon_pickler = pickle.Pickler(fichier_ecriture)
                    mon_pickler.dump(dico_graphes)

def recup_structure_new_data(taille_ext) :
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
        
    with open("all_aminor.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            all_aminor = mon_depickler.load()
            
            liste_pb = []
            dico_graphes = {}
            for cle in all_aminor.keys() :
                #print(cle)
                #if cle == "4u4r" :
                    print("petit rat")
                    with open("Graphs/%s.pickle"%cle, 'rb') as fichier_tout :
                        mon_depickler_graphes = pickle.Unpickler(fichier_tout)
                        graphe = mon_depickler_graphes.load()
                        compter_graphe = 1
                        for graph_motif in all_aminor[cle] :
                            #if compter_graphe == 24 :
                            
                        
                                tab_id = []
                                positions_ajoutees = []
                                graphe_a_minor_motif = nx.MultiDiGraph()
                                graphe_a_minor_motif, positions_ajoutees = stockage_motif(graphe_a_minor_motif, graph_motif, positions_ajoutees)
                                
                                print(graphe_a_minor_motif.nodes.data())
                                #exit()
                                graphe_a_minor = nx.MultiDiGraph()
                                for noeud, data in graphe_a_minor_motif.nodes(data=True) :
    
                                    for noeud2, data2 in graph_motif.nodes(data=True) :
    #                                     print(data["position"][0])
    #                                     print(noeud2)
                                        if (data["num_ch"], data["position"][0]) == noeud2 :
                                            data2.update({'chaine' : [noeud],  'motif' : noeud, 'type' : 10+noeud, 'poids' : 1})
                                
                                            graphe_a_minor.add_node(noeud2, **data2)
                                                                    
                                print(graphe_a_minor.nodes.data())
                                
                                
                                
                                for u,v,data in graph_motif.edges(data=True) :
                                    if u in graphe_a_minor.nodes() and v in graphe_a_minor.nodes() :
                                        data.update({'motif' : True})
                                        graphe_a_minor.add_edge(u,v,**data)
                                
                                print(graphe_a_minor.nodes.data())
                                #graphe_a_minor = copy.deepcopy(graph_motif)
                                
                                
                                
    #                             for noeud in graph_motif.nodes() :
    #                                 pos = occ["a_minor"][j]
    #                             
    #                                 if pos not in graphe_a_minor.nodes() :
    #                                     #tab_id.append(pos)
    #                                     attrs = graphes[num].nodes[pos]
    #                                     attrs.update({'chaine' : [j+1], 'motif' : j+1})
    #     
    #                                     graphe_a_minor.add_node(pos, **attrs)
    #                                 for voisin in graphes[num][pos] :
    #                                     
    #                                     deja_vu = False
    #                                     for voisin_deja_vu in graphe_a_minor.successors(pos) :
    #                                         for edge in graphe_a_minor[pos][voisin_deja_vu] :
    #                                             if voisin == voisin_deja_vu and (voisin in graphes[num].successors(pos) and graphe_a_minor[pos][voisin_deja_vu][edge]["label"] == graphes[num].edges[pos, voisin]["label"]) or (voisin in graphes[num].predecessors(pos) and graphe_a_minor[pos][voisin_deja_vu][edge]["label"] == graphes[num].edges[voisin, pos]["label"]) :
    #                                                 deja_vu = True
    #                                     for voisin_deja_vu in graphe_a_minor.predecessors(pos) :
    #                                         for edge in graphe_a_minor[voisin_deja_vu][pos] :
    #                                             if voisin == voisin_deja_vu and (voisin in graphes[num].successors(pos) and graphe_a_minor[voisin_deja_vu][pos][edge]["label"] == graphes[num].edges[pos, voisin]["label"]) or (voisin in graphes[num].predecessors(pos) and graphe_a_minor[voisin_deja_vu][pos][edge]["label"] == graphes[num].edges[voisin, pos]["label"]) :
    #                                                 deja_vu = True
    #                                                 
    #                                     if deja_vu == False :
    #                                         if voisin not in graphe_a_minor.nodes() :
    #                                             attrs = graphes[num].nodes[voisin]
    #                                             attrs.update({'chaine' : [j+1], 'motif' : False})
    #                                             graphe_a_minor.add_node(voisin, **attrs)
    #                                             
    #                                         if voisin in graphes[num].successors(pos) :
    #                                             graphe_a_minor.add_edge(pos, voisin, **graphes[num].edges[pos, voisin])
    #                                             if graphes[num].edges[pos, voisin]["label"] != 'B53' :
    #                                                 graphe_a_minor.add_edge(voisin, pos, **graphes[num].edges[voisin, pos])
    #                                         else :
    #                                             graphe_a_minor.add_edge(voisin, pos, **graphes[num].edges[voisin, pos])
    #                                             if graphes[num].edges[voisin, pos]["label"] != 'B53' :
    #                                                 graphe_a_minor.add_edge(pos, voisin **graphes[num].edges[pos, voisin])
    #                                                 
    #                                         if j+1 not in graphe_a_minor.nodes[voisin]["chaine"] :
    #                                             graphe_a_minor.nodes[voisin]["chaine"].append(j+1)
    #                             
                                compteur_art = 'a'
                                for j in range(1, 5) :
                                        pos = graphe_a_minor_motif.nodes[j]["position"][0]
                                        #print(pos)
                                        #print(graphe.number_of_nodes())
                                        for i in range(0,taille_ext) :
                                            #print(i)
    #                                         if i == 7 :
    #                                                 print("petit rat")
    #                                                 print(pos-i)
                                            if j == 2 or j == 3 :
                                                if pos - i > 0  and (graphe_a_minor_motif.nodes[j]["num_ch"], pos-i) in graphe.nodes() :
                                                    
                                                    num_pos = (graphe_a_minor_motif.nodes[j]["num_ch"], pos-i)
                                                    print("bouh")
                                                    print(num_pos)
                                                    if num_pos not in graphe_a_minor.nodes() :
                                                        #tab_id.append(pos - i)
                                                        print("pap")
                                                        attrs = graphe.nodes[num_pos]
                                                        typ = type_sommet(graphe[num_pos], num_pos, graphe_a_minor, graphe, 4)
                                                        attrs.update({'chaine' : [j], 'motif' : False, 'type' : typ, 'poids' : 1})
                                                        graphe_a_minor.add_node(num_pos, **attrs)
#                                                         print("hihi")
                                                        if typ == 0 :
                                                            data = {'fr3d' : -1, 'nt' : -1, 'real_nt' : -1, 'chaine' : [j], 'motif' : False, 'type' : -1, "poids" : 1}
                                                            graphe_a_minor.add_node((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), **data)
                                                            graphe_a_minor.add_edge(num_pos, (graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), label='0', near=False, motif=False)
                                                            graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), num_pos, label='0', near=False, motif=False)
                                                            compteur_art += str(1)
                                                    else : ## on a deja ajoute le noeud, il faut peut-etre changer son type
                                                        typ = type_sommet(graphe[num_pos], num_pos, graphe_a_minor, graphe, 4)
                                                        ok = False
                                                        print("pip")
                                                        for noeud, data in graphe_a_minor.nodes(data=True) :
                                                            if noeud == num_pos :
                                                                vieux_type = data["type"]
                                                                print(data["type"])
                                                                if data["type"] != typ and data["type"] not in [11,12,13,14,15] :
                                                                    data["type"] = typ
                                                                    ok = True
                                                                    print("hihih")
                                                                    print(j)
                                                                    print(graphe_a_minor[noeud])
                                                                    print(data["chaine"])
                                                                if data["type"] not in [11, 12, 13, 14, 15] :
                                                                    if j not in data["chaine"] or vieux_type in [None, -1] :
                                                                        data["chaine"] = [j]
                                                                        print("hahahaha")
                                                                    elif j not in data["chaine"] :
                                                                        data["chaine"].append(j)
                                                                    
                                                        
                                                        if ok and typ == 0 :
                                                                    data = {'fr3d' : -1, 'nt' : -1, 'real_nt' : -1, 'chaine' : [j], 'motif' : False, 'type' : -1, "poids" : 1}
                                                                    graphe_a_minor.add_node((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), **data)
                                                                    graphe_a_minor.add_edge(num_pos, (graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), label='0', near=False, motif=False)
                                                                    graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), num_pos, label='0', near=False, motif=False)
                                                                    compteur_art += str(1)
                              
            #                                         
            #                                         if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
            #                                             print(pos-i)
            #                                             print(graphes[num][pos-i])
            #                                             print(pos+i)
            #                                             print(graphes[num][pos+i])
                                                    for voisin in graphe[num_pos] :
                                                        #if not (voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') and not (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53') : 
                                                            print("coucou")
                                                            print(num_pos)
                                                        
                                                            deja_vu = False
                                                            for voisin_deja_vu in graphe_a_minor.successors(num_pos) :
                                                                for edge in graphe_a_minor[num_pos][voisin_deja_vu] :
                                                                    if voisin == voisin_deja_vu and ((voisin in graphe.successors(num_pos) and graphe_a_minor[num_pos][voisin_deja_vu][edge]["label"] == graphe.edges[num_pos, voisin]["label"]) or (voisin in graphe.predecessors(num_pos) and graphe_a_minor[num_pos][voisin_deja_vu][edge]["label"] == graphe.edges[voisin, num_pos]["label"])) :
                                                                        deja_vu = True
                                                                        print("rapoulou")
                                                                        print(voisin)
                                                                        print(graphe_a_minor.nodes[voisin]["chaine"])
                                                                        if j not in graphe_a_minor.nodes[voisin]["chaine"] and graphe_a_minor.nodes[voisin]["type"] not in [11, 12, 13, 14, 15]  :
                                                                            if graphe_a_minor.nodes[voisin]["type"] in [None, -1] and ((voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') or (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53')) :
                                                                                graphe_a_minor.nodes[voisin]["chaine"] = [j]
                                                                            elif graphe_a_minor.nodes[voisin]["type"] in [None, -1] :                       
                                                                                graphe_a_minor.nodes[voisin]["chaine"].append(j)                                                         
                                                                                
                                                            for voisin_deja_vu in graphe_a_minor.predecessors(num_pos) :
                                                                for edge in graphe_a_minor[voisin_deja_vu][num_pos] :
                                                                    
                                                                    if voisin == voisin_deja_vu and ((voisin in graphe.successors(num_pos) and graphe_a_minor[voisin_deja_vu][num_pos][edge]["label"] == graphe.edges[num_pos, voisin]["label"]) or (voisin in graphe.predecessors(num_pos) and graphe_a_minor[voisin_deja_vu][num_pos][edge]["label"] == graphe.edges[voisin, num_pos]["label"])) :
                                                                        deja_vu = True
                                                                        print("rapoulou")
                                                                        print(voisin)
                                                                        print(graphe_a_minor.nodes[voisin]["chaine"])
                                                                        if j not in graphe_a_minor.nodes[voisin]["chaine"] and graphe_a_minor.nodes[voisin]["type"] not in [11, 12, 13, 14, 15]  :
                                                                            if graphe_a_minor.nodes[voisin]["type"] in [None, -1] and ((voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') or (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53')) :
                                                                                graphe_a_minor.nodes[voisin]["chaine"] = [j]
                                                                            elif graphe_a_minor.nodes[voisin]["type"] in [None, -1] :
                                                                                graphe_a_minor.nodes[voisin]["chaine"].append(j)
#                                                             print(voisin)
#                                                             print(deja_vu)            
                                                            if deja_vu == False :                                                 
                                                                if voisin not in graphe_a_minor.nodes() :
                                                                    attrs = graphe.nodes[voisin]
                                                                    attrs.update({'chaine' : [j], 'motif' : False, 'type' : None, 'poids' : 1})
                                                                    graphe_a_minor.add_node(voisin, **attrs)
                                                                else :
                                                                    print("rapoulou")
                                                                    print(voisin)
                                                                    print(graphe_a_minor.nodes[voisin]["chaine"])
                                                                    if j not in graphe_a_minor.nodes[voisin]["chaine"] and graphe_a_minor.nodes[voisin]["type"] not in [11, 12, 13, 14, 15]  :
                                                                            if graphe_a_minor.nodes[voisin]["type"] in [None, -1] and ((voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') or (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53')) :
                                                                                graphe_a_minor.nodes[voisin]["chaine"] = [j]
                                                                            elif graphe_a_minor.nodes[voisin]["type"] in [None, -1] :                       
                                                                                graphe_a_minor.nodes[voisin]["chaine"].append(j)
                                                                
                                                                
                                                                if voisin in graphe.successors(num_pos) :
                                                                        if graphe.edges[num_pos, voisin]["label"] == 'CWW' and not ((graphe.nodes[num_pos]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[num_pos]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                                                                            graphe.edges[num_pos, voisin]["label"] = 'CWWn'
                                                                            graphe.edges[voisin, num_pos]["label"] = 'CWWn'
                                                                        graphe.edges[num_pos, voisin].update({'motif' : False})
                                                                        graphe_a_minor.add_edge(num_pos, voisin, **graphe.edges[num_pos, voisin])
                                                                        if graphe.edges[num_pos, voisin]["label"] != 'B53' :
                                                                            graphe.edges[voisin, num_pos].update({'motif' : False})
                                                                            graphe_a_minor.add_edge(voisin, num_pos, **graphe.edges[voisin, num_pos])
                                                                else :
                                                                        if graphe.edges[voisin, num_pos]["label"] == 'CWW' and not ((graphe.nodes[num_pos]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[num_pos]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                                                                            graphe.edges[voisin, num_pos]["label"] = 'CWWn'
                                                                            graphe.edges[num_pos, voisin]["label"] = 'CWWn'
                                                                        graphe.edges[voisin, num_pos].update({'motif' : False})
                                                                        graphe_a_minor.add_edge(voisin, num_pos, **graphe.edges[voisin, num_pos])
                                                                        if graphe.edges[voisin, num_pos]["label"] != 'B53' :
                                                                            graphe.edges[num_pos, voisin].update({'motif' : False})
                                                                            graphe_a_minor.add_edge(num_pos, voisin, **graphe.edges[num_pos, voisin])
                                                                            
                                                    if (num_pos, (num_pos[0], num_pos[1]+1)) in graphe_a_minor.edges() : 
                                                        vu = False
                                                        for edge in graphe_a_minor[num_pos][(num_pos[0], num_pos[1]+1)] :
                                                            if graphe_a_minor[num_pos][(num_pos[0], num_pos[1]+1)][edge]["label"] == 'B53'  :
                                                                vu = True
                                                        if vu == False :
                                                            print("ripili")
                                                            if j not in graphe_a_minor.nodes[num_pos]["chaine"] or graphe_a_minor.nodes[num_pos]["chaine"][0] != j : 
                                                                ch = [j]
                                                                ch.extend([x for x in graphe_a_minor.nodes[num_pos]["chaine"] if x != j])
                                                                graphe_a_minor.nodes[num_pos]["chaine"] = list(ch)
                                                            if j not in graphe_a_minor.nodes[(num_pos[0], num_pos[1]+1)]["chaine"] or graphe_a_minor.nodes[(num_pos[0], num_pos[1]+1)]["chaine"][0] != j: 
                                                                ch = [j]
                                                                ch.extend([x for x in graphe_a_minor.nodes[(num_pos[0], num_pos[1]+1)]["chaine"] if x != j])
                                                                graphe_a_minor.nodes[(num_pos[0], num_pos[1]+1)]["chaine"] = list(ch)
                                                            graphe_a_minor.add_edge(num_pos, (num_pos[0], num_pos[1]+1), label= 'B53', near=False, motif=False)                
            #                                                 if j+1 not in graphe_a_minor.nodes[voisin]["chaine"] :
            #                                                     graphe_a_minor.nodes[voisin]["chaine"].append(j+1)      
            #                                             if ((pos-i, voisin) in graphes[num].edges() and graphes[num].edges[pos-i,voisin]["label"] != 'B53') or ((voisin, pos-i) in graphes[num].edges() and graphes[num].edges[voisin, pos-i]["label"] != 'B53') :
            #                                                 if voisin not in tab_id :
            #                                                     tab_id.append(voisin)
                                                    
        #                                             if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
        #                                                 print("pos + i : " + str(pos+i))
                                                        
                                                
                                            else :        
                                                if pos + i < graphe.number_of_nodes() and (graphe_a_minor_motif.nodes[j]["num_ch"], pos+i) in graphe.nodes()   :  
                                                    num_pos = (graphe_a_minor_motif.nodes[j]["num_ch"], pos+i)
                                                    print(num_pos)
                                                    print(graphe[num_pos])
#                                                     print("ouhou")
#                                                     for noeud, data in graphe_a_minor.nodes(data=True) :
#                                                         print(noeud, data)
                                                    #print(graphe_a_minor.nodes.data())
                                                    if num_pos not in graphe_a_minor.nodes() :
                                                        #tab_id.append(pos + i)
                                                        attrs = graphe.nodes[num_pos]
                                                        typ = type_sommet(graphe[num_pos], num_pos, graphe_a_minor, graphe, 4)
                                                        attrs.update({'chaine' : [j], 'motif' : False, 'type' : typ, 'poids' : 1})
                                                        graphe_a_minor.add_node(num_pos, **attrs)
#                                                         print("ahaha")
                                                        if typ == 0 :
                                                            data = {'fr3d' : -1, 'nt' : -1, 'real_nt' : -1, 'chaine' : [j], 'motif' : False, 'type' : -1, 'poids' : 1}
                                                            graphe_a_minor.add_node((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), **data)
                                                            graphe_a_minor.add_edge(num_pos, (graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), label='0', near=False, motif=False)
                                                            graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), num_pos, label='0', near=False, motif=False)
                                                            compteur_art += str(1)
                                                    else : ## on a deja ajoute le noeud, il faut peut-etre changer son type
                                                        print("ahaha")
                                                        typ = type_sommet(graphe[num_pos], num_pos, graphe_a_minor, graphe, 4)
                                                        ok = False
                                                        print(typ)
                                                        for noeud, data in graphe_a_minor.nodes(data=True) :
                                                            if noeud == num_pos :
                                                                print(data["type"])
                                                                vieux_type = data["type"]
                                                                if data["type"] != typ and data["type"] not in [11, 12, 13, 14, 15]:
                                                                    data["type"] = typ
                                                                    ok = True
                                                                if data["type"] not in [11, 12, 13, 14, 15] :
                                                                    if j not in data["chaine"] or vieux_type in [None, -1] :
                                                                        data["chaine"] = [j]
                                                                    elif j not in data["chaine"] :
                                                                        data["chaine"].append(j)

                                                                
                                                                print(data["chaine"])
                                                        if ok and typ == 0 :
                                                                data = {'fr3d' : -1, 'nt' : -1, 'real_nt' : -1, 'chaine' : [j], 'motif' : False, 'type' : -1, "poids" : 1}
                                                                graphe_a_minor.add_node((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), **data)
                                                                graphe_a_minor.add_edge(num_pos, (graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), label='0', near=False, motif=False)
                                                                graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[j]["num_ch"], compteur_art), num_pos, label='0', near=False, motif=False)
                                                                compteur_art += str(1)
        
                                                     
                                                    for voisin in graphe[num_pos] :
                                                            print("poupou")
                                                            print(graphe_a_minor.nodes())
                                                        #if i < 9 and not (voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') and not (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53') : 
#                                                             print(voisin)
#                                                             print(graphe.edges[num_pos, voisin]["label"])
                                                            deja_vu = False
                                                            for voisin_deja_vu in graphe_a_minor.successors(num_pos) :
                                                                for edge in graphe_a_minor[num_pos][voisin_deja_vu] :
                                                                    if voisin == voisin_deja_vu and ((voisin in graphe.successors(num_pos) and graphe_a_minor[num_pos][voisin_deja_vu][edge]["label"] == graphe.edges[num_pos, voisin]["label"]) or (voisin in graphe.predecessors(num_pos) and graphe_a_minor[num_pos][voisin_deja_vu][edge]["label"] == graphe.edges[voisin, num_pos]["label"])) :
                                                                        deja_vu = True
                                                                        print("ramousnif")
                                                                        print(voisin)
                                                                        #print()
                                                                        if j not in graphe_a_minor.nodes[voisin]["chaine"] and graphe_a_minor.nodes[voisin]["type"] not in [11, 12, 13, 14, 15]  :
                                                                            if graphe_a_minor.nodes[voisin]["type"] in [None, -1] and ((voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') or (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53')) :
                                                                                graphe_a_minor.nodes[voisin]["chaine"] = [j]
                                                                            elif graphe_a_minor.nodes[voisin]["type"] in [None, -1] :                       
                                                                                graphe_a_minor.nodes[voisin]["chaine"].append(j)
                                                                        if pos+i == 867 : 
                                                                            print("petit rat")
                                                                            print(voisin)
                                                                            print(voisin_deja_vu)
                                                                        
                                                            for voisin_deja_vu in graphe_a_minor.predecessors(num_pos) :
                                                                for edge in graphe_a_minor[voisin_deja_vu][num_pos] :
                                                                    if voisin == voisin_deja_vu and ((voisin in graphe.successors(num_pos) and graphe_a_minor[voisin_deja_vu][num_pos][edge]["label"] == graphe.edges[num_pos, voisin]["label"]) or (voisin in graphe.predecessors(num_pos) and graphe_a_minor[voisin_deja_vu][num_pos][edge]["label"] == graphe.edges[voisin, num_pos]["label"])) :
                                                                        deja_vu = True
                                                                        print("ramous")
                                                                        print(voisin)
                                                                        if j not in graphe_a_minor.nodes[voisin]["chaine"] and graphe_a_minor.nodes[voisin]["type"] not in [11, 12, 13, 14, 15]  :
                                                                            if graphe_a_minor.nodes[voisin]["type"] in [None, -1] and ((voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') or (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53')) :
                                                                                graphe_a_minor.nodes[voisin]["chaine"] = [j]
                                                                            elif graphe_a_minor.nodes[voisin]["type"] in [None, -1]:                       
                                                                                graphe_a_minor.nodes[voisin]["chaine"].append(j)
        #                                                                 if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
        #                                                                     print("petit rat")
                                                            
                                                            if pos+i == 867 : 
                                                                    print("ramou")
                                                                    print(voisin)
                                                                    print(deja_vu)
                                                            
        #                                                     if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
        #                                                         print(deja_vu)            
                                                            if deja_vu == False and (i < taille_ext -1 or \
                                                                                      ((num_pos,voisin) in graphe.edges() and graphe.edges[num_pos, voisin]["label"] != 'B53') or ((voisin, num_pos) in graphe.edges() and graphe.edges[voisin, num_pos]["label"] != 'B53')) :
                                                               
                                                                
                                                                if voisin not in graphe_a_minor.nodes() :
                                                                    attrs = graphe.nodes[voisin]
                                                                    attrs.update({'chaine' : [j], 'motif' : False, 'type' : None, 'poids' : 1})
                                                                    graphe_a_minor.add_node(voisin, **attrs)
                                                                else :
                                                                    if j not in graphe_a_minor.nodes[voisin]["chaine"] and graphe_a_minor.nodes[voisin]["type"] not in [11, 12, 13, 14, 15]  :
                                                                        if graphe_a_minor.nodes[voisin]["type"] in [None, -1] and ((voisin in graphe.successors(num_pos) and graphe.edges[num_pos, voisin]["label"] == 'B53') or (voisin in graphe.predecessors(num_pos) and graphe.edges[voisin, num_pos]["label"] == 'B53')) :
                                                                                graphe_a_minor.nodes[voisin]["chaine"] = [j]
                                                                        elif graphe_a_minor.nodes[voisin]["type"] in [None, -1] :                       
                                                                                graphe_a_minor.nodes[voisin]["chaine"].append(j)
                                                                if voisin in graphe.successors(num_pos) :
                                                                        if graphe.edges[num_pos, voisin]["label"] == 'CWW' and not ((graphe.nodes[num_pos]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[num_pos]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                                                                            graphe.edges[num_pos, voisin]["label"] = 'CWWn'
                                                                            graphe.edges[voisin, num_pos]["label"] = 'CWWn'
                                                                        graphe.edges[num_pos, voisin].update({'motif' : False})
                                                                        graphe_a_minor.add_edge(num_pos, voisin, **graphe.edges[num_pos, voisin])
                                                                        if graphe.edges[num_pos, voisin]["label"] != 'B53' :
                                                                            graphe.edges[voisin, num_pos].update({'motif' : False})
                                                                            graphe_a_minor.add_edge(voisin, num_pos, **graphe.edges[voisin, num_pos])
                                                                else :
                                                                        if graphe.edges[voisin, num_pos]["label"] == 'CWW' and not ((graphe.nodes[num_pos]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[num_pos]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[num_pos]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[num_pos]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                                                                            graphe.edges[voisin, num_pos]["label"] = 'CWWn'
                                                                            graphe.edges[num_pos, voisin]["label"] = 'CWWn'
                                                                        graphe.edges[voisin, num_pos].update({'motif' : False})
                                                                        graphe_a_minor.add_edge(voisin, num_pos, **graphe.edges[voisin, num_pos])
                                                                        if graphe.edges[voisin, num_pos]["label"] != 'B53' :
                                                                            graphe.edges[num_pos, voisin].update({'motif' : False})
                                                                            graphe_a_minor.add_edge(num_pos, voisin **graphe.edges[num_pos, voisin])
            #                                                     if j+1 not in graphe_a_minor.nodes[voisin]["chaine"] :
            #                                                         graphe_a_minor.nodes[voisin]["chaine"].append(j+1)
                                                            
            #                                         if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
            #                                             print("petit rat")
                                                    
                                                            
                                                    if ((num_pos[0], num_pos[1]-1), num_pos) in graphe_a_minor.edges() :
#                                                         print("ouhou")
#                                                         print(num_pos)
                                                        vu = False
                                                        for edge in graphe_a_minor[(num_pos[0], num_pos[1]-1)][num_pos] :
                                                            if graphe_a_minor[(num_pos[0], num_pos[1]-1)][num_pos][edge]["label"] == 'B53' :
                                                                print("ehh")
                                                                vu = True
                                                        if vu == False :
#                                                             print("ripili")
#                                                             print(num_pos)
                                                            if j not in graphe_a_minor.nodes[num_pos]["chaine"] or graphe_a_minor.nodes[num_pos]["chaine"][0] != j : 
                                                                ch = [j]
                                                                ch.extend([x for x in graphe_a_minor.nodes[num_pos]["chaine"] if x !=j])
                                                                graphe_a_minor.nodes[num_pos]["chaine"] = list(ch)
                                                            if j not in graphe_a_minor.nodes[(num_pos[0], num_pos[1]-1)]["chaine"] or graphe_a_minor.nodes[(num_pos[0], num_pos[1]-1)]["chaine"][0] != j : 
                                                                ch = [j]
                                                                ch.extend([x for x in graphe_a_minor.nodes[(num_pos[0], num_pos[1]-1)]["chaine"] if x != j])
                                                                graphe_a_minor.nodes[(num_pos[0], num_pos[1]-1)]["chaine"] = list(ch)
                                                            graphe_a_minor.add_edge((num_pos[0], num_pos[1]-1), num_pos, label= 'B53', near=False, motif=False)
                                        
                                
        #                                                 if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
    #                                                     print("gros rat")
    #                                                     print(pos+i-1)
                                                
    #                                             if ((pos+i, voisin) in graphes[num].edges() and graphes[num].edges[pos+i,voisin]["label"] != 'B53') or ((voisin, pos+i) in graphes[num].edges() and graphes[num].edges[voisin, pos+i]["label"] != 'B53') :
    #                                                 if voisin not in tab_id :
    #                                                     tab_id.append(voisin)
    #                         if occ["num_PDB"] == '5DM6' and occ["num_ch"] == 'X' and occ["num_motif"] == 328 and occ["num_occ"] == 2 :
    #                             print(occ["a_minor"])
    #                             print(tab_id)
    #                             print(len(tab_id))
                                
                            #print(tab_id)
                            
                            #graphe_a_minor = graphes[num].subgraph(tab_id)
                            #print(graphe_a_minor.nodes.data())
                                print(graphe_a_minor.nodes())
                                if not graphe_a_minor.has_edge((graphe_a_minor_motif.nodes[1]["num_ch"], graphe_a_minor_motif.nodes[1]["position"][0]), (graphe_a_minor_motif.nodes[2]["num_ch"], graphe_a_minor_motif.nodes[2]["position"][0])) :
                                    graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[1]["num_ch"], graphe_a_minor_motif.nodes[1]["position"][0]), (graphe_a_minor_motif.nodes[2]["num_ch"], graphe_a_minor_motif.nodes[2]["position"][0]), label="CSS", near=False, motif=True)
                                    graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[2]["num_ch"], graphe_a_minor_motif.nodes[2]["position"][0]), (graphe_a_minor_motif.nodes[1]["num_ch"], graphe_a_minor_motif.nodes[1]["position"][0]), label="CSS", near=False, motif=True)
                                if not graphe_a_minor.has_edge((graphe_a_minor_motif.nodes[3]["num_ch"], graphe_a_minor_motif.nodes[3]["position"][0]), (graphe_a_minor_motif.nodes[4]["num_ch"], graphe_a_minor_motif.nodes[4]["position"][0])) :
                                    graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[3]["num_ch"], graphe_a_minor_motif.nodes[3]["position"][0]), (graphe_a_minor_motif.nodes[4]["num_ch"], graphe_a_minor_motif.nodes[4]["position"][0]), label="CSS", near=False, motif=True)
                                    graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[4]["num_ch"], graphe_a_minor_motif.nodes[4]["position"][0]), (graphe_a_minor_motif.nodes[3]["num_ch"], graphe_a_minor_motif.nodes[3]["position"][0]), label="CSS", near=False, motif=True)
                                if not graphe_a_minor.has_edge((graphe_a_minor_motif.nodes[1]["num_ch"], graphe_a_minor_motif.nodes[1]["position"][0]), (graphe_a_minor_motif.nodes[5]["num_ch"], graphe_a_minor_motif.nodes[5]["position"][0])) :
                                    graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[1]["num_ch"], graphe_a_minor_motif.nodes[1]["position"][0]), (graphe_a_minor_motif.nodes[5]["num_ch"], graphe_a_minor_motif.nodes[5]["position"][0]), label="TSS", near=False, motif=True)
                                    graphe_a_minor.add_edge((graphe_a_minor_motif.nodes[5]["num_ch"], graphe_a_minor_motif.nodes[5]["position"][0]), (graphe_a_minor_motif.nodes[1]["num_ch"], graphe_a_minor_motif.nodes[1]["position"][0]), label="TSS", near=False, motif=True)
                                    
                                dico_graphes.update({(cle, compter_graphe) : (graphe_a_minor.copy(), graphe_a_minor_motif.copy())})
                                print(len(dico_graphes))
                                compter_graphe += 1
                                
                #print(dico_graphes)     
                #print(graphes)   
    with open("grands_graphes_new_data_taille_%s_test.pickle"%taille_ext, 'wb') as fichier_ecriture :
                    mon_pickler = pickle.Pickler(fichier_ecriture)
                    mon_pickler.dump(dico_graphes)


def renommage(graphe):
    mapping = {}
    compteur = 6
    for noeud, data in graphe.nodes(data=True) :
        if data["motif"] != False :
            mapping.update({noeud : data["chaine"][0]})
        else :
            mapping.update({noeud : compteur})
            compteur += 1
    print(mapping)
    new_graphe = nx.relabel_nodes(graphe, mapping)
    
    nx.set_node_attributes(new_graphe, [], "position")
    nx.set_node_attributes(new_graphe, -1, "num_ch")
    for cle in mapping :
        new_graphe.nodes[mapping[cle]]["num_ch"] = cle[0]
        new_graphe.nodes[mapping[cle]]["position"] = [cle[1]]
        
    return new_graphe
            

if __name__ == '__main__':
    for i in range(4, 5) :
        recup_structure_new_data(i)

    with open("grands_graphes_new_data_taille_4_test.pickle", 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        dico_graphes = mon_depickler.load()
        
        dico_graphes_renommes = {}
        
        for cle in dico_graphes.keys() :
            new_graphe = renommage(dico_graphes[cle][0])
            dico_graphes_renommes.update({cle : new_graphe})
            
            print(new_graphe.nodes.data())
            #exit()
        
        with open("grands_graphes_new_data_taille_4_renommes.pickle", 'wb') as fichier_2 :
            mon_pickler = pickle.Pickler(fichier_2)
            mon_pickler.dump(dico_graphes_renommes)
        
        print(len(dico_graphes))
        print(list(dico_graphes.keys())[0])
        print(dico_graphes[list(dico_graphes.keys())[0]])
        for u,v,data in dico_graphes[list(dico_graphes.keys())[0]][0].edges(data=True) :
            print(u,v,data)
        #print(dico_graphes[list(dico_graphes.keys())[0]][0].edges.data())
#         print(len(dico_graphes))            
#         print(dico_graphes[('5FDU', '1A', 25, 78)].nodes.data())
#         print(dico_graphes[('5FDU', '1A', 25, 78)].edges.data())
#         
#         print(dico_graphes[('5FDU', '1A', 25, 78)][490])
             
        
    
                    
                        