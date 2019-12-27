'''
Created on 9 oct. 2018

@author: Coline Gi

Construction des extensions à partir des graphes globaux de structure
(version CaRNAval)
'''
import pickle
import networkx as nx
from recup_data.environnement_liaisons_non_can import nb_canoniques_short_range,\
    nb_non_canoniques_short_range, nb_non_cov
from recup_data.constantes import EXTENSION_PATH_TAILLE


liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']

''' ajout d'un sommet au graphe d'extension et d'aretes si necessaire  '''
def ajout_sommet(G, compteur, compteur_tige, position, typ, poids, label, positions_ajoutees, int_tige, valeur_debut, *args, **kwargs):
    long_range = kwargs.get("long_range", False)
    
    voisin_chaine =  kwargs.get("voisin_chaine", -1)
    voisin_voisin =  kwargs.get("voisin_voisin", -1)
    voisin_none = False
    a_voisin_voisin = False
#     if len(pos_liaisons) > 0 :
#         if typ == 1 :
#             if len(pos_liaisons[0]) > 1 : ## il faut decouper le sommet en deux voire plus car il est lie a plusieurs sommets de l autre chaine
#                 for elt in pos_liaisons[0] :
#                     G.add_node(compteur, position=position, type=typ, poids=G.nodes[elt]["poids"], chaine=[valeur_debut], pos_liaisons=pos_liaisons[1])
#                     
#                     if valeur_debut not in G.nodes[elt]["pos_liaisons"] and G.nodes[elt]["type"] == 1 :
#                         G.nodes[elt]["pos_liaisons"].append(valeur_debut)
#                         
#                     G.add_edge(compteur, elt, label="CWW", long_range=False)   
#                     G.add_edge(elt, compteur, label="CWW", long_range=False)
#                     
#                     compteur += 1
#                 compteur = compteur-1
#             else :
#                 G.add_node(compteur, position=position, type=typ, poids=poids, chaine=[valeur_debut], pos_liaisons=pos_liaisons[1]) 
#                 
#                 if valeur_debut not in G.nodes[pos_liaisons[0]]["pos_liaisons"] and G.nodes[pos_liaisons[0]]["type"] == 1 :
#                     G.nodes[elt]["pos_liaisons"].append(valeur_debut)
#                     G.add_edge(compteur, elt, label="CWW", long_range=False)
#                     G.add_edge(elt, compteur, label="CWW", long_range=False)
#         else :
#             G.add_node(compteur, position=position, type=typ, poids=poids, chaine=[valeur_debut], pos_liaisons=[])
#     else : 
#         G.add_node(compteur, position=position, type=typ, poids=poids, chaine=[valeur_debut], pos_liaisons=pos_liaisons)
    
    
    
    G.add_node(compteur, position=position, type=typ, poids=poids, chaine=[valeur_debut])
    
    if typ == 1 : ## il faut peut etre rajouter des aretes si on a deja mis le sommet de la deuxieme chaine
#         if compteur == 14 :
#             print("ramou")
#             print(voisin_chaine)
#         if compteur == 24 :
#             print("ramous")
#             print(G.nodes.data())
#             print(G.edges.data())
#             print(voisin_chaine)

        noeud_existe = False
        for noeud, data in G.nodes(data=True) :
            if voisin_chaine >= data["position"][0] and voisin_chaine <= data["position"][1] :
                noeud_existe = noeud
        if noeud_existe != False :
            G.add_edge(compteur, noeud_existe, label="CWW", long_range=False)
            G.add_edge(noeud_existe, compteur, label="CWW", long_range=False)
        
        
        
        ##cas ou une liaison can avec un sommet None
        if voisin_voisin != False : 
            noeud_existe = False
            for noeud, data in G.nodes(data=True) :
                if voisin_voisin >= data["position"][0] and voisin_voisin <= data["position"][1] :
                    noeud_existe = noeud
            if noeud_existe :
                G.add_edge(compteur, noeud_existe, label="CWW", long_range=False)
                G.add_edge(noeud_existe, compteur, label="CWW", long_range=False)
            else :
                G.add_node(compteur+1, position=(voisin_voisin, voisin_voisin), type=None, poids = 1, chaine=[valeur_debut])
                G.add_edge(compteur, compteur+1, label="CWW", long_range=False)
                G.add_edge(compteur+1, compteur, label="CWW", long_range=False)
                
                voisin_none = voisin_voisin
            a_voisin_voisin = True
        
    if label == "B53" :
        if G.nodes[compteur]["position"][0] < G.nodes[compteur_tige]["position"][0] :
            if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[compteur]["position"][1]) == 1 :
                if (compteur, compteur_tige) not in G.edges() :
                    G.add_edge(compteur, compteur_tige, label="B53", long_range=False)
        else :
            if abs(G.nodes[compteur]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                if (compteur_tige, compteur) not in G.edges() :
                    G.add_edge(compteur_tige, compteur, label="B53", long_range=False)
    else :
        if len(label) == 3 :
            label_inverse = label[0]+label[2]+label[1]
        else :
            label_inverse = label[0]+label[2]+label[1] + label[3]
        G.add_edge(compteur_tige, compteur, label=label, long_range=long_range)
        G.add_edge(compteur, compteur_tige, label=label_inverse, long_range=long_range)
    for i in range(position[0], position[1]+1) :
        positions_ajoutees.append(i)
    
    print("ramousnif")
    print(voisin_voisin)
    print(voisin_none)
    return voisin_none, a_voisin_voisin
# def tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees) :
#     valeur_debut = compteur_tige
#    
# #     if valeur_debut == 2 or valeur_debut == 1:
# #         new_position = G.node[valeur_debut]["position"][0]-i*int_tige
# #     else :
# #         new_position = G.node[valeur_debut]["position"][0]+i*int_tige 
# 
#     new_position = G.node[valeur_debut]["position"][0]
#     
#     while valeur_debut < 5 and new_position > 0 and new_position < graphes[nom_cle].number_of_nodes() :
#         #print(i)
#         voisins = graphes[nom_cle][new_position]
#         if nom_cle == ('4V88', 'A6') :
#             print(new_position)
#             print("voisins")
#         for voisin in voisins :
# #             if nom_cle == ('4V88', 'A6') :
# #                 print(voisin)
# #                 print(positions_ajoutees)
#             if voisin not in positions_ajoutees :
#                 
#                 G.add_node(compteur, position = (voisin, voisin), type=None, poids=1, chaine=[valeur_debut])
#                 G.add_edge(compteur, valeur_debut, label=graphes[nom_cle].edges[new_position,voisin]["label"])
#                 G.add_edge(valeur_debut, compteur, label=graphes[nom_cle].edges[new_position,voisin]["label"])
#                 compteur = compteur + 1
#                 positions_ajoutees.append(voisin)
#             else :
#                 if nom_cle == ('4V88', 'A6') :
#                     print(voisin)
#                     print(G.edges())
#                 for u, pos in G.nodes(data="position") :
#                     if pos[0] <= voisin and pos[1] >= voisin  :
#                         if (u,valeur_debut) not in G.edges() and (valeur_debut, u) not in G.edges():
#                             G.add_edge(u, valeur_debut, label=graphes[nom_cle].edges[new_position,voisin]["label"])
#                             G.add_edge(valeur_debut, u, label=graphes[nom_cle].edges[new_position,voisin]["label"])
#                         
#                             if valeur_debut not in G.nodes[u]["chaine"] :
#                                 G.nodes[u]["chaine"].append(valeur_debut)
#         
#         valeur_debut += 1               
#         new_position = G.node[valeur_debut]["position"][0]
#         
#     
#     return G

''' renvoie le type du sommet new_position
voisins : l'ensemble de ses voisins
G : le graphe d'extension qui est en train d'etre cree
graphe_base : le graphe de structure dont on part pour creer le graphe d'extension '''
def type_sommet(voisins, new_position, G, graphe_base):
    
    if new_position == G.nodes[1]["position"][0] :
        return 11 
    if new_position == G.nodes[2]["position"][0] :
        return 12 
    if new_position == G.nodes[3]["position"][0] :
        return 13
    if new_position == G.nodes[4]["position"][0] :
        return 14
    if new_position == G.nodes[5]["position"][0] :
        return 15 
    
    type_sommet_actuel = -1
    if len(voisins) == 0 :
        type_sommet_actuel = 0
    if len(voisins) == 1 : ##Pas de voisin en dehors de la sequence
        for voisin in voisins :
            if graphe_base[new_position][voisin]["label"] == "B53": 
                type_sommet_actuel = 0
            else :
                type_sommet_actuel = 3
    elif len(voisins) == 2 : ## Un autre voisin en dehors de la sequence
        for voisin in voisins :
            label_voisin = graphe_base[new_position][voisin]["label"]
            if label_voisin != "B53" :
                if label_voisin == 'CWW' and graphe_base[new_position][voisin]["long_range"] == False and ((graphe_base.nodes[new_position]["nt"] == 'A' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'A') or (graphe_base.nodes[new_position]["nt"] == 'C' and graphe_base.nodes[voisin]["nt"] == 'G') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'C') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                ##sommet de type 0
                    type_sommet_actuel = 1       
                    #compteur_tige2 = voisin
                else : ##sommet de type 3

                    type_sommet_actuel = 3
    else : ##Plus d'un autre voisin en dehors de la sequence
        for voisin in voisins : #Recherche d'une liaison can de tige
            label_voisin = graphe_base[new_position][voisin]["label"]
            if label_voisin != "B53" :
                if label_voisin == 'CWW' and graphe_base[new_position][voisin]["long_range"] == False and ((graphe_base.nodes[new_position]["nt"] == 'A' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'A') or (graphe_base.nodes[new_position]["nt"] == 'C' and graphe_base.nodes[voisin]["nt"] == 'G') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'C') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                    type_sommet_actuel = 2
                    #compteur_tige2 = voisin
        if type_sommet_actuel == -1 :
            type_sommet_actuel = 3
            
    return type_sommet_actuel

''' extension d'une chaine du motif
G : graphe d'extension en cours de creation
graphes : ensemble des graphes globaux stockes dans un dictionnaire
nom_cle : nom de la cle du graphe global courant
compteur : numero du dernier sommet ajoute a G
compteur_tige : numero du dernier sommet de la chaine ajoute a G
positions_ajoutees : liste des positions du graphe global deja ajoutees a G
int_tige : 1 ou -1 en fonction du sens dans lequel on doit se deplacer sur la chaine (sens de la sequence ou sens inverse)
taille_ext : taille de l'extension
positions_chaines : positions dans le graphe global courant des sommets appartenant a l'une des chaines'''
def extension_tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines):
    #print(compteur)
    valeur_debut = compteur_tige
    i = 0
    poids_sommet = 1
    type_sommet_prec = None
    type_sommet_actuel = None 
    type_sommet_voisin = -1
    type_sommet_voisin_prec = -1
   
    if valeur_debut == 2 or valeur_debut == 1:
        new_position = G.node[valeur_debut]["position"][0]-i*int_tige
    else :
        new_position = G.node[valeur_debut]["position"][0]+i*int_tige 

    #new_position = G.node[valeur_debut]["position"][0]                    
#     if nom_cle == ('4V88', 'A6') and occ["num_motif"] == 17 and occ["num_occ"] == 55 :
#                         print(G.edges.data())
    
    position_prec = new_position
    position_prec_groupe = new_position
    
    voisin_chaine = -1
    voisin_chaine_prec = -1
    chaine_sommet_voisin = -1
    chaine_sommet_voisin_prec = -1
    voisin_dans_chaine = False
    voisin_voisin = False
    voisin_voisin_prec = False
    
    
    while i < taille_ext and new_position > 0 and new_position <= graphes[nom_cle].number_of_nodes() :
        if i == 0 or new_position not in positions_ajoutees[:4] :## pour ne pas revenir sur le motif
            
            voisins = graphes[nom_cle][new_position]
 
            type_sommet_actuel = type_sommet(voisins, new_position, G, graphes[nom_cle])

            if type_sommet_actuel == 1 :
                for voisin in voisins :
                    label_voisin = graphes[nom_cle][new_position][voisin]["label"]
                    if label_voisin != "B53" :
                        if label_voisin == 'CWW' and graphes[nom_cle][new_position][voisin]["long_range"] == False : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                            voisins_du_voisin = non_cov_graphe_base(voisin, graphes[nom_cle])
                            
                            voisin_dans_chaine = False
                            voisin_voisin_prec = voisin_voisin
                            #voisin_voisin = False
                            for k in range(4) :
      
                                if voisin in positions_chaines[k] : ## pour savoir si le voisin appartient a une des chaines ou pas 
                                    
                                    type_sommet_voisin_prec = type_sommet_voisin
                                    type_sommet_voisin = type_sommet(graphes[nom_cle][voisin],voisin, G, graphes[nom_cle])
                                    chaine_sommet_voisin_prec = chaine_sommet_voisin
                                    chaine_sommet_voisin = k

#                                     if compteur == 6 :
#                                         print("ramou")
#                                         print(voisin_chaine)
#                                         print(chaine_sommet_voisin)
#                                         print(chaine_sommet_voisin_prec)
                                    
                                    voisin_dans_chaine = True
                                
                                for elt in voisins_du_voisin :

                                    if elt != new_position :
                                        if elt in positions_chaines[k] :
                                            voisin_voisin = voisin
                            #print(voisin_voisin)
                            
                            voisin_chaine_prec = voisin_chaine
                            voisin_chaine = voisin
                            if new_position == 1506 :
                                print("gros rat")
                                print(voisin_chaine_prec)
                            
                            if voisin_dans_chaine == False :
                                chaine_sommet_voisin_prec = chaine_sommet_voisin
                                chaine_sommet_voisin = -1
                                type_sommet_voisin_prec = type_sommet_voisin
                                type_sommet_voisin = -1
                            else :
                                voisin_voisin = voisin_voisin_prec              
            if new_position not in positions_ajoutees :
                
                if type_sommet_actuel == 0 or type_sommet_actuel == 1 :  
                    
                    if type_sommet_actuel == type_sommet_prec \
                     and (type_sommet_actuel == 0 or  \
                          (type_sommet_actuel == 1 and \
                            (voisin_voisin == False and \
                             (type_sommet_voisin_prec == type_sommet_voisin and (type_sommet_voisin == 1 or type_sommet_voisin == 0)) and \
                             (chaine_sommet_voisin == chaine_sommet_voisin_prec) or (voisin_voisin == False) and \
                             (type_sommet_voisin == -1 and type_sommet_voisin_prec == -1) ) and \
                            ((voisin_dans_chaine == False and (abs(voisin_chaine-voisin_chaine_prec) < 10 or voisin_chaine == -1 or voisin_chaine_prec == -1)) or \
                            (voisin_dans_chaine == True and (abs(voisin_chaine-voisin_chaine_prec) == 1 or voisin_chaine == -1 or voisin_chaine_prec == -1))))): 
                        poids_sommet += 1
                        
    #                 elif type_sommet_actuel == 1 and type_sommet_voisin_prec != type_sommet_voisin and type_sommet_voisin != None  : ## il faut ajouter le sommet precedent
    #                     
    #                     
                    else :
                        
                        if (type_sommet_prec == 1 or type_sommet_prec == 0) :  ## on va ajouter le sommet precedent car il est de type 0 ou 1   
                            
                            
                            deja_vu = False
                            for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                                if k in positions_ajoutees :
                                    deja_vu = True
                                    
                            #if position_prec not in positions_ajoutees :
                            if deja_vu == False :
    #                             if compteur == 38 :
    #                                 print("ramousnif")
    #                                 print(chaine)
                                
                                if position_prec - position_prec_groupe <= 0 :
                                    if type_sommet_actuel != 1 :
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                    else : 
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                    
                                else : 
                                    if type_sommet_actuel != 1 :
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                    else :
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                                if new_position == 1506 :
                                    print("ramous")
                                    print(G.edges.data())
                                    print(compteur)
                                    print(voisin_chaine_prec)
                                    print(voisin_voisin)
                                    print(voisin_none)
                                compteur_tige = compteur
                                if voisin_none != False:
                                    positions_ajoutees.append(voisin_none)
                                    compteur = compteur+2
                                else :
                                    compteur = compteur+1
                                    
                                if a_voisin_voisin :
                                    voisin_voisin = False
                                
                            else : ## le sommet d'avant existe deja mais il est peut etre en groupe maintenant
                                
                                if poids_sommet != 1 :
                                    #print("petit rat")
                                    num_sommet_prec = -1
                                    #print("ramou")
                                    for noeud in G.nodes() :
                                        for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                                            if  k <= G.nodes[noeud]["position"][1] and k >= G.nodes[noeud]["position"][0] :
                                                num_sommet_prec = noeud ## retrouver le numero du sommet
                                                if valeur_debut not in G.nodes[noeud]["chaine"] :
                                                    G.nodes[noeud]["chaine"].append(valeur_debut)
                                                    
                                                for voisin in G[noeud] :
                                                    if G.nodes[voisin]["type"] == None :
                                                        G.nodes[voisin]["chaine"].append(valeur_debut)
                                    if G.nodes[num_sommet_prec]["type"] == None :
                                        G.nodes[num_sommet_prec]["type"] = type_sommet_prec
                                        del(G.nodes[noeud]["chaine"][:])
                                        G.nodes[noeud]["chaine"].append(valeur_debut)
                                    if num_sommet_prec == -1 :
                                        print("probleme6") 
                                        
                                    if position_prec - position_prec_groupe <= 0 :
                                        G.nodes[num_sommet_prec]["position"] = (min(position_prec, G.nodes[num_sommet_prec]["position"][0]), max(position_prec_groupe, G.nodes[num_sommet_prec]["position"][1]))
                                    else :
                                        G.nodes[num_sommet_prec]["position"] = (min(position_prec_groupe, G.nodes[num_sommet_prec]["position"][0]), max(position_prec, G.nodes[num_sommet_prec]["position"][1]))
                                    G.nodes[num_sommet_prec]["poids"] = G.nodes[num_sommet_prec]["position"][1] - G.nodes[num_sommet_prec]["position"][0] + 1
                                    
    #                                 if type_sommet_prec == 1 :
    #                                     
    #                                     if len(chaine) > 0 :
    #                                         for elt in chaine[1] :
    #                                             if elt not in G.nodes[num_sommet_prec]["pos_liaisons"] :
    #                                                 G.nodes[num_sommet_prec]["pos_liaisons"].append(elt)
    #                                         for elt in chaine[0] :
    #                                             if valeur_debut not in G.nodes[elt]["pos_liaisons"] and G.nodes[elt]["type"] == 1  :
    #                                                 G.nodes[elt]["pos_liaisons"].append(valeur_debut)
    #                                                 G.add_edge(num_sommet_prec, elt, label="CWW", long_range=False)
    #                                                 G.add_edge(elt, num_sommet_prec, label="CWW", long_range=False)
    #                                     else : 
    #                                         G.nodes[num_sommet_prec]["pos_liaisons"].append(chaine)
                                    
                                    if G.nodes[num_sommet_prec]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                                        if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_sommet_prec]["position"][1]) == 1 :
                                            if (num_sommet_prec, compteur_tige) not in G.edges() :
                                                G.add_edge(num_sommet_prec, compteur_tige, label="B53", long_range=False)
                                    else :
                                        if abs(G.nodes[num_sommet_prec]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                            if (compteur_tige, num_sommet_prec) not in G.edges() :
                                                G.add_edge(compteur_tige, num_sommet_prec, label="B53", long_range=False)
                                    
                                    compteur_tige = num_sommet_prec 
                        position_prec_groupe = new_position    
                        type_sommet_prec = type_sommet_actuel
                        poids_sommet = 1  
                        
                        
                        
                         
                            
                else :
                    if poids_sommet != 1 or (type_sommet_prec == 0 or type_sommet_prec == 1 and poids_sommet == 1): ## on va ajouter le sommet precedent car il est de type 0 ou 1
                        
                        deja_vu = False
                        for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                            if k in positions_ajoutees :
                                deja_vu = True
                              
                        #if position_prec not in positions_ajoutees : ## le sommet d'avant sur la sequence de type 0 ou 1 n'etait pas deja vu
                        if deja_vu == False :
                            if position_prec - position_prec_groupe <= 0 :
                                if type_sommet_actuel != 1   :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                else : 
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                
                            else : 
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                else :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                            
                            compteur_tige = compteur
                            if voisin_none != False:

                                positions_ajoutees.append(voisin_none)
                                compteur = compteur+2
                            else :
                                compteur = compteur+1
                                
                                
                            if a_voisin_voisin :
                                voisin_voisin = False
                        else :
                            if poids_sommet != 1 :
                                #print("petit rat")
                                num_sommet_prec = -1
                                #print("ramou")
                                for noeud in G.nodes() :
                                    for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                                        if  k <= G.nodes[noeud]["position"][1] and k >= G.nodes[noeud]["position"][0] :
                                            num_sommet_prec = noeud ## retrouver le numero du sommet
                                            if valeur_debut not in G.nodes[noeud]["chaine"] :
                                                G.nodes[noeud]["chaine"].append(valeur_debut)
                                                
                                            for voisin in G[noeud] :
                                                if G.nodes[voisin]["type"] == None :
                                                    G.nodes[voisin]["chaine"].append(valeur_debut)
                                if G.nodes[num_sommet_prec]["type"] == None :
                                    G.nodes[num_sommet_prec]["type"] = type_sommet_prec
                                    del(G.nodes[noeud]["chaine"][:])
                                    G.nodes[noeud]["chaine"].append(valeur_debut)
                                if num_sommet_prec == -1 :
                                    print("probleme6") 
                                    
                                if position_prec - position_prec_groupe <= 0 :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec, G.nodes[num_sommet_prec]["position"][0]), max(position_prec_groupe, G.nodes[num_sommet_prec]["position"][1]))
                                else :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec_groupe, G.nodes[num_sommet_prec]["position"][0]), max(position_prec, G.nodes[num_sommet_prec]["position"][1]))
                                G.nodes[num_sommet_prec]["poids"] = G.nodes[num_sommet_prec]["position"][1] - G.nodes[num_sommet_prec]["position"][0] + 1
                                
    #                                 if type_sommet_prec == 1 :
    #                                     
    #                                     if len(chaine) > 0 :
    #                                         for elt in chaine[1] :
    #                                             if elt not in G.nodes[num_sommet_prec]["pos_liaisons"] :
    #                                                 G.nodes[num_sommet_prec]["pos_liaisons"].append(elt)
    #                                         for elt in chaine[0] :
    #                                             if valeur_debut not in G.nodes[elt]["pos_liaisons"] and G.nodes[elt]["type"] == 1  :
    #                                                 G.nodes[elt]["pos_liaisons"].append(valeur_debut)
    #                                                 G.add_edge(num_sommet_prec, elt, label="CWW", long_range=False)
    #                                                 G.add_edge(elt, num_sommet_prec, label="CWW", long_range=False)
    #                                     else : 
    #                                         G.nodes[num_sommet_prec]["pos_liaisons"].append(chaine)
                                
                                if G.nodes[num_sommet_prec]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                                    if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_sommet_prec]["position"][1]) == 1 :
                                        if (num_sommet_prec, compteur_tige) not in G.edges() :
                                            G.add_edge(num_sommet_prec, compteur_tige, label="B53", long_range=False)
                                else :
                                    if abs(G.nodes[num_sommet_prec]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                        if (compteur_tige, num_sommet_prec) not in G.edges() :
                                            G.add_edge(compteur_tige, num_sommet_prec, label="B53", long_range=False)
                                
                                compteur_tige = num_sommet_prec
    #                     else :
    #                         num_voisin_seq = -1
    #                         #print(position_prec)
    #                         for noeud in G.nodes() :
    #                             if position_prec <= G.nodes[noeud]["position"][1] and position_prec >= G.nodes[noeud]["position"][0] :
    #                                 num_voisin_seq = noeud
    #                         if num_voisin_seq == -1 :
    #                             print("probleme")
    #                         #G.add_edge(num_voisin_seq, compteur_tige, label="B53")
                    poids_sommet = 1
                    ajout_sommet(G, compteur, compteur_tige, (new_position, new_position), type_sommet_actuel, 1, "B53", positions_ajoutees, int_tige, valeur_debut)
                    
                    
                    compteur_tige = compteur
                    compteur = compteur+1
                    
                    ## ajout de chacun des voisins du sommet actuel car il est de type 2 ou 3
                    for voisin in voisins : 
                        label_voisin = graphes[nom_cle][new_position][voisin]["label"]
                        if label_voisin != "B53" :
                            if voisin not in positions_ajoutees :
                                if label_voisin == 'CWW' and not ((graphes[nom_cle].nodes[new_position]["nt"] == 'A' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'A') or (graphes[nom_cle].nodes[new_position]["nt"] == 'C' and graphes[nom_cle].nodes[voisin]["nt"] == 'G') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'C') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'G'))  : 
                                    label_voisin = 'CWWn'
                                ajout_sommet(G, compteur, compteur_tige, (voisin,voisin), None, 1, label_voisin, positions_ajoutees, int_tige, valeur_debut, long_range = graphes[nom_cle][new_position][voisin]["long_range"] )
                                compteur = compteur+1
                            else :
                                num_sommet = -1
                                #print("ramou")
                                for noeud in G.nodes() :
                                    if  voisin <= G.nodes[noeud]["position"][1] and voisin >= G.nodes[noeud]["position"][0] :
                                        num_sommet = noeud ## retrouver le numero du sommet 
                                if num_sommet == -1 :
                                    print("probleme2")
                                deja_fait = False
                                for v in G[num_sommet] :
                                    for edge in G[num_sommet][v] :
                                        if v == compteur_tige and G[num_sommet][v][edge]["label"] != "B53" : 
                                            deja_fait = True 
                                if deja_fait == False :
                                    if label_voisin == 'CWW' and not ((graphes[nom_cle].nodes[new_position]["nt"] == 'A' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'A') or (graphes[nom_cle].nodes[new_position]["nt"] == 'C' and graphes[nom_cle].nodes[voisin]["nt"] == 'G') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'C') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'G'))  : 
                                        label_voisin = 'CWWn'
                                    if len(label_voisin) == 3 :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1]
                                    else :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1] + label_voisin[3]
                                    G.add_edge(num_sommet, compteur_tige, label=label_voisin_inverse, long_range = graphes[nom_cle][new_position][voisin]["long_range"])
                                    G.add_edge(compteur_tige, num_sommet, label=label_voisin, long_range = graphes[nom_cle][new_position][voisin]["long_range"])
                    type_sommet_prec = type_sommet_actuel 
    
            else : ## le sommet actuel existe deja
                
                num_sommet = -1
                #print(G.nodes.data())
                #print(positions_ajoutees)
                #print(new_position)
                for noeud in G.nodes() :
                    if new_position <= G.nodes[noeud]["position"][1] and new_position >= G.nodes[noeud]["position"][0]  :
                        num_sommet = noeud ## retrouver le numero du sommet 
                        if valeur_debut not in G.nodes[noeud]["chaine"] :
                            G.nodes[noeud]["chaine"].append(valeur_debut)
                            
                        for voisin in G[noeud] :
                            if G.nodes[voisin]["type"] == None :
                                G.nodes[voisin]["chaine"].append(valeur_debut)
                                
                        if G.nodes[noeud]["type"] == None :
                            G.nodes[noeud]["type"] = type_sommet_actuel
                            del(G.nodes[noeud]["chaine"][:])
                            G.nodes[noeud]["chaine"].append(valeur_debut)
                        #else :
                        #    print("probleme")
    #                     if G.nodes[noeud]["type"] == 0 :
    #                         print(new_position)
    #                         print(position_prec_groupe)
    #                         print(num_sommet)
    #                         print(poids_sommet)
    #                         
    #                         
    #                         print(G.nodes.data())
    #                         print("ramousnif")
                if num_sommet == -1 :
                    print("probleme5")
                    #print(G.nodes.data())
                    #print(new_position)
                
                                                          
                #print(num_sommet)
                #print(compteur_tige)
                #print(compteur)
                if type_sommet_actuel == 0 or type_sommet_actuel == 1 :  
                    if type_sommet_actuel == type_sommet_prec \
                    and (type_sommet_actuel == 0 or  \
                          (type_sommet_actuel == 1 and \
                            ((voisin_voisin == False) and \
                              (type_sommet_voisin_prec == type_sommet_voisin and (type_sommet_voisin == 1 or type_sommet_voisin == 0)) and \
                             (chaine_sommet_voisin == chaine_sommet_voisin_prec) or (voisin_voisin == False) and \
                             (type_sommet_voisin == -1 and type_sommet_voisin_prec == -1) ) and \
                            ((voisin_dans_chaine == False and (abs(voisin_chaine-voisin_chaine_prec) < 10 or voisin_chaine == -1 or voisin_chaine_prec == -1)) or \
                            (voisin_dans_chaine == True and (abs(voisin_chaine-voisin_chaine_prec) == 1 or voisin_chaine == -1 or voisin_chaine_prec == -1))))): 
                        poids_sommet += 1
                    else :
                        
                        deja_vu = False
                        for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                            if k in positions_ajoutees :
                                deja_vu = True
                        #if position_prec not in positions_ajoutees :
                        if deja_vu == False :
                            if position_prec - position_prec_groupe <= 0 :
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                else : 
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                
                            else : 
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                else :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                            
                            compteur_tige = compteur
                            if voisin_none != False:

                                positions_ajoutees.append(voisin_none)
                                compteur = compteur+2
                            else :
                                compteur = compteur+1
                                
                                
                            if a_voisin_voisin :
                                voisin_voisin = False
                        else : ## le sommet d'avant existe deja mais il est peut etre en groupe maintenant
    #                         print("petit rat")
    #                         print(new_position)
                            
                            if poids_sommet != 1 :
                                num_sommet_prec = -1
                                    #print("ramou")
                                for noeud in G.nodes() :
                                    for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                                        if  k <= G.nodes[noeud]["position"][1] and k >= G.nodes[noeud]["position"][0] :
                                            num_sommet_prec = noeud ## retrouver le numero du sommet 
                                            if valeur_debut not in G.nodes[noeud]["chaine"] :
                                                G.nodes[noeud]["chaine"].append(valeur_debut)
                                                
                                            for voisin in G[noeud] :
                                                if G.nodes[voisin]["type"] == None :
                                                    G.nodes[voisin]["chaine"].append(valeur_debut)
                                if G.nodes[num_sommet_prec]["type"] == None :
                                    G.nodes[num_sommet_prec]["type"] = type_sommet_prec
                                    del(G.nodes[noeud]["chaine"][:])
                                    G.nodes[noeud]["chaine"].append(valeur_debut)
                                if num_sommet_prec == -1 :
                                    print("probleme7") 
                                        
                                if position_prec - position_prec_groupe <= 0 :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec, G.nodes[num_sommet_prec]["position"][0]), max(position_prec_groupe, G.nodes[num_sommet_prec]["position"][1]))
                                else :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec_groupe, G.nodes[num_sommet_prec]["position"][0]), max(position_prec, G.nodes[num_sommet_prec]["position"][1]))
                                G.nodes[num_sommet_prec]["poids"] = G.nodes[num_sommet_prec]["position"][1] - G.nodes[num_sommet_prec]["position"][0] + 1
                                
#                                 if nom_cle == ('2XD0', 'V') :
#                                     print('ramou')
#                                     print(compteur_tige)
#                                     print(num_sommet)
#                                     print(num_sommet_prec)
#                                     print(G.edges.data())
                                
                                if compteur_tige != num_sommet_prec :
                                    if G.nodes[num_sommet_prec]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                                        if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_sommet_prec]["position"][1]) == 1 :
                                            if (num_sommet_prec, compteur_tige) not in G.edges() :
                                                G.add_edge(num_sommet_prec, compteur_tige, label="B53", long_range=False)
                                    else :
                                        if abs(G.nodes[num_sommet_prec]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                            if (compteur_tige, num_sommet_prec) not in G.edges() :
                                                G.add_edge(compteur_tige, num_sommet_prec, label="B53", long_range=False)
                                
#                                 if nom_cle == ('2XD0', 'V') :
#                                     print('ramou')
#                                     print(compteur_tige)
#                                     print(num_sommet)
#                                     print(num_sommet_prec)
#                                     print(G.edges.data())
                                
                                compteur_tige = num_sommet_prec
                                
    #                             if type_sommet_prec == 1 :
    #                                 if len(chaine) > 0 :
    #                                     for elt in chaine[1] :
    #                                         if elt not in G.nodes[num_sommet_prec]["pos_liaisons"] :
    #                                             G.nodes[num_sommet_prec]["pos_liaisons"].append(elt)
    #                                     for elt in chaine[0] :
    #                                         if valeur_debut not in G.nodes[elt]["pos_liaisons"] and G.nodes[elt]["type"] == 1  :
    #                                             G.nodes[elt]["pos_liaisons"].append(valeur_debut)
    #                                             G.add_edge(num_sommet_prec, elt, label="CWW", long_range=False)
    #                                             G.add_edge(elt, num_sommet_prec, label="CWW", long_range=False)
    #                                 else : 
    #                                     G.nodes[num_sommet_prec]["pos_liaisons"].append(chaine)
                                
                        if G.nodes[num_sommet]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                            if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_sommet]["position"][1]) == 1 :
                                if (num_sommet, compteur_tige) not in G.edges() :
                                    G.add_edge(num_sommet, compteur_tige, label="B53", long_range=False)
                        else :
                            if abs(G.nodes[num_sommet]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                if (compteur_tige, num_sommet) not in G.edges() :
                                    G.add_edge(compteur_tige, num_sommet, label="B53", long_range=False)
                        compteur_tige = num_sommet
                            
    #                     if position_prec - position_prec_groupe <= 0 :
    #                         G.nodes[num_sommet]["position"] = (position_prec, position_prec_groupe)
    #                     else :
    #                         G.nodes[num_sommet]["position"] = (position_prec_groupe, position_prec)
                            
                        position_prec_groupe = new_position                       
                        type_sommet_prec = type_sommet_actuel
                        poids_sommet = 1
    
                        
                else :
    
                    if poids_sommet != 1 or (type_sommet_prec == 0 or type_sommet_prec == 1 and poids_sommet == 1):                        
                        deja_vu = False
                        for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                            if k in positions_ajoutees :
                                deja_vu = True
                        
                        #if position_prec not in positions_ajoutees : ## le sommet d'avant sur la sequence de type 0 ou 1 n'etait pas deja vu
                        if deja_vu == False :
    #                         if compteur == 11 :
    #                             print("ramou")
    #                             print(voisin_chaine_prec)
                            if position_prec - position_prec_groupe <= 0 :
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                else : 
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                
                            else : 
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
                                else :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec)
                            
                            compteur_tige = compteur
                            if voisin_none != False:

                                positions_ajoutees.append(voisin_none)
                                compteur = compteur+2
                            else :
                                compteur = compteur+1
                                
                                
                            if a_voisin_voisin :
                                voisin_voisin = False
                        else :
                            num_voisin_seq = -1
                           #print(position_prec)
                            for noeud in G.nodes() :
                                for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                                    if  k <= G.nodes[noeud]["position"][1] and k >= G.nodes[noeud]["position"][0] :
                                        num_voisin_seq = noeud
                                        if valeur_debut not in G.nodes[noeud]["chaine"] :
                                            G.nodes[noeud]["chaine"].append(valeur_debut)
                                            
                                        for voisin in G[noeud] :
                                            if G.nodes[voisin]["type"] == None :
                                                G.nodes[voisin]["chaine"].append(valeur_debut)
                                        
                            if num_voisin_seq == -1 :
                                print("probleme3")
                                #print(new_position)
                                #print(position_prec)
                                #print(G.nodes.data())
                            #G.add_edge(num_voisin_seq, compteur_tige, label="B53")
                            
                            if position_prec - position_prec_groupe <= 0 :
                                G.nodes[num_voisin_seq]["position"] = (min(position_prec, G.nodes[num_voisin_seq]["position"][0]), max(position_prec_groupe, G.nodes[num_voisin_seq]["position"][1]))
                            else :
                                G.nodes[num_voisin_seq]["position"] = (min(position_prec_groupe, G.nodes[num_voisin_seq]["position"][0]), max(position_prec, G.nodes[num_voisin_seq]["position"][1]))
                            G.nodes[num_voisin_seq]["poids"] = G.nodes[num_voisin_seq]["position"][1] - G.nodes[num_voisin_seq]["position"][0] +1
    #                 
                            if G.nodes[num_voisin_seq]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                                if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_voisin_seq]["position"][1]) == 1 :
                                    if (num_voisin_seq, compteur_tige) not in G.edges() :
                                        G.add_edge(num_voisin_seq, compteur_tige, label="B53", long_range=False)
                            else :
                                if abs(G.nodes[num_voisin_seq]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                    if (compteur_tige, num_voisin_seq) not in G.edges() :
                                        G.add_edge(compteur_tige, num_voisin_seq, label="B53", long_range=False)
                            
                            compteur_tige = num_voisin_seq
    #                         if type_sommet_prec == 1 :
    #                             if len(chaine) > 0 :
    #                                 for elt in chaine[1] :
    #                                     if elt not in G.nodes[num_voisin_seq]["pos_liaisons"] :
    #                                         G.nodes[num_voisin_seq]["pos_liaisons"].append(elt)
    #                                 for elt in chaine[0] :
    #                                     if valeur_debut not in G.nodes[elt]["pos_liaisons"] and G.nodes[elt]["type"] == 1 :
    #                                         G.nodes[elt]["pos_liaisons"].append(valeur_debut)
    #                                         G.add_edge(num_voisin_seq, elt, label="CWW", long_range=False)
    #                                         G.add_edge(elt, num_voisin_seq, label="CWW", long_range=False)
    #                             else : 
    #                                 G.nodes[num_voisin_seq]["pos_liaisons"].append(chaine)
                        
                    #if nom_cle == ('1FJG', 'A') :
                    #    print(new_position)
                    #    print(num_sommet)
                    #    print(compteur_tige)
                    
                    if i > 0 :
                        if G.nodes[num_sommet]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                            if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_sommet]["position"][1]) == 1 :
                                if (num_sommet, compteur_tige) not in G.edges() :
                                    G.add_edge(num_sommet, compteur_tige, label="B53", long_range=False)
                        else :
                            if abs(G.nodes[num_sommet]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                if (compteur_tige, num_sommet) not in G.edges() :
                                    G.add_edge(compteur_tige, num_sommet, label="B53", long_range=False)
                        G.nodes[num_sommet]["poids"] = 1
                        compteur_tige = num_sommet
                    
                    #if nom_cle == ('1FJG', 'A') :
                    #    print(new_position)
                    #    print(compteur_tige)
                    #    print(G.nodes.data())
                    
                    ## ajout de chacun des voisins du sommet actuel car il est de type 2 ou 3
                    for voisin in voisins :
                        label_voisin = graphes[nom_cle][new_position][voisin]["label"]
                        if label_voisin != "B53" :
                            if voisin not in positions_ajoutees :
                                if label_voisin == 'CWW' and not ((graphes[nom_cle].nodes[new_position]["nt"] == 'A' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'A') or (graphes[nom_cle].nodes[new_position]["nt"] == 'C' and graphes[nom_cle].nodes[voisin]["nt"] == 'G') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'C') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'G'))  : 
                                    label_voisin = 'CWWn'
                                ajout_sommet(G, compteur, compteur_tige, (voisin,voisin), None, 1, label_voisin, positions_ajoutees, int_tige, valeur_debut, long_range=graphes[nom_cle][new_position][voisin]["long_range"])
                                compteur = compteur+1
                            else :
                                num_sommet = -1
                                #print(nom_cle)
                                #print(compteur_tige)
                                #print(G.nodes.data())
                                for noeud in G.nodes() :
    #                                 print(G.nodes[noeud])
                                    if voisin <= G.nodes[noeud]["position"][1] and voisin >= G.nodes[noeud]["position"][0] :
                                        num_sommet = noeud ## retrouver le numero du sommet  
                                if num_sommet == -1 :
                                    print("probleme4")
                                deja_fait = False
                                for v in G[num_sommet] :
                                    for edge in G[num_sommet][v] :
                                        if v == compteur_tige and G[num_sommet][v][edge]["label"] != "B53" : 
                                            deja_fait = True 
                                if deja_fait == False :
                                    if label_voisin == 'CWW' and not ((graphes[nom_cle].nodes[new_position]["nt"] == 'A' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'A') or (graphes[nom_cle].nodes[new_position]["nt"] == 'C' and graphes[nom_cle].nodes[voisin]["nt"] == 'G') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'C') or (graphes[nom_cle].nodes[new_position]["nt"] == 'G' and graphes[nom_cle].nodes[voisin]["nt"] == 'U') or (graphes[nom_cle].nodes[new_position]["nt"] == 'U' and graphes[nom_cle].nodes[voisin]["nt"] == 'G'))  : 
                                        label_voisin = 'CWWn'
                                    
                                    if len(label_voisin) == 3 :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1]
                                    else :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1] + label_voisin[3]
                                    G.add_edge(num_sommet, compteur_tige, label=label_voisin_inverse, long_range=graphes[nom_cle][new_position][voisin]["long_range"])
                                    G.add_edge(compteur_tige, num_sommet, label=label_voisin, long_range=graphes[nom_cle][new_position][voisin]["long_range"])
                    type_sommet_prec = type_sommet_actuel
                    poids_sommet = 1 
            
            
            if voisin_dans_chaine == False :
                chaine_sommet_voisin = -1
                chaine_sommet_voisin_prec = -1
                type_sommet_voisin_prec = -1
                type_sommet_voisin = -1 
            if type_sommet_actuel != 1 :
                voisin_chaine_prec = voisin_chaine
                voisin_chaine = -1             
                     
            i = i+1
            
            if new_position == 1506 :
                print(G.nodes.data())
                print(G.edges.data())
                print("\n")
        
        
            position_prec = new_position
            if valeur_debut == 2 or valeur_debut == 1 :
                new_position = G.node[valeur_debut]["position"][0]-i*int_tige
            else : #normalement 4 ou 3
                new_position = G.node[valeur_debut]["position"][0]+i*int_tige
        else :
            break
            
#     print(chaine)

    print("petit rat")
    print(new_position)
    print(G.nodes.data())
    print(G.edges.data())
    print(voisin_voisin_prec)
    print(voisin_voisin)
    if poids_sommet != 1 or (type_sommet_actuel == 0 or type_sommet_actuel == 1 and poids_sommet == 1):

        deja_vu = False
        for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
            if k in positions_ajoutees :
                deja_vu = True
        #if position_prec not in positions_ajoutees : 
        if deja_vu == False :
            if position_prec - position_prec_groupe <= 0 :
                voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec,position_prec_groupe), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)

            else : 
                voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe,position_prec), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine)
            print("petit petit rat")
            print(G.edges.data())
            
            compteur_tige = compteur
            if voisin_none != False:
                positions_ajoutees.append(voisin_none)
                compteur = compteur+2
            else :
                compteur = compteur+1
            
            if a_voisin_voisin : 
                voisin_voisin = False
        else : 
            num_sommet = -1
            for noeud in G.nodes() :
                for k in range(min(position_prec_groupe, position_prec), max(position_prec, position_prec_groupe)+1) :
                    if  k <= G.nodes[noeud]["position"][1] and k >= G.nodes[noeud]["position"][0] :
                        num_sommet = noeud
                        if valeur_debut not in G.nodes[noeud]["chaine"] :
                            G.nodes[noeud]["chaine"].append(valeur_debut)
                            
                        for voisin in G[noeud] :
                            if G.nodes[voisin]["type"] == None :
                                G.nodes[voisin]["chaine"].append(valeur_debut)
            if G.nodes[num_sommet]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_sommet]["position"][1]) == 1 :
                    if (num_sommet, compteur_tige) not in G.edges() :
                        G.add_edge(num_sommet, compteur_tige, label="B53", long_range=False)
            else :
                if abs(G.nodes[num_sommet]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                    if (compteur_tige, num_sommet) not in G.edges() :
                        G.add_edge(compteur_tige, num_sommet, label="B53", long_range=False)
                
            if position_prec - position_prec_groupe <= 0 :
#                 print(num_sommet)
                G.nodes[num_sommet]["position"] = (min(position_prec, G.nodes[num_sommet]["position"][0]), max(position_prec_groupe, G.nodes[num_sommet]["position"][1]))
            else :
                G.nodes[num_sommet]["position"] = (min(position_prec_groupe, G.nodes[num_sommet]["position"][0]), max(position_prec, G.nodes[num_sommet]["position"][1]))
            G.nodes[num_sommet]["poids"] = G.nodes[num_sommet]["position"][1] - G.nodes[num_sommet]["position"][0] +1   
            
#             if type_sommet_actuel == 1 :
# #                 print(chaine)
# #                 print(num_sommet)
#                 if len(chaine) > 0 :
#                     for elt in chaine[1] :
#                         if elt not in G.nodes[num_sommet]["pos_liaisons"] :
#                             G.nodes[num_sommet]["pos_liaisons"].append(elt)
#                     for elt in chaine[0] :
#                         if valeur_debut not in G.nodes[elt]["pos_liaisons"] and G.nodes[elt]["type"] == 1 :
#                             G.nodes[elt]["pos_liaisons"].append(valeur_debut)
#                             G.add_edge(num_sommet, elt, label="CWW", long_range=False)
#                             G.add_edge(elt, num_sommet, label="CWW", long_range=False)
#                 else : 
#                     G.nodes[num_sommet]["pos_liaisons"].append(chaine)               
#     
             
    return G

''' ajout au graphe d'extensions d'un sommet et d'une arete (qui les relie) pour chaque sommet de type 0 (pour pouvoir compter uniquement le nombre d'aretes pendant la comparaison)  '''
def ajout_aretes_artificielles(G):
    a_ajoute = []
    compteur = 1
    for noeud, label in G.nodes(data="type") :
        if label == 0 or label == 1 :
            a_voisin = False
            for voisin in G[noeud] :
                for edge in G[noeud][voisin] :
                    if G[noeud][voisin][edge]["label"] != 'B53' :
                        a_voisin = True
            if a_voisin == False :
                a_ajoute.append(noeud)
        compteur = noeud
    for elt in a_ajoute :
        compteur += 1
        print(compteur)
        G.add_node(compteur, position=(-1,-1), type=-1, poids=G.nodes[elt]["poids"], chaine=G.nodes[elt]["chaine"])
        G.add_edge(compteur, elt, label='0', long_range=None)
        G.add_edge(elt, compteur, label='0', long_range=None)
    
    a_enlever = []
    for noeud, data in G.nodes(data=True) :
        if data["type"] == -1 and len(G[noeud]) == 0 : ## si un noeud artificiel ne sert plus a rien on l enleve
            a_enlever.append(noeud)
            
    for elt in a_enlever :
        G.remove_node(elt)
    
    return G

''' en cas de probleme de liaisons B53 qui manquent, on les rajoute au graphe d'extension '''
def ajout_liaisons_B53_qui_manquent(G):
    for noeud, pos in G.nodes(data="position") :
        if G.nodes[noeud]["type"] != -1 and G.nodes[noeud]["type"] != None and noeud != 5:
            a_voisin_B53 = False
            for voisin in G[noeud] :
                for edge in G[noeud][voisin] :
                        if G[noeud][voisin][edge]["label"] == 'B53' :
                            if pos[1] + 1 == G.nodes[voisin]["position"][0] :
                                a_voisin_B53 = True
            le_bon_voisin = -1
            if a_voisin_B53 == False :
                for noeud_autre, pos_autre in G.nodes(data="position") :
                    if G.nodes[noeud_autre]["type"] != -1 and G.nodes[noeud_autre]["type"] != None and noeud_autre != noeud and noeud_autre != 5 :
                        if pos_autre[0] - 1 == pos[1] :
                            le_bon_voisin = noeud_autre
                            
            if le_bon_voisin != -1 : 
                G.add_edge(noeud, le_bon_voisin, label='B53', long_range=False)
    return G

''' renvoie les voisins du sommet passe en parametre, lies a lui par des liaisons non-cov (numero du voisin + numero de l'arete)
dans le graphe d'extension'''
def non_cov(numero_sommet, graphe):
    non_cov = []
    for voisin in graphe[numero_sommet] :
        for edge in graphe[numero_sommet][voisin] :
            if graphe[numero_sommet][voisin][edge]["label"] != 'B53' :
                non_cov.append((voisin, edge))
    return non_cov

''' renvoie les voisins du sommet passe en parametre, lies a lui par des liaisons non-cov (numero du voisin + numero de l'arete)
dans le graphe de structure de depart'''
def non_cov_graphe_base(numero_sommet, graphe_mol):
    non_cov = []
    for voisin in graphe_mol[numero_sommet] :
        if graphe_mol.edges[numero_sommet,voisin]["label"] != 'B53' :
            non_cov.append(voisin)
    return non_cov

''' renvoie le voisin (s'il existe, renvoie -1 sinon) du sommet passe en parametre, lie a lui par une liaison covalente dans le sens de la sequence
dans le graphe d'extension '''
def voisin_cov_sens(numero_sommet, graphe):
    for voisin in graphe[numero_sommet] :
        for edge in graphe[numero_sommet][voisin] :
            if graphe[numero_sommet][voisin][edge]["label"] == 'B53' :
                return (voisin, edge)
    return -1

''' renvoie le voisin (s'il existe, renvoie -1 sinon) du sommet passe en parametre, lie a lui par une liaison covalente dans le sens inverse de la sequence
dans le graphe d'extension '''
def voisin_cov_nonsens(numero_sommet, graphe):
    for voisin in graphe.predecessors(numero_sommet) :
        for edge in graphe[voisin][numero_sommet] :
            if graphe[voisin][numero_sommet][edge]["label"] == 'B53' :
                return (voisin, edge)
    return -1

''' regroupe dans le graphe d'extension les liaisons non canoniques short range se trouvant a l'interieur d'helices de liaisons canoniques,
dans le sommet correspondant a l'helice de liaisons canoniques  '''
def regroupement_liaisons_short_range(G, graphe_mol): ## pour des helices avec des non canoniques (parfaites ou non)
    print("ramou")
    print(G.nodes.data())
    print(G.edges.data())
    nx.set_node_attributes(G, 0, "nb_non_can")
    a_voir = []
    a_enlever_noeud = []
    for noeud, data in G.nodes(data=True) :
        if noeud not in a_enlever_noeud :
            ##Eligibilite du sommet pour un regroupement##
            if data["type"] == 3 : ## on sait que ce sommet est implique dans au moins une liaison non canonique ou une liaison canonique long range
                voisin_non_cov = non_cov(noeud, G)
                if len(voisin_non_cov) == 1 and G.nodes[voisin_non_cov[0][0]]["type"] in [None, -1] and G[noeud][voisin_non_cov[0][0]][voisin_non_cov[0][1]]["long_range"] == False : ## on a mnt une seule liaison qui nest pas long range donc cest une liaison non canonique
                    if len(non_cov(voisin_non_cov[0][0], G)) == 1 : ## le voisin du sommet na pas dautre liaison non covalente
                        print("candidat : %s"%noeud)
                        ##Eligibilite du voisinage pour un regroupement##
                        groupe1 = []
                        nb_non_can = 1
                        
                        voisin_cov = voisin_cov_sens(noeud, G)
                        vieux_voisin = (-1,-1)
                        vieux_voisin_non_cov = non_cov_graphe_base(G.nodes[noeud]["position"][0], graphe_mol)
                        while voisin_cov != vieux_voisin and voisin_cov != -1 :
                            if (G.nodes[voisin_cov[0]]["type"] == 1 and vieux_voisin == (-1,-1)) or (G.nodes[voisin_cov[0]]["type"] in [1,3] and vieux_voisin != (-1,-1))  : ##le premier nt rencontre doit etre implique dans une liaison canonique, pour les autres, ce peut etre une liaison can ou non can
                                    vieux_voisin = voisin_cov
                                    voisin_non_cov = non_cov(voisin_cov[0], G)
                                    
                                    voisin_voisin_non_cov = non_cov_graphe_base(G.nodes[voisin_cov[0]]["position"][0], graphe_mol)
                                    print(voisin_voisin_non_cov)
                                    if len(voisin_non_cov) == 1 and G.nodes[voisin_non_cov[0][0]]["type"] in [None, -1] and G[voisin_cov[0]][voisin_non_cov[0][0]][voisin_non_cov[0][1]]["long_range"] != True \
                                    and len(non_cov(voisin_non_cov[0][0], G)) == 1 \
                                    and (vieux_voisin_non_cov == -1 or abs(voisin_voisin_non_cov[0] - vieux_voisin_non_cov[0]) < 10) :
                                        #print("ramousnif")
                                        groupe1.append(voisin_cov) 
                                        print(voisin_cov)
                                        print(G.nodes[voisin_cov[0]]["type"])
                                        print(G.nodes[voisin_cov[0]]["type"] == 3)
                                        if G.nodes[voisin_cov[0]]["type"]  == 3 :
                                            nb_non_can += 1
                                        voisin_cov = voisin_cov_sens(voisin_cov[0], G)
                                        vieux_voisin_non_cov = voisin_voisin_non_cov
                            else :
                                vieux_voisin = voisin_cov
                                   
                            print(voisin_cov)
                            print(vieux_voisin)
                        print("rmou")
                        print(groupe1)
                        if len(groupe1) > 0 and G.nodes[groupe1[len(groupe1)-1][0]]["type"] == 3  :
                            print("rmou")
                            ## on ne va enlever le dernier que si le nt suivant (qui n est pas dans l extension)  n est pas dans  l helice
                            voisin = non_cov_graphe_base(G.nodes[groupe1[len(groupe1)-1][0]]["position"][0], graphe_mol)
                            pos_suivant = G.nodes[groupe1[len(groupe1)-1][0]]["position"][1] + 1
                            can_short_range_new = nb_canoniques_short_range(pos_suivant, graphe_mol)
                            non_can_short_range_new, type_short_range_non_can = nb_non_canoniques_short_range(pos_suivant, graphe_mol)
                            sommet_non_cov = nb_non_cov(pos_suivant, graphe_mol)
                            
                            ok = False
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_short_range_new in d["position"] and d["type"] != None :
                                    dedans = True
                            
                            if can_short_range_new != -1 and len(non_can_short_range_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques_short_range(can_short_range_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques_short_range(can_short_range_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_short_range_new, graphe_mol) == 1)) and (abs(voisin[0] - can_short_range_new) < 10 or abs(voisin[0] - can_short_range_new) < 10) :
                                    ok = True
                            if not ok :
                                groupe1.remove(groupe1[len(groupe1)-1])
                                nb_non_can -= 1   
                                
                        if len(groupe1) == 0 and voisin_cov_sens(noeud, G) == -1 :
                            ## si la liaison non can est en bout de chaine, il faut quand meme verifier le nt d apres
                            voisin = non_cov_graphe_base(G.nodes[noeud]["position"][0], graphe_mol)
                            pos_suivant = G.nodes[noeud]["position"][1] +1
                            can_short_range_new = nb_canoniques_short_range(pos_suivant, graphe_mol)
                            non_can_short_range_new, type_short_range_non_can = nb_non_canoniques_short_range(pos_suivant, graphe_mol)
                            sommet_non_cov = nb_non_cov(pos_suivant, graphe_mol)
                            
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_short_range_new in d["position"] and d["type"] != None :
                                    dedans = True
                                    
                            if can_short_range_new != -1 and len(non_can_short_range_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques_short_range(can_short_range_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques_short_range(can_short_range_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_short_range_new, graphe_mol) == 1)) and (abs(voisin[0] - can_short_range_new) < 10 or abs(voisin[0] - can_short_range_new) < 10) :
   
                                        groupe1.append(0)
                        
                        
                        groupe2 = []
                        voisin_cov = voisin_cov_nonsens(noeud, G)
                        vieux_voisin = (-1,-1)
                        vieux_voisin_non_cov = non_cov_graphe_base(G.nodes[noeud]["position"][0], graphe_mol)
                        while voisin_cov != vieux_voisin and voisin_cov != -1 :
                            
                            if (G.nodes[voisin_cov[0]]["type"] == 1 and vieux_voisin == (-1,-1)) or (G.nodes[voisin_cov[0]]["type"] in [1,3] and vieux_voisin != (-1,-1))  : ##le premier nt rencontre doit etre implique dans une liaison canonique, pour les autres, ce peut etre une liaison can ou non can
                                vieux_voisin = voisin_cov
                                voisin_non_cov = non_cov(voisin_cov[0], G)
                            
                                
                                voisin_voisin_non_cov = non_cov_graphe_base(G.nodes[voisin_cov[0]]["position"][len(G.nodes[voisin_cov[0]]["position"])-1], graphe_mol)
                                
                                if len(voisin_non_cov) == 1 and  G.nodes[voisin_non_cov[0][0]]["type"] in [None, -1] and G[voisin_cov[0]][voisin_non_cov[0][0]][voisin_non_cov[0][1]]["long_range"] != True \
                                and len(non_cov(voisin_non_cov[0][0], G)) == 1 \
                                and (vieux_voisin_non_cov == -1 or abs(voisin_voisin_non_cov[0] - vieux_voisin_non_cov[0]) < 10) :
                                    groupe2.append(voisin_cov) 
                                    if G.nodes[voisin_cov[0]]["type"] == 3 :
                                        nb_non_can += 1
                                    voisin_cov = voisin_cov_nonsens(voisin_cov[0], G)
                                    vieux_voisin_non_cov = voisin_voisin_non_cov
                            else :
                                vieux_voisin = voisin_cov
                        print("rmou")
                        print(groupe2)
                        if len(groupe2) > 0 and G.nodes[groupe2[len(groupe2)-1][0]]["type"] == 3  :
                            
                            ## on ne va enlever le dernier que si le nt suivant (qui n est pas dans l extension)  n est pas dans  l helice
                            voisin = non_cov_graphe_base(G.nodes[groupe2[len(groupe2)-1][0]]["position"][0], graphe_mol)
                            pos_suivant = G.nodes[groupe2[len(groupe2)-1][0]]["position"][0] -1
                            can_short_range_new = nb_canoniques_short_range(pos_suivant, graphe_mol)
                            non_can_short_range_new, type_short_range_non_can = nb_non_canoniques_short_range(pos_suivant, graphe_mol)
                            sommet_non_cov = nb_non_cov(pos_suivant, graphe_mol)
                            
                            ok = False
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_short_range_new in d["position"] and d["type"] != None :
                                    dedans = True
                            
                            if can_short_range_new != -1 and len(non_can_short_range_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques_short_range(can_short_range_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques_short_range(can_short_range_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_short_range_new, graphe_mol) == 1)) and (abs(voisin[0] - can_short_range_new) < 10 or abs(voisin[0] - can_short_range_new) < 10) :
                                    ok = True
                                
                            if not ok :
                                groupe2.remove(groupe2[len(groupe2)-1])
                                nb_non_can -= 1 
                        
                        if len(groupe2) == 0 and voisin_cov_nonsens(noeud, G) == -1 :
                            print("gros rat")
                            ## si la liaison non can est en bout de chaine, il faut quand meme verifier le nt d apres
                            voisin = non_cov_graphe_base(G.nodes[noeud]["position"][0], graphe_mol)
                            pos_suivant = G.nodes[noeud]["position"][0] -1
                            can_short_range_new = nb_canoniques_short_range(pos_suivant, graphe_mol)
                            non_can_short_range_new, type_short_range_non_can = nb_non_canoniques_short_range(pos_suivant, graphe_mol)
                            sommet_non_cov = nb_non_cov(pos_suivant, graphe_mol)
                            
                            print(G.nodes[noeud]["position"][0])
                            print(voisin)
                            print(can_short_range_new)
                            print(non_can_short_range_new)
                            print(sommet_non_cov)
                            #print(nb_non_cov(can_short_range_new, graphe_mol))
                            print()
                            
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_short_range_new in d["position"] and d["type"] != None :
                                    dedans = True
                            
                            if can_short_range_new != -1 and len(non_can_short_range_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques_short_range(can_short_range_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques_short_range(can_short_range_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_short_range_new, graphe_mol) == 1)) and (abs(voisin[0] - can_short_range_new) < 10 or abs(voisin[0] - can_short_range_new) < 10) :

                                    groupe2.append(0)
                                    
                        print(groupe1)
                        print(groupe2)   
                        print(nb_non_can)        
                        if len(groupe1) > 0 and len(groupe2) > 0 :
                            for elt in groupe1 :
                                if elt != 0 :
                                    a_enlever_noeud.append(elt[0])
                            for elt in groupe2 :
                                if elt != 0 :
                                    a_enlever_noeud.append(elt[0])
                            
                            a_voir.append((noeud, list(groupe1), list(groupe2), nb_non_can))
                            
    
    print(a_voir)
    for elt in a_voir : 
        print(elt)
        noeud = elt[0]
        groupe1 = elt[1]
        groupe2 = elt[2]
        nb_non_can = elt[3]
        pos_2 = G.nodes[noeud]["position"][0]
        pos_1 = G.nodes[noeud]["position"][1]
        voisin_non_cov = non_cov(noeud, G)
        print(voisin_non_cov)
        G.remove_node(voisin_non_cov[0][0])
        if len(groupe1) > 0 and groupe1[0] != 0:
            suivant = voisin_cov_sens(groupe1[len(groupe1)-1][0], G)
            pos_1 = G.nodes[groupe1[len(groupe1)-1][0]]["position"][1]
            if suivant != -1 :
                G.add_edge(noeud, suivant[0], label='B53', long_range=False)
        if len(groupe2) > 0 and groupe2[0] != 0:
            prec = voisin_cov_nonsens(groupe2[len(groupe2)-1][0], G)
            pos_2 = G.nodes[groupe2[len(groupe2)-1][0]]["position"][0]
            if prec != -1 :
                G.add_edge(prec[0], noeud, label='B53', long_range=False)
        print(noeud)
        print(G.nodes[noeud])
        poids = G.nodes[noeud]["poids"]
        for g1 in groupe1 :
            if g1 != 0 :
                poids += G.nodes[g1[0]]["poids"]
                if G.nodes[g1[0]]["type"] == 3 :
                    G.remove_node(non_cov(g1[0], G)[0][0])
                G.remove_node(g1[0])
            
            
        for g2 in groupe2 :
            if g2 != 0 :
                poids += G.nodes[g2[0]]["poids"]
                if G.nodes[g2[0]]["type"] == 3 :
                    G.remove_node(non_cov(g2[0], G)[0][0])
                G.remove_node(g2[0])
            
        G.nodes[noeud]["poids"] = poids
        G.nodes[noeud]["nb_non_can"] = nb_non_can
        G.nodes[noeud]["type"] = 1
        G.nodes[noeud]["position"] = (pos_2, pos_1)
        print(G.nodes.data())
                                
    return G                  
    
    
def obtenir_extension(taille_ext) :                
    with open("fichiers_pickle/a-minor_test2.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        tab_aminor = mon_depickler.load()
        
        with open("graphs_2.92.pickle", 'rb') as fichier_tout :
            mon_depickler_graphes = pickle.Unpickler(fichier_tout)
            graphes = mon_depickler_graphes.load()
            
            with open("fichier_test_graphes_version3_avec_boucles_test_digraph_new_"+str(taille_ext)+"cww_non_can_new.txt", "w") as fichier :
                 
                compteur_nb = 0
                compteur_nb_2 = 0
                for occ in tab_aminor :
                    #print(occ)
                    
                    le_bon = False
                    pas_bon = False
                    for elt in liste :
                        if occ["num_PDB"]+"_"+ occ["num_ch"]+"_"+ str(occ["num_motif"])+"_"+ str(occ["num_occ"]) == elt :
                            pas_bon = True
    #                 if occ["num_PDB"] == '5J7L' and occ["num_ch"] == 'DA' and occ["num_motif"] == 191 and occ["num_occ"] == 4 :
    #                         print("ramounsnif 3")
    #                         print(pas_bon)       
                    if occ["num_PDB"] == '5FDU' and occ["num_ch"] == '1A' and occ["num_motif"] == 48 and occ["num_occ"] == 25 :
                        le_bon = True
    
                    if pas_bon == False :## and le_bon :
                        
                        print(str(occ["num_PDB"]) + "_" + str(occ["num_ch"]) + "_" + str(occ["num_motif"]) + "_" + str(occ["num_occ"]))
                        
                        positions_chaines = [[],[],[],[]]
                         
                        G = nx.MultiDiGraph()
                        i = 1
                        positions_ajoutees = []
                        for elt in occ["a_minor"] :
                            G.add_node(i, position = (elt,elt), type = i+10, poids=1, chaine=[i])
                            positions_ajoutees.append(elt)
                            
                            i = i+1
                        G.add_edge(1,2, label="CSS", long_range = True)
                        G.add_edge(2,1, label="CSS", long_range = True)
                        if(G.node[1]["position"][0] < G.node[3]["position"][0]) :
                            G.add_edge(1,3, label="B53", long_range = False)
                        else :
                            G.add_edge(3,1, label="B53", long_range = False)
                        G.add_edge(1,5, label="TSS", long_range = True)
                        G.add_edge(5,1, label="TSS", long_range = True)
                        if(G.node[2]["position"][0] < G.node[4]["position"][0]) :
                            G.add_edge(2,4, label="B53", long_range = False)
                        else :
                            G.add_edge(4,2, label="B53", long_range = False)
                        G.add_edge(2,5, label="CWW", long_range = False)
                        G.add_edge(5,2, label="CWW", long_range = False)
                        G.add_edge(3,4, label="CSS", long_range = True)
                        G.add_edge(4,3, label="CSS", long_range = True)
                        compteur = 6
                        compteur_tige = 2
                         
    #                     print(G.nodes.data())
                        #print(G.edges.data())
                        #fichier.write(str(graphes[('1U9S', 'A')].nodes.data())+"\n")
                        #fichier.write(str(graphes[('1U9S', 'A')].edges.data()))
                        nom_cle = (occ["num_PDB"], occ["num_ch"])
                        #print(str(nom_cle) + str(occ["num_motif"]) + str(occ["num_occ"])) 
                         
    #                     if nom_cle == ('4V88', 'A6') and occ["num_motif"] == 17 and occ["num_occ"] == 55 :
    #                         print(nom_cle)
    #                         print(occ["num_motif"])
    #                         print(occ["num_occ"])
                         
                        ##Garder toutes les positions de la sequence 
                        if G.node[2]["position"][0] - G.node[4]["position"][0] < 0 : ## ordre de haut en bas
                            for i in range(taille_ext) :
                                if occ["a_minor"][1]-i > 0 :
                                    positions_chaines[1].append(occ["a_minor"][1]-i)
                            for i in range(taille_ext) :
                                if occ["a_minor"][3]+i < graphes[nom_cle].number_of_nodes() :
                                    positions_chaines[3].append(occ["a_minor"][3]+i)
                        else : ## ordre de bas en haut
                            for i in range(taille_ext) :
                                if occ["a_minor"][1]+i < graphes[nom_cle].number_of_nodes() :
                                    positions_chaines[1].append(occ["a_minor"][1]+i)
                            for i in range(taille_ext) :
                                if occ["a_minor"][3]-i > 0 :
                                    positions_chaines[3].append(occ["a_minor"][3]-i)
                                    
                                    
                        if G.node[1]["position"][0] - G.node[3]["position"][0] < 0 : ## ordre de haut en bas
                            for i in range(taille_ext) :
                                if occ["a_minor"][0]-i > 0 :
                                    positions_chaines[0].append(occ["a_minor"][0]-i)
                            for i in range(taille_ext) :
                                if occ["a_minor"][2]+i < graphes[nom_cle].number_of_nodes() :
                                    positions_chaines[2].append(occ["a_minor"][2]+i)
                        else : ## ordre de bas en haut
                            for i in range(taille_ext) :
                                if occ["a_minor"][0]+i < graphes[nom_cle].number_of_nodes() :
                                    positions_chaines[0].append(occ["a_minor"][0]+i)
                            for i in range(taille_ext) :
                                if occ["a_minor"][2]-i > 0 :
                                    positions_chaines[2].append(occ["a_minor"][2]-i)
                         
                          
                        ##Connaitre l'ordre des positions sur la sequence
                        int_tige = 0
                        if G.node[2]["position"][0] - G.node[4]["position"][0] < 0 : ## ordre de haut en bas
                            int_tige = 1
                        else : ## ordre de bas en haut
                            int_tige = -1
                        #compteur_tige2 = G.node[5]["position"][0] + int_tige
                         
                        ##Initialisation de la boucle 1
                        #===================================================================
                        G = extension_tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines)
                        #print(G.nodes.data())
                        #print(positions_ajoutees)
                        #if nom_cle == ('1FJG', 'A') :
                        #    print(G.nodes.data())
                        #    print(G.edges.data())
                         
                        compteur = G.number_of_nodes()+1
                             
                        #compteur_tige2 = G.node[5]["position"][0] - int_tige
                        compteur_tige = 4
                         
    #                     print(nom_cle)
    #                     print(occ["num_motif"])
    #                     print(occ["num_occ"])
                         
                        G = extension_tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines)  
     
                        compteur_tige = 1 
                        compteur = G.number_of_nodes()+1
                        if G.node[1]["position"][0] - G.node[3]["position"][0] < 0 : ## ordre de haut en bas
                            int_tige = 1
                        else : ## ordre de bas en haut
                            int_tige = -1
                          
                         
                          
                        G = extension_tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines) 
                         
     
                        compteur_tige = 3
                        compteur = G.number_of_nodes()+1
                         
                        G = extension_tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines) 
                         
                        
                        
                        G = ajout_liaisons_B53_qui_manquent(G)
                         
                        print(G.nodes.data())
                        print(G.edges.data()) 
                        G = ajout_aretes_artificielles(G)
                        G = regroupement_liaisons_short_range(G, graphes[nom_cle])
                       
                        G = ajout_aretes_artificielles(G)
                        print(G.nodes.data())
                        print(G.edges.data())
                        print(G.number_of_edges())
                       
                        #print(occ)
    #                     print(G.nodes.data())
    #                     print(G.edges.data())
                        
    #                     compteur = G.number_of_nodes()+1
    #                     compteur_tige = 1
                        #G = tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees) 
#                         for noeud, data in G.nodes(data=True) :
#                             if noeud not in G2.nodes or data["type"] != G2.nodes[noeud]["type"]:
#                                 print(noeud)
                        #print(G.nodes.data())
                        #print(G.edges.data()) 
    
                        fichier.write(str(nom_cle) + " " + str(occ["num_motif"]) + " " + str(occ["num_occ"]) + "\n")
                        for noeud, data in G.nodes(data=True) :
                            fichier.write(str(noeud) + " " + str(data) +"\n")
                                    
                        for u,v,data in G.edges(data=True) :
                            fichier.write(str((u,v)) + " " + str(data)+"\n")
# # #                         
# # #     #        
                        with open(EXTENSION_PATH_TAILLE%taille_ext+"fichier_{}_3.pickle".format(str(occ["num_PDB"]) + "_" + str(occ["num_ch"]) + "_" + str(occ["num_motif"]) + "_" + str(occ["num_occ"])), "wb") as fichier_sortie :
                            mon_pickler = pickle.Pickler(fichier_sortie)
                            mon_pickler.dump(G)
                        
#                         for noeud, data in G.nodes(data=True) :
#                             print(str(noeud) + " " + str(data) +"\n")
#                                     
#                         for u,v,data in G.edges(data=True) :
#                             print(str((u,v)) + " " + str(data)+"\n")    
                    
            
            #pos=nx.spring_layout(G)
        
            #===================================================================

if __name__ == '__main__':
    for i in range(1,11) :
        obtenir_extension(i)
#     obtenir_extension(10)