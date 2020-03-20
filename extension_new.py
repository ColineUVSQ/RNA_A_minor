'''
Created on 2 sept. 2019

@author: coline

Adaptation de la construction des extensions aux nouvelles donnees (toutes PDB)

'''
import networkx as nx
import pickle
from recup_data.constantes import NEW_EXTENSION_PATH_TAILLE
from recup_data.environnement_liaisons_non_can import nb_canoniques_short_range,\
    nb_non_canoniques_short_range
    
''' ajout d'un sommet au graphe d'extension et d'aretes si necessaire  '''
def ajout_sommet(G, compteur, compteur_tige, position, typ, poids, label, positions_ajoutees, int_tige, valeur_debut, *args, **kwargs):
    near = kwargs.get("long_range", False)
    num_ch = kwargs.get("num_ch", -1)
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
    
#     print("position")
#     print(position)
    G.add_node(compteur, position=position, type=typ, poids=poids, chaine=[valeur_debut], num_ch = num_ch)
    
#     print("gros rat")
#     print(voisin_chaine)
#     print(voisin_voisin)
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
#             print(voisin_chaine)
#             print(data["position"][0])
            
            if voisin_chaine >= data["position"][0] and voisin_chaine <= data["position"][1] :
                noeud_existe = noeud
        if noeud_existe != False :
            G.add_edge(compteur, noeud_existe, label="CWW", near = near)
            G.add_edge(noeud_existe, compteur, label="CWW", near = near)
        
        
        
        ##cas ou une liaison can avec un sommet None
        if voisin_voisin != False : 
            noeud_existe = False
            for noeud, data in G.nodes(data=True) :
                if voisin_voisin >= data["position"][0] and voisin_voisin <= data["position"][1] :
                    noeud_existe = noeud
            if noeud_existe :
                deja = False
                if (compteur,noeud_existe) in G.edges() :
                    for edge in G[compteur][noeud_existe] :
                        if G[compteur][noeud_existe][edge]["label"] == 'CWW' :
                            deja = True
                if (compteur, noeud_existe) not in G.edges() or not deja :        
                    G.add_edge(compteur, noeud_existe, label="CWW", near = near)
                    G.add_edge(noeud_existe, compteur, label="CWW", near = near)
            else :
                G.add_node(compteur+1, position=(voisin_voisin, voisin_voisin), type=None, poids = 1, chaine=[valeur_debut], num_ch = num_ch)
                G.add_edge(compteur, compteur+1, label="CWW", near = near)
                G.add_edge(compteur+1, compteur, label="CWW", near = near)

                voisin_none = (num_ch, voisin_voisin)
            a_voisin_voisin = True
        
    if label == "B53" :
        if G.nodes[compteur]["position"][0] < G.nodes[compteur_tige]["position"][0] :
            if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[compteur]["position"][1]) == 1 :
                if (compteur, compteur_tige) not in G.edges() :
                    G.add_edge(compteur, compteur_tige, label="B53", near=False)
        else :
            if abs(G.nodes[compteur]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                if (compteur_tige, compteur) not in G.edges() :
                    G.add_edge(compteur_tige, compteur, label="B53", near=False)
    else :
        if len(label) == 3 :
            label_inverse = label[0]+label[2]+label[1]
        else :
            label_inverse = label[0]+label[2]+label[1] + label[3]
        G.add_edge(compteur_tige, compteur, label=label, near=near)
        G.add_edge(compteur, compteur_tige, label=label_inverse, near=near)
    for i in range(position[0], position[1]+1) :
        positions_ajoutees.append((num_ch,i))
    
#     print("ramousnif")
#     print(voisin_voisin)
#     print(voisin_none)
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

''' renvoie la liste des types near des liaisons CWW canoniques stockees dans un sommet de type 1 (pos) '''
def retrouver_attr_near_type_1(graphe, pos):
    type_near = -1
    for voisin in graphe.successors(pos) :
        if graphe.edges[pos,voisin]["label"] == 'CWW' :
            type_near = graphe.edges[pos,voisin]["near"]
                
    return type_near
            

''' renvoie le type du sommet new_position
voisins : l'ensemble de ses voisins
G : le graphe d'extension qui est en train d'etre cree
graphe_base : le graphe de structure dont on part pour creer le graphe d'extension
taille_motif : le nombre de nucleotides contenant le motif qu'on veut etendre (5 pour A-minor) '''
def type_sommet(voisins, new_position, G, graphe_base, taille_motif):
    
    for i in range(1,taille_motif) :
        if new_position == G.nodes[i]["position"][0] :
            return i+10
    
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
                #if label_voisin == 'CWW' and ((graphe_base.nodes[new_position]["nt"] == 'A' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'A') or (graphe_base.nodes[new_position]["nt"] == 'C' and graphe_base.nodes[voisin]["nt"] == 'G') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'C') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                if label_voisin == 'CWW' :
                ##sommet de type 0
                    type_sommet_actuel = 1       
                    #compteur_tige2 = voisin
                else : ##sommet de type 3

                    type_sommet_actuel = 3
    else : ##Plus d'un autre voisin en dehors de la sequence
        for voisin in voisins : #Recherche d'une liaison can de tige
            label_voisin = graphe_base[new_position][voisin]["label"]
            if label_voisin != "B53" :
                #if label_voisin == 'CWW' and ((graphe_base.nodes[new_position]["nt"] == 'A' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'A') or (graphe_base.nodes[new_position]["nt"] == 'C' and graphe_base.nodes[voisin]["nt"] == 'G') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'C') or (graphe_base.nodes[new_position]["nt"] == 'G' and graphe_base.nodes[voisin]["nt"] == 'U') or (graphe_base.nodes[new_position]["nt"] == 'U' and graphe_base.nodes[voisin]["nt"] == 'G'))  : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                if label_voisin == 'CWW' :
                    type_sommet_actuel = 2
                    #compteur_tige2 = voisin
        if type_sommet_actuel == -1 :
            type_sommet_actuel = 3
            
    return type_sommet_actuel


''' 02/09/19
idem que extension tige
mais adaptee a la forme des nouvelles donnees
(sans prendre en compte les long_range et en considerant les a minor entree chaines)'''
def extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, taille_motif):
    #print(compteur)
    valeur_debut = compteur_tige
    i = 0
    poids_sommet = 1
    type_sommet_prec = None
    type_sommet_actuel = None 
    type_sommet_voisin = -1
    type_sommet_voisin_prec = -1
   
    if valeur_debut == num_ch_1[0] or valeur_debut == num_ch_2[0] : ## ou juste if int_tige == -1 :
        new_position = (G.node[valeur_debut]["num_ch"], G.node[valeur_debut]["position"][0]-i)
    else :
        new_position = (G.node[valeur_debut]["num_ch"], G.node[valeur_debut]["position"][0]+i) 

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
    num_ch_sommet_prec = -1
    num_ch_sommet = -1
    
    
    while i < taille_ext and new_position[1] > 0 and new_position[1] < graphe.number_of_nodes() and new_position[0] == G.node[valeur_debut]["num_ch"] and new_position in graphe.nodes():
        print(new_position)
        print(G.nodes.data())
        print(G.edges.data())
        if new_position == ('A', 6) :
                                    print("petit rat")
                                    print(type_sommet_actuel)
                                    print(type_sommet_prec)
                                    print(num_ch_sommet)
                                    print(num_ch_sommet_prec)
                                    print(voisin_chaine_prec)
                                    print(G.edges.data())
                                    print(positions_ajoutees)
                                    print(i)
        if i == 0 or new_position not in positions_ajoutees[:taille_motif-1] :## pour ne pas revenir sur le motif
            
            #print(graphe.nodes.data())
            voisins = graphe[new_position]
 
            type_sommet_actuel = type_sommet(voisins, new_position, G, graphe, taille_motif)
            
            num_ch_sommet_prec = num_ch_sommet
            num_ch_sommet = new_position[0]

            if type_sommet_actuel == 1 :
               
                for voisin in voisins :
                    label_voisin = graphe[new_position][voisin]["label"]
                    if label_voisin != "B53" :
                        if label_voisin == 'CWW' : #and ((voisin <= compteur_tige2 and voisin > compteur_tige2-5) or (voisin >= compteur_tige2 and voisin < compteur_tige2+5) or (voisin <= G.node[valeur_debut]["position"][0] and voisin > G.node[valeur_debut]["position"][0]-10) or (voisin >= G.node[valeur_debut]["position"][0] and voisin < G.node[valeur_debut]["position"][0]+10)) :
                            voisins_du_voisin = non_cov_graphe_base(voisin, graphe)
                            
                            voisin_dans_chaine = False
                            voisin_voisin_prec = voisin_voisin
                            #voisin_voisin = False
                            for k in range(4) :
#                                 print(voisin)
#                                 print(positions_chaines)
                                if voisin in positions_chaines[k] : ## pour savoir si le voisin appartient a une des chaines ou pas 
                                    
                                    type_sommet_voisin_prec = type_sommet_voisin
                                    type_sommet_voisin = type_sommet(graphe[voisin],voisin, G, graphe, taille_motif)
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
                                            voisin_voisin = voisin[1]
                            #print(voisin_voisin)
                            
                            voisin_chaine_prec = voisin_chaine
                            voisin_chaine = voisin[1]
                            
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
                            (voisin_voisin != False and \
                             (type_sommet_voisin_prec == type_sommet_voisin and (type_sommet_voisin == 1 or type_sommet_voisin == 0)) and \
                             (chaine_sommet_voisin == chaine_sommet_voisin_prec) or ((voisin_voisin == False) and \
                             (type_sommet_voisin == -1 and type_sommet_voisin_prec == -1)) ) and \
                            (num_ch_sommet == num_ch_sommet_prec or num_ch_sommet_prec == -1 or num_ch_sommet == -1 ) and \
                            ((voisin_dans_chaine == False and (abs(voisin_chaine-voisin_chaine_prec) < 10 or voisin_chaine == -1 or voisin_chaine_prec == -1)) or \
                            (voisin_dans_chaine == True and (abs(voisin_chaine-voisin_chaine_prec) == 1 or voisin_chaine == -1 or voisin_chaine_prec == -1))))): 
                        poids_sommet += 1
                        
                        
                        
    #                 elif type_sommet_actuel == 1 and type_sommet_voisin_prec != type_sommet_voisin and type_sommet_voisin != None  : ## il faut ajouter le sommet precedent
    #                     
    #                     
                    else :
                        
                        if (type_sommet_prec == 1 or type_sommet_prec == 0) :  ## on va ajouter le sommet precedent car il est de type 0 ou 1   
                            
                            type_near = []
                            deja_vu = False
                            for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
                                if (position_prec_groupe[0],k) in positions_ajoutees :
                                    deja_vu = True 
                                near = retrouver_attr_near_type_1(graphe, (G.node[valeur_debut]["num_ch"], k))
                                if near not in type_near :
                                    type_near.append(near)
                                
                                    
                            #if position_prec not in positions_ajoutees :
                            if deja_vu == False :
    #                             if compteur == 38 :
    #                                 print("ramousnif")
    #                                 print(chaine)
                                
                                if position_prec[1] - position_prec_groupe[1] <= 0 :
                                    if type_sommet_actuel != 1 :
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                    else : 
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                    
                                else : 
                                    if type_sommet_actuel != 1 :
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                    else :
                                        voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                                
                                if position_prec[0] != position_prec_groupe[0] : ## differents num de chaine
                                    print("probleme18464")
                                

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
                                        
                                    if position_prec[1] - position_prec_groupe[1] <= 0 :
                                        G.nodes[num_sommet_prec]["position"] = (min(position_prec[1], G.nodes[num_sommet_prec]["position"][0]), max(position_prec_groupe[1], G.nodes[num_sommet_prec]["position"][1]))
                                    else :
                                        G.nodes[num_sommet_prec]["position"] = (min(position_prec_groupe[1], G.nodes[num_sommet_prec]["position"][0]), max(position_prec[1], G.nodes[num_sommet_prec]["position"][1]))
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
                                                G.add_edge(num_sommet_prec, compteur_tige, label="B53", near=False)
                                    else :
                                        if abs(G.nodes[num_sommet_prec]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                            if (compteur_tige, num_sommet_prec) not in G.edges() :
                                                G.add_edge(compteur_tige, num_sommet_prec, label="B53", near=False)
                                    
                                    compteur_tige = num_sommet_prec 
                        position_prec_groupe = new_position    
                        type_sommet_prec = type_sommet_actuel
                        poids_sommet = 1  
                        
                        
                        
                         
                            
                else :
                    
                    if poids_sommet != 1 or (type_sommet_prec == 0 or type_sommet_prec == 1 and poids_sommet == 1): ## on va ajouter le sommet precedent car il est de type 0 ou 1
                        
                        
                        
                        type_near = []
                        deja_vu = False
                        for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
                            if (position_prec_groupe[0],k) in positions_ajoutees :
                                deja_vu = True
                                
                            near = retrouver_attr_near_type_1(graphe, (G.node[valeur_debut]["num_ch"], k))
                            if near not in type_near :
                                type_near.append(near)
                              
                        #if position_prec not in positions_ajoutees : ## le sommet d'avant sur la sequence de type 0 ou 1 n'etait pas deja vu
                        if deja_vu == False :
                            if position_prec[1] - position_prec_groupe[1] <= 0 :
                                if type_sommet_actuel != 1   :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                else : 
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                
                            else : 
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                else :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                            
                            if position_prec[0] != position_prec_groupe[0] : ## differents num de chaine
                                print("probleme48523")
                            
                            compteur_tige = compteur
                            if voisin_none != False:
#                                 print("petit voisin")
#                                 print(voisin_none)
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
                                    for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
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
                                    
                                if position_prec[1] - position_prec_groupe[1] <= 0 :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec[1], G.nodes[num_sommet_prec]["position"][0]), max(position_prec_groupe[1], G.nodes[num_sommet_prec]["position"][1]))
                                else :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec_groupe[1], G.nodes[num_sommet_prec]["position"][0]), max(position_prec[1], G.nodes[num_sommet_prec]["position"][1]))
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
                                            G.add_edge(num_sommet_prec, compteur_tige, label="B53", near=False)
                                else :
                                    if abs(G.nodes[num_sommet_prec]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                        if (compteur_tige, num_sommet_prec) not in G.edges() :
                                            G.add_edge(compteur_tige, num_sommet_prec, label="B53", near=False)
                                
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
                    ajout_sommet(G, compteur, compteur_tige, (new_position[1], new_position[1]), type_sommet_actuel, 1, "B53", positions_ajoutees, int_tige, valeur_debut, num_ch=new_position[0])
                    
                    
                    compteur_tige = compteur
                    compteur = compteur+1
                    
                    ## ajout de chacun des voisins du sommet actuel car il est de type 2 ou 3
                    for voisin in voisins : 
                        
                        
                        label_voisin = graphe[new_position][voisin]["label"]
                        if label_voisin != "B53" :
                            if voisin not in positions_ajoutees :
                                
#                                 if label_voisin == 'CWW' and not ((graphe.nodes[new_position]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[new_position]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : 
#                                     label_voisin = 'CWWn'
                                ajout_sommet(G, compteur, compteur_tige, (voisin[1],voisin[1]), None, 1, label_voisin, positions_ajoutees, int_tige, valeur_debut, near = graphe[new_position][voisin]["near"], num_ch = voisin[0])
                                compteur = compteur+1
                            else :
                                num_sommet = -1
                                #print("ramou")
                                for noeud in G.nodes() :
                                    if  voisin[1] <= G.nodes[noeud]["position"][1] and voisin[1] >= G.nodes[noeud]["position"][0] :
                                        num_sommet = noeud ## retrouver le numero du sommet 
                                if num_sommet == -1 :
                                    print("probleme2")
                                deja_fait = False
                                for v in G[num_sommet] :
                                    for edge in G[num_sommet][v] :
                                        if v == compteur_tige and G[num_sommet][v][edge]["label"] != "B53" : 
                                            deja_fait = True 
                                if deja_fait == False :
#                                     if label_voisin == 'CWW' and not ((graphe.nodes[new_position]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[new_position]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : 
#                                         label_voisin = 'CWWn'
                                    if len(label_voisin) == 3 :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1]
                                    else :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1] + label_voisin[3]
                                    G.add_edge(num_sommet, compteur_tige, label=label_voisin_inverse, near = graphe[new_position][voisin]["near"])
                                    G.add_edge(compteur_tige, num_sommet, label=label_voisin, near = graphe[new_position][voisin]["near"])
                    type_sommet_prec = type_sommet_actuel 
    
            else : ## le sommet actuel existe deja
                
                num_sommet = -1
                #print(G.nodes.data())
                #print(positions_ajoutees)
                #print(new_position)
                for noeud in G.nodes() :
                    if new_position[1] <= G.nodes[noeud]["position"][1] and new_position[1] >= G.nodes[noeud]["position"][0]  :
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
                            ((voisin_voisin != False) and \
                              (type_sommet_voisin_prec == type_sommet_voisin and (type_sommet_voisin == 1 or type_sommet_voisin == 0)) and \
                             (chaine_sommet_voisin == chaine_sommet_voisin_prec) or ((voisin_voisin == False) and \
                             (type_sommet_voisin == -1 and type_sommet_voisin_prec == -1)) ) and \
                            (num_ch_sommet == num_ch_sommet_prec or num_ch_sommet_prec == -1 or num_ch_sommet == -1 ) and \
                            ((voisin_dans_chaine == False and (abs(voisin_chaine-voisin_chaine_prec) < 10 or voisin_chaine == -1 or voisin_chaine_prec == -1)) or \
                            (voisin_dans_chaine == True and (abs(voisin_chaine-voisin_chaine_prec) == 1 or voisin_chaine == -1 or voisin_chaine_prec == -1))))): 
                        poids_sommet += 1
                        
                    else :
                        
                        type_near = []
                        deja_vu = False
                        for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
                            if (position_prec_groupe[0], k) in positions_ajoutees :
                                deja_vu = True
                                
                            near = retrouver_attr_near_type_1(graphe, (G.node[valeur_debut]["num_ch"], k))
                            if near not in type_near :
                                type_near.append(near)
                        #if position_prec not in positions_ajoutees :
                        if new_position == ('CA', 1801) :
                                print("tout petit rat")
                                print(compteur)
                                print(compteur_tige)
                                print(num_sommet)
                                print(position_prec_groupe)
                                print(positions_ajoutees)
                                print(deja_vu)
                        if deja_vu == False :
                            if position_prec[1] - position_prec_groupe[1] <= 0 :
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                else : 
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                
                            else : 
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                else :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                            
                            if position_prec[0] != position_prec_groupe[0] :
                                print("probleme253646")
                            
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
                                    for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
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
                                        
                                if position_prec[1] - position_prec_groupe[1] <= 0 :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec[1], G.nodes[num_sommet_prec]["position"][0]), max(position_prec_groupe[1], G.nodes[num_sommet_prec]["position"][1]))
                                else :
                                    G.nodes[num_sommet_prec]["position"] = (min(position_prec_groupe[1], G.nodes[num_sommet_prec]["position"][0]), max(position_prec[1], G.nodes[num_sommet_prec]["position"][1]))
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
                                                G.add_edge(num_sommet_prec, compteur_tige, label="B53", near=False)
                                    else :
                                        if abs(G.nodes[num_sommet_prec]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                            if (compteur_tige, num_sommet_prec) not in G.edges() :
                                                G.add_edge(compteur_tige, num_sommet_prec, label="B53", near=False)
                                
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
                                    G.add_edge(num_sommet, compteur_tige, label="B53", near=False)
                        else :
                            if abs(G.nodes[num_sommet]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                if (compteur_tige, num_sommet) not in G.edges() :
                                    G.add_edge(compteur_tige, num_sommet, label="B53", near=False)
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
                        type_near = []
                        for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
                            if (position_prec_groupe[0],k) in positions_ajoutees :
                                deja_vu = True
                                
                            near = retrouver_attr_near_type_1(graphe, (G.node[valeur_debut]["num_ch"], k))
                            if near not in type_near :
                                type_near.append(near)
                        
                        #if position_prec not in positions_ajoutees : ## le sommet d'avant sur la sequence de type 0 ou 1 n'etait pas deja vu
                        if deja_vu == False :
    #                         if compteur == 11 :
    #                             print("ramou")
    #                             print(voisin_chaine_prec)
                            if position_prec[1] - position_prec_groupe[1] <= 0 :
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                else : 
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                
                            else : 
                                if type_sommet_actuel != 1 :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
                                else :
                                    voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin_prec, voisin_chaine = voisin_chaine_prec, num_ch = position_prec[0], near = type_near)
                            
                            if position_prec[0] != position_prec_groupe[0] :
                                print("probleme144654")
                            
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
                                for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
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
                            
                            if position_prec[1] - position_prec_groupe[1] <= 0 :
                                G.nodes[num_voisin_seq]["position"] = (min(position_prec[1], G.nodes[num_voisin_seq]["position"][0]), max(position_prec_groupe[1], G.nodes[num_voisin_seq]["position"][1]))
                            else :
                                G.nodes[num_voisin_seq]["position"] = (min(position_prec_groupe[1], G.nodes[num_voisin_seq]["position"][0]), max(position_prec[1], G.nodes[num_voisin_seq]["position"][1]))
                            G.nodes[num_voisin_seq]["poids"] = G.nodes[num_voisin_seq]["position"][1] - G.nodes[num_voisin_seq]["position"][0] +1
    #                 
                            if G.nodes[num_voisin_seq]["position"][0] < G.nodes[compteur_tige]["position"][0] :
                                if abs(G.nodes[compteur_tige]["position"][0] - G.nodes[num_voisin_seq]["position"][1]) == 1 :
                                    if (num_voisin_seq, compteur_tige) not in G.edges() :
                                        G.add_edge(num_voisin_seq, compteur_tige, label="B53", near=False)
                            else :
                                if abs(G.nodes[num_voisin_seq]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                    if (compteur_tige, num_voisin_seq) not in G.edges() :
                                        G.add_edge(compteur_tige, num_voisin_seq, label="B53", near=False)
                            
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
                                    G.add_edge(num_sommet, compteur_tige, label="B53", near=False)
                        else :
                            if abs(G.nodes[num_sommet]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                                if (compteur_tige, num_sommet) not in G.edges() :
                                    G.add_edge(compteur_tige, num_sommet, label="B53", near=False)
                        G.nodes[num_sommet]["poids"] = 1
                        compteur_tige = num_sommet
                    
                    #if nom_cle == ('1FJG', 'A') :
                    #    print(new_position)
                    #    print(compteur_tige)
                    #    print(G.nodes.data())
                    
                    ## ajout de chacun des voisins du sommet actuel car il est de type 2 ou 3
                    for voisin in voisins :
                        label_voisin = graphe[new_position][voisin]["label"]
                        if label_voisin != "B53" :
                            if voisin not in positions_ajoutees :
#                                 if label_voisin == 'CWW' and not ((graphe.nodes[new_position]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[new_position]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : 
#                                     label_voisin = 'CWWn'
                                ajout_sommet(G, compteur, compteur_tige, (voisin[1],voisin[1]), None, 1, label_voisin, positions_ajoutees, int_tige, valeur_debut, near=graphe[new_position][voisin]["near"], num_ch=voisin[0])
                                compteur = compteur+1
                            else :
                                num_sommet = -1
                                #print(nom_cle)
                                #print(compteur_tige)
                                #print(G.nodes.data())
                                for noeud in G.nodes() :
    #                                 print(G.nodes[noeud])
#                                     print(voisin)
#                                     print(noeud)
                                    if voisin[1] <= G.nodes[noeud]["position"][1] and voisin[1] >= G.nodes[noeud]["position"][0] :
                                        num_sommet = noeud ## retrouver le numero du sommet  
                                if num_sommet == -1 :
                                    print("probleme4")
                                deja_fait = False
                                for v in G[num_sommet] :
                                    for edge in G[num_sommet][v] :
                                        if v == compteur_tige and G[num_sommet][v][edge]["label"] != "B53" : 
                                            deja_fait = True 
                                if deja_fait == False :
#                                     if label_voisin == 'CWW' and not ((graphe.nodes[new_position]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[new_position]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[new_position]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[new_position]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))  : 
#                                         label_voisin = 'CWWn'
                                    
                                    if len(label_voisin) == 3 :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1]
                                    else :
                                        label_voisin_inverse = label_voisin[0]+label_voisin[2]+label_voisin[1] + label_voisin[3]
                                    G.add_edge(num_sommet, compteur_tige, label=label_voisin_inverse, near=graphe[new_position][voisin]["near"])
                                    G.add_edge(compteur_tige, num_sommet, label=label_voisin, near=graphe[new_position][voisin]["near"])
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
            

        
        
            position_prec = new_position
            if valeur_debut == num_ch_1[0] or valeur_debut == num_ch_2[0] : # ou juste int_tige == -1
                new_position = (G.node[valeur_debut]["num_ch"], G.node[valeur_debut]["position"][0]-i)
            else : #normalement 4 ou 3
                new_position = (G.node[valeur_debut]["num_ch"], G.node[valeur_debut]["position"][0]+i)
        else :
            break
            
#     print(chaine)

#     print("petit rat")
#     print(new_position)
#     print(G.nodes.data())
#     print(G.edges.data())
#     print(voisin_voisin_prec)
#     print(voisin_voisin)
    if poids_sommet != 1 or (type_sommet_actuel == 0 or type_sommet_actuel == 1 and poids_sommet == 1):
        
        type_near = []
        deja_vu = False
        for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
            if (position_prec_groupe[0], k) in positions_ajoutees :
                deja_vu = True
                
            near = retrouver_attr_near_type_1(graphe, (G.node[valeur_debut]["num_ch"], k))
            if near not in type_near :
                type_near.append(near)
        #if position_prec not in positions_ajoutees : 
        if deja_vu == False :
            if position_prec[1] - position_prec_groupe[1] <= 0 :
                voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec[1],position_prec_groupe[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)

            else : 
                voisin_none, a_voisin_voisin = ajout_sommet(G, compteur, compteur_tige, (position_prec_groupe[1],position_prec[1]), type_sommet_prec, poids_sommet, "B53", positions_ajoutees, int_tige, valeur_debut, voisin_voisin = voisin_voisin, voisin_chaine = voisin_chaine, num_ch = position_prec[0], near = type_near)
            
            if position_prec[0] != position_prec_groupe[0] :
                print("probleme526")
            
#             print("petit petit rat")
#             print(G.edges.data())
            
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
                for k in range(min(position_prec_groupe[1], position_prec[1]), max(position_prec[1], position_prec_groupe[1])+1) :
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
                        G.add_edge(num_sommet, compteur_tige, label="B53", near=False)
            else :
                if abs(G.nodes[num_sommet]["position"][0] - G.nodes[compteur_tige]["position"][1]) == 1 :
                    if (compteur_tige, num_sommet) not in G.edges() :
                        G.add_edge(compteur_tige, num_sommet, label="B53", near=False)
                
            if position_prec[1] - position_prec_groupe[1] <= 0 :
#                 print(num_sommet)
                G.nodes[num_sommet]["position"] = (min(position_prec[1], G.nodes[num_sommet]["position"][0]), max(position_prec_groupe[1], G.nodes[num_sommet]["position"][1]))
            else :
                G.nodes[num_sommet]["position"] = (min(position_prec_groupe[1], G.nodes[num_sommet]["position"][0]), max(position_prec[1], G.nodes[num_sommet]["position"][1]))
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
def ajout_aretes_artificielles(G, version):
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
        if version == "regroupement" :
            G.add_node(compteur, position=(-1,-1), type=-1, poids=G.nodes[elt]["poids"], chaine=G.nodes[elt]["chaine"], num_ch=G.nodes[elt]["num_ch"], nb_non_can=0)
        else:
            G.add_node(compteur, position=(-1,-1), type=-1, poids=G.nodes[elt]["poids"], chaine=G.nodes[elt]["chaine"], num_ch=G.nodes[elt]["num_ch"])            
        G.add_edge(compteur, elt, label='0', near=None)
        G.add_edge(elt, compteur, label='0', near=None)
    
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
                G.add_edge(noeud, le_bon_voisin, label='B53', near=False)
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

''' renvoie le voisin du sommet passe en parametre, lie a lui par une liaison canonique, renvoie -1 s'il n'existe pas '''
def nb_canoniques(numero_sommet, graphe_mol): ## normalement, pas possible qu'il y en ait plus qu'une
    nb_can = 0
    voisin_can = -1
#     print("ramous")
    if numero_sommet in graphe_mol.nodes() :
        for voisin in graphe_mol[numero_sommet] :
    #         print(voisin)
            if graphe_mol.edges[numero_sommet, voisin]["label"] == 'CWW' and \
                                ((graphe_mol.nodes[numero_sommet]["nt"] == 'A' and graphe_mol.nodes[voisin]["nt"] == 'U') or \
                                (graphe_mol.nodes[numero_sommet]["nt"] == 'U' and graphe_mol.nodes[voisin]["nt"] == 'A') or \
                                (graphe_mol.nodes[numero_sommet]["nt"] == 'C' and graphe_mol.nodes[voisin]["nt"] == 'G') or \
                                (graphe_mol.nodes[numero_sommet]["nt"] == 'G' and graphe_mol.nodes[voisin]["nt"] == 'C') or \
                                (graphe_mol.nodes[numero_sommet]["nt"] == 'G' and graphe_mol.nodes[voisin]["nt"] == 'U') or \
                                (graphe_mol.nodes[numero_sommet]["nt"] == 'U' and graphe_mol.nodes[voisin]["nt"] == 'G'))  : 
                nb_can += 1
                voisin_can = voisin
            

    if nb_can > 1 :
        print("bizarre")
    return voisin_can

''' renvoie la liste des voisins du sommet passe en parametre, lies a lui par une liaison non canonique, renvoie une liste vide s'il n'y en a pas
(renvoie des doublets : (numero du voisin, type de liaison non canonique) '''
def nb_non_canoniques(numero_sommet, graphe_mol):
    nb_non_can = 0
    voisin_non_can = []
    type_non_can = []
#     print(type(graphe_mol))
#     print(type(numero_sommet))
    if numero_sommet in graphe_mol.nodes() :
        for voisin in graphe_mol[numero_sommet] :
            if graphe_mol.edges[numero_sommet, voisin]["label"] != 'B53' and \
            not (graphe_mol.edges[numero_sommet, voisin]["label"] == 'CWW' and 
                  ((graphe_mol.nodes[numero_sommet]["nt"] == 'A' and graphe_mol.nodes[voisin]["nt"] == 'U') or \
                  (graphe_mol.nodes[numero_sommet]["nt"] == 'U' and graphe_mol.nodes[voisin]["nt"] == 'A') or \
                  (graphe_mol.nodes[numero_sommet]["nt"] == 'C' and graphe_mol.nodes[voisin]["nt"] == 'G') or \
                  (graphe_mol.nodes[numero_sommet]["nt"] == 'G' and graphe_mol.nodes[voisin]["nt"] == 'C') or \
                    (graphe_mol.nodes[numero_sommet]["nt"] == 'G' and graphe_mol.nodes[voisin]["nt"] == 'U') or \
                   (graphe_mol.nodes[numero_sommet]["nt"] == 'U' and graphe_mol.nodes[voisin]["nt"] == 'G'))):
    
    #         if graphe_mol.edges[numero_sommet, voisin]["long_range"] == False and \
    #         graphe_mol.edges[numero_sommet, voisin]["label"] != 'B53' and \
    #         graphe_mol.edges[numero_sommet, voisin]["label"] == 'CWW' and not \
    #               ((graphe_mol.nodes[numero_sommet]["nt"] == 'A' and graphe_mol.nodes[voisin]["nt"] == 'U') or \
    #               (graphe_mol.nodes[numero_sommet]["nt"] == 'U' and graphe_mol.nodes[voisin]["nt"] == 'A') or \
    #               (graphe_mol.nodes[numero_sommet]["nt"] == 'C' and graphe_mol.nodes[voisin]["nt"] == 'G') or \
    #               (graphe_mol.nodes[numero_sommet]["nt"] == 'G' and graphe_mol.nodes[voisin]["nt"] == 'C') or \
    #                 (graphe_mol.nodes[numero_sommet]["nt"] == 'G' and graphe_mol.nodes[voisin]["nt"] == 'U') or \
    #                (graphe_mol.nodes[numero_sommet]["nt"] == 'U' and graphe_mol.nodes[voisin]["nt"] == 'G')):
                        nb_non_can += 1
                        voisin_non_can.append(voisin)
                        type_non_can.append(graphe_mol.edges[numero_sommet, voisin]["label"])
           
    return voisin_non_can, type_non_can

''' renvoie le nombre de liaisons non cov qu'a le sommet passe en parametre
dans le graphe de structure de depart '''
def nb_non_cov(numero_sommet, graphe_mol):
    nb_non_cov = 0
    if numero_sommet in graphe_mol.nodes() :
        for voisin in graphe_mol[numero_sommet] :
            if graphe_mol.edges[numero_sommet, voisin]["label"] != 'B53' :
                nb_non_cov += 1
    return nb_non_cov

''' regroupe dans le graphe d'extension les liaisons non canoniques se trouvant a l'interieur d'helices de liaisons canoniques,
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
                if len(voisin_non_cov) == 1 and G.nodes[voisin_non_cov[0][0]]["type"] in [None, -1]  : ## on a mnt une seule liaison qui nest pas long range donc cest une liaison non canonique
                    if len(non_cov(voisin_non_cov[0][0], G)) == 1 : ## le voisin du sommet na pas dautre liaison non covalente
                        print("candidat : %s"%noeud)
                        ##Eligibilite du voisinage pour un regroupement##
                        groupe1 = []
                        nb_non_can = 1
                         
                        voisin_cov = voisin_cov_sens(noeud, G)
                        vieux_voisin = (-1,-1)
                        vieux_voisin_non_cov = non_cov_graphe_base((G.nodes[noeud]["num_ch"], G.nodes[noeud]["position"][0]), graphe_mol)
                        while voisin_cov != vieux_voisin and voisin_cov != -1 :
                            if (G.nodes[voisin_cov[0]]["type"] == 1 and vieux_voisin == (-1,-1)) or (G.nodes[voisin_cov[0]]["type"] in [1,3] and vieux_voisin != (-1,-1))  : ##le premier nt rencontre doit etre implique dans une liaison canonique, pour les autres, ce peut etre une liaison can ou non can
                                    vieux_voisin = voisin_cov
                                    voisin_non_cov = non_cov(voisin_cov[0], G)
                                     
                                    voisin_voisin_non_cov = non_cov_graphe_base((G.nodes[voisin_cov[0]]["num_ch"], G.nodes[voisin_cov[0]]["position"][0]), graphe_mol)
                                    print(voisin_voisin_non_cov)
                                    if len(voisin_non_cov) == 1 and G.nodes[voisin_non_cov[0][0]]["type"] in [None, -1]\
                                    and len(non_cov(voisin_non_cov[0][0], G)) == 1 \
                                    and (vieux_voisin_non_cov == -1 or (abs(voisin_voisin_non_cov[0][1] - vieux_voisin_non_cov[0][1]) < 10 and abs(voisin_voisin_non_cov[0][1] - vieux_voisin_non_cov[0][1]) > 0) ) :
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
                            voisin = non_cov_graphe_base((G.nodes[groupe1[len(groupe1)-1][0]]["num_ch"], G.nodes[groupe1[len(groupe1)-1][0]]["position"][0]), graphe_mol)
                            pos_suivant = G.nodes[groupe1[len(groupe1)-1][0]]["position"][1] + 1
                            num_ch_suivant = G.nodes[groupe1[len(groupe1)-1][0]]["num_ch"]
                            can_new = nb_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            non_can_new, type_non_can = nb_non_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            sommet_non_cov = nb_non_cov((num_ch_suivant, pos_suivant), graphe_mol)
                             
                            ok = False
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_new != -1 and can_new[1] in d["position"] and d["type"] != None :
                                    dedans = True
                             
                            if can_new != -1 and len(non_can_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques(can_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques(can_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_new, graphe_mol) == 1)) and (abs(voisin[0][1] - can_new[1]) < 10 and abs(voisin[0][1] - can_new[1]) > 0 ) :
                                    ok = True
                            if not ok :
                                groupe1.remove(groupe1[len(groupe1)-1])
                                nb_non_can -= 1   
                                 
                        if len(groupe1) == 0 and voisin_cov_sens(noeud, G) == -1 :
                            ## si la liaison non can est en bout de chaine, il faut quand meme verifier le nt d apres
                            voisin = non_cov_graphe_base((G.nodes[noeud]["num_ch"], G.nodes[noeud]["position"][0]), graphe_mol)
                            pos_suivant = G.nodes[noeud]["position"][1] +1
                            num_ch_suivant = G.nodes[noeud]["num_ch"]
                            can_new = nb_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            print(can_new)
                            non_can_new, type_non_can = nb_non_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            sommet_non_cov = nb_non_cov((num_ch_suivant, pos_suivant), graphe_mol)
                             
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_new != -1 and can_new[1] in d["position"] and d["type"] != None :
                                    dedans = True
                                     
                            if can_new != -1 and len(non_can_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques(can_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques(can_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_new, graphe_mol) == 1)) and (abs(voisin[0][1] - can_new[1]) < 10 and abs(voisin[0][1] - can_new[1]) > 0) :
    
                                        groupe1.append(0)
                         
                         
                        groupe2 = []
                        voisin_cov = voisin_cov_nonsens(noeud, G)
                        vieux_voisin = (-1,-1)
                        vieux_voisin_non_cov = non_cov_graphe_base((G.nodes[noeud]["num_ch"], G.nodes[noeud]["position"][0]), graphe_mol)
                        while voisin_cov != vieux_voisin and voisin_cov != -1 :
                             
                            if (G.nodes[voisin_cov[0]]["type"] == 1 and vieux_voisin == (-1,-1)) or (G.nodes[voisin_cov[0]]["type"] in [1,3] and vieux_voisin != (-1,-1))  : ##le premier nt rencontre doit etre implique dans une liaison canonique, pour les autres, ce peut etre une liaison can ou non can
                                vieux_voisin = voisin_cov
                                voisin_non_cov = non_cov(voisin_cov[0], G)
                             
                                 
                                voisin_voisin_non_cov = non_cov_graphe_base((G.nodes[voisin_cov[0]]["num_ch"], G.nodes[voisin_cov[0]]["position"][len(G.nodes[voisin_cov[0]]["position"])-1]), graphe_mol)
                                 
                                if len(voisin_non_cov) == 1 and  G.nodes[voisin_non_cov[0][0]]["type"] in [None, -1] \
                                and len(non_cov(voisin_non_cov[0][0], G)) == 1 \
                                and (vieux_voisin_non_cov == -1 or (abs(voisin_voisin_non_cov[0][1] - vieux_voisin_non_cov[0][1]) < 10 and abs(voisin_voisin_non_cov[0][1] - vieux_voisin_non_cov[0][1]) > 0) ) :
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
                            voisin = non_cov_graphe_base((G.nodes[groupe2[len(groupe2)-1][0]]["num_ch"], G.nodes[groupe2[len(groupe2)-1][0]]["position"][0]), graphe_mol)
                            pos_suivant = G.nodes[groupe2[len(groupe2)-1][0]]["position"][0] -1
                            num_ch_suivant = G.nodes[groupe2[len(groupe2)-1][0]]["num_ch"]
                            can_new = nb_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            non_can_new, type_non_can = nb_non_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            sommet_non_cov = nb_non_cov((num_ch_suivant, pos_suivant), graphe_mol)
                             
                            ok = False
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_new != -1 and can_new[1] in d["position"] and d["type"] != None :
                                    dedans = True
                             
                            if can_new != -1 and len(non_can_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques(can_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques(can_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_new, graphe_mol) == 1)) and (abs(voisin[0][1] - can_new[1]) < 10 and abs(voisin[0][1] - can_new[1]) > 0) :
                                    ok = True
                                 
                            if not ok :
                                groupe2.remove(groupe2[len(groupe2)-1])
                                nb_non_can -= 1 
                         
                        if len(groupe2) == 0 and voisin_cov_nonsens(noeud, G) == -1 :
                            print("gros rat")
                            ## si la liaison non can est en bout de chaine, il faut quand meme verifier le nt d apres
                            voisin = non_cov_graphe_base((G.nodes[noeud]["num_ch"],G.nodes[noeud]["position"][0]), graphe_mol)
                            pos_suivant = G.nodes[noeud]["position"][0] -1
                            num_ch_suivant = G.nodes[noeud]["num_ch"]
                            can_new = nb_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            non_can_new, type_non_can = nb_non_canoniques((num_ch_suivant, pos_suivant), graphe_mol)
                            sommet_non_cov = nb_non_cov((num_ch_suivant, pos_suivant), graphe_mol)
                             
                            print(G.nodes[noeud]["position"][0])
                            print(voisin)
                            print(can_new)
                            print(non_can_new)
                            print(sommet_non_cov)
                            #print(nb_non_cov(can_short_range_new, graphe_mol))
                            print()
                             
                            dedans = False
                            for n, d in G.nodes(data = True) :
                                if can_new != -1 and can_new[1] in d["position"] and d["type"] != None :
                                    dedans = True
                             
                            if can_new != -1 and len(non_can_new) == 0 and sommet_non_cov == 1 :
                                if (not dedans or (nb_canoniques(can_new, graphe_mol) != -1  \
                                and len(nb_non_canoniques(can_new, graphe_mol)[0]) == 0 and \
                                nb_non_cov(can_new, graphe_mol) == 1)) and (abs(voisin[0][1] - can_new[1]) < 10 and abs(voisin[0][1] - can_new[1]) > 0) :
 
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
                G.add_edge(noeud, suivant[0], label='B53', near=False)
        if len(groupe2) > 0 and groupe2[0] != 0:
            prec = voisin_cov_nonsens(groupe2[len(groupe2)-1][0], G)
            pos_2 = G.nodes[groupe2[len(groupe2)-1][0]]["position"][0]
            if prec != -1 :
                G.add_edge(prec[0], noeud, label='B53', near=False)
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

''' correspondance des nucleotides de l'occurrence (stockes dans graphe_motif) avec les numeros des sommets du graphe d'extension (stockes dans G) pour pouvoir ensuite etendre aux voisins
cas generique pour tous les motifs '''
def stockage_motif_generique(G, graphe_motif, positions_ajoutees, graphe_modele):
    dico_correspondance = {}
    
    for i in range(1, graphe_modele.number_of_nodes()+1) :
        dico_correspondance.update({i : -1})
    
    diff = True
    while -1 in dico_correspondance.values() and diff :
        diff = False
        for noeud_motif in graphe_motif.nodes() :
            if noeud_motif not in dico_correspondance.values() :
                nb_corres = 0
                noeud_corres = -1
                for noeud_modele in graphe_modele.nodes() :
                    if dico_correspondance[noeud_modele] == -1 :
                        nb_voisins_idem = 0
                       
                        for voisin_motif in graphe_motif[noeud_motif] :
                            existe = False
                            voisin_voisin_pas_bon = True
                            for voisin_modele in graphe_modele[noeud_modele] :
                                for edge in graphe_modele[noeud_modele][voisin_modele] :
                                    if graphe_motif.edges[noeud_motif, voisin_motif]["label"] == graphe_modele[noeud_modele][voisin_modele][edge]["label"] :
                                        existe = True
                            
                                
                                        nb_voisins_voisin_idem = 0
                                        for voisin_voisin_motif in graphe_motif[voisin_motif] :
                                            existe_voisin = False
                                            for voisin_voisin_modele in graphe_modele[voisin_modele] :
                                                for edge in graphe_modele[voisin_modele][voisin_voisin_modele] :
                                                    if graphe_motif.edges[voisin_motif, voisin_voisin_motif]["label"] == graphe_modele[voisin_modele][voisin_voisin_modele][edge]["label"] :
                                                        existe_voisin = True
                                            if existe_voisin :
                                                nb_voisins_voisin_idem += 1
                                                
                                        if nb_voisins_voisin_idem == len(graphe_motif[voisin_motif]) and nb_voisins_voisin_idem == len(graphe_modele[voisin_modele]) :
                                            voisin_voisin_pas_bon = False
                                    
                            if existe and not voisin_voisin_pas_bon:
                                nb_voisins_idem += 1
#                         print(nb_voisins_idem)
#                         print(len(graphe_motif[noeud_motif]))
#                         print(len(graphe_modele[noeud_modele]))
#                         print(voisin_voisin_pas_bon)
                        if nb_voisins_idem == len(graphe_motif[noeud_motif]) and nb_voisins_idem == len(graphe_modele[noeud_modele]) and not voisin_voisin_pas_bon:
                            nb_corres += 1 
                            noeud_corres = noeud_modele
                        
                if nb_corres == 1 :
                    dico_correspondance[noeud_corres] = noeud_motif
                    diff = True
                elif nb_corres == 0 :
                    print(noeud_motif)
                    print(graphe_motif[noeud_motif])
                    print("bizarre")
                    exit()
        print(dico_correspondance)
        
    while -1 in dico_correspondance.values() :
        for noeud_modele in graphe_modele.nodes() :
            if dico_correspondance[noeud_modele] != -1 :
                for voisin_modele in graphe_modele[noeud_modele] :
                    if dico_correspondance[voisin_modele] == -1 :
                        for edge in graphe_modele[noeud_modele][voisin_modele] :
                            for voisin_motif in graphe_motif[dico_correspondance[noeud_modele]] :
                                if graphe_motif.edges[dico_correspondance[noeud_modele], voisin_motif]["label"] == graphe_modele[noeud_modele][voisin_modele][edge]["label"] :
                                    dico_correspondance[voisin_modele] = voisin_motif 
        
        
        print(dico_correspondance)                  

    for i in range(1, graphe_modele.number_of_nodes()+1) :
        G.add_node(i, position = (dico_correspondance[i][1],dico_correspondance[i][1]), type = graphe_modele.nodes[i]["type"], poids=graphe_modele.nodes[i]["poids"], chaine=list(graphe_modele.nodes[i]["chaine"]), num_ch = dico_correspondance[i][0])
        positions_ajoutees.append(dico_correspondance[i])
        
    for u,v,data in graphe_motif.edges(data=True) :
        if u in dico_correspondance.values() and v in dico_correspondance.values() :
            G.add_edge(list(dico_correspondance.keys())[list(dico_correspondance.values()).index(u)], list(dico_correspondance.keys())[list(dico_correspondance.values()).index(v)], **data )
            
    return G, positions_ajoutees
                   
''' correspondance des nucleotides de l'occurrence (stockes dans graphe_motif) avec les numeros des sommets du graphe d'extension (stockes dans G) pour pouvoir ensuite etendre aux voisins
cas du A-minor'''
def stockage_motif(G, graphe_motif, positions_ajoutees):
    dico_correspondance = {1:-1, 2:-1, 3:-1, 4:-1, 5:-1}
    for noeud in graphe_motif.nodes() :
        b53 = False
        cww = False
        css = False
        tss = False
        for voisin in graphe_motif[noeud] :
            if graphe_motif.edges[noeud, voisin]["label"] == 'B53' : 
                b53 = True
            if graphe_motif.edges[noeud, voisin]["label"] == 'CWW' : 
                cww = True
            if graphe_motif.edges[noeud, voisin]["label"] == 'CSS' : 
                css = True
            if graphe_motif.edges[noeud, voisin]["label"] == 'TSS' : 
                tss = True
                
        if tss and css and not b53 and not cww :
            dico_correspondance[1] = noeud
        elif css and cww and b53 and not tss :
            dico_correspondance[2] = noeud
        elif css and b53 and not cww and not tss :
            dico_correspondance[3] = noeud
        elif css and cww and not b53 and not tss :
            dico_correspondance[4] = noeud
        elif tss and cww and not b53 and not css :
            dico_correspondance[5] = noeud
            
    for i in range(1, 6) :
        G.add_node(i, position = (dico_correspondance[i][1],dico_correspondance[i][1]), type = i+10, poids=1, chaine=[i], num_ch = dico_correspondance[i][0])
        positions_ajoutees.append(dico_correspondance[i])
        
    for u,v,data in graphe_motif.edges(data=True) :
        if u in dico_correspondance.values() and v in dico_correspondance.values() :
            G.add_edge(list(dico_correspondance.keys())[list(dico_correspondance.values()).index(u)], list(dico_correspondance.keys())[list(dico_correspondance.values()).index(v)], **data  )
            
    return G, positions_ajoutees
            
''' 10/03/20 !! :) '''
def obtenir_extension_un_elt(elt, graphe,taille_ext):
    with open("all_aminor.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        all_aminor = mon_depickler.load()
        print(elt)
        liste_pb = []
        for cle in all_aminor.keys() :
            
            if cle == elt[0] :
                    print(cle)
                    compteur_nb = 0
                    compteur_nb_2 = 0
                    compter_graphe = 1
                    for graph_motif in all_aminor[cle] :
                        if compter_graphe == int(elt[1]) :
                            print(compter_graphe)
#                             for noeud, data in graphe.nodes(data=True) :
#                                 print(str(noeud) + ' ' + str(data) + "\n")
                        
                            positions_chaines = [[],[],[],[]]
                             
                            G = nx.MultiDiGraph()
                            i = 1
                            positions_ajoutees = []
                            print(graph_motif.edges.data())
                            #G, positions_ajoutees = stockage_motif_generique(G, graph_motif, positions_ajoutees, motif_7left)
                            G, positions_ajoutees = stockage_motif(G, graph_motif, positions_ajoutees)

                            
                            print(graph_motif.edges.data())
                            print(G.edges.data())
                            print(positions_ajoutees)
                            
                            #compteur = motif_7left.number_of_nodes()+1
                            compteur = 6
                            
                            print(G.nodes.data())
                            print(G.edges.data())
        #                     print(G.nodes.data())
                            #print(G.edges.data())
                            #fichier.write(str(graphes[('1U9S', 'A')].nodes.data())+"\n")
                            #fichier.write(str(graphes[('1U9S', 'A')].edges.data()))
                            #nom_cle = (occ["num_PDB"], occ["num_ch"])
                            #print(str(nom_cle) + str(occ["num_motif"]) + str(occ["num_occ"])) 
                             
        #                     if nom_cle == ('4V88', 'A6') and occ["num_motif"] == 17 and occ["num_occ"] == 55 :
        #                         print(nom_cle)
        #                         print(occ["num_motif"])
        #                         print(occ["num_occ"])
                            
                            num_ch_1 = []
                            num_ch_2 = []
                            for noeud, data in G.nodes(data=True) :
                                print(data)
                                #if data["chaine"] == [1,3] :#or data["chaine"] == [3] :
                                if data["chaine"] == [1] or data["chaine"] == [3] :
                                    num_ch_1.insert(0, noeud)
                                if data["chaine"] == [2] or data["chaine"] == [4] :
                                    num_ch_2.append(noeud)
                            compteur_tige = num_ch_2[0]   
                            
                            print(num_ch_1)   
                            print(num_ch_2)
                            ##Garder toutes les positions de la sequence
                            if len(num_ch_1) == 2 :
                                if G.node[num_ch_1[0]]["position"][0] - G.node[num_ch_1[1]]["position"][0] < 0 : ## ordre de haut en bas
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_1[0]]["position"][0]-i > 0 and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i) in graphe.nodes():
                                            positions_chaines[0].append((G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i))
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_1[1]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_1[1]]["num_ch"], G.node[num_ch_1[1]]["position"][0]+i) in graphe.nodes() :
                                            positions_chaines[2].append((G.node[num_ch_1[1]]["num_ch"],G.node[num_ch_1[1]]["position"][0]+i))
                                else : ## ordre de bas en haut
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_1[0]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]+i) in graphe.nodes() :
                                            positions_chaines[0].append((G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]+i))
                                    for i in range(taille_ext) :
                                        if G.node[3]["position"][0]-i > 0 and (G.node[3]["num_ch"], G.node[3]["position"][0]-i) in graphe.nodes():
                                            positions_chaines[2].append((G.node[3]["num_ch"], G.node[3]["position"][0]-i))
                            else :
                                for i in range(taille_ext) :
                                    if G.node[num_ch_1[0]]["position"][0]-i > 0 and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i) in graphe.nodes():
                                        positions_chaines[0].append((G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i))
                                for i in range(taille_ext) :
                                    if G.node[num_ch_1[0]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]+i) in graphe.nodes() :
                                        positions_chaines[2].append((G.node[num_ch_1[0]]["num_ch"],G.node[num_ch_1[0]]["position"][0]+i))
          
                                        
                            if G.node[num_ch_2[0]]["position"][0] - G.node[num_ch_2[1]]["position"][0] > 0 : ## ordre de bas en haut
                                for i in range(taille_ext) :
                                    if G.node[num_ch_2[1]]["position"][0]-i > 0 and (G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]-i) in graphe.nodes() :
                                        positions_chaines[3].append((G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]-i))
                                for i in range(taille_ext) :
                                    if G.node[num_ch_2[0]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_2[0]]["num_ch"], G.node[num_ch_2[0]]["position"][0]+i) in graphe.nodes() :
                                        positions_chaines[1].append((G.node[num_ch_2[0]]["num_ch"], G.node[num_ch_2[0]]["position"][0]+i))
                            else : ## ordre de bas en haut
                                print("gros rat")
                                for i in range(taille_ext) :
                                    if G.node[num_ch_2[1]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]+i) in graphe.nodes() :
                                        positions_chaines[3].append((G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]+i))
                                for i in range(taille_ext) :
                                    if G.node[num_ch_2[0]]["position"][0]-i > 0 and  (G.node[num_ch_2[0]]["num_ch"], G.node[2]["position"][0]-i) in graphe.nodes() :
                                        positions_chaines[1].append((G.node[num_ch_2[0]]["num_ch"], G.node[num_ch_2[0]]["position"][0]-i))
                             
                            print(positions_chaines)  
                            ##Connaitre l'ordre des positions sur la sequence
                            int_tige = 0
                            if G.node[num_ch_2[0]]["position"][0] - G.node[num_ch_2[1]]["position"][0] < 0 : ## ordre de haut en bas
                                int_tige = 1
                            else : ## ordre de bas en haut
                                int_tige = -1
                            #compteur_tige2 = G.node[5]["position"][0] + int_tige
                             
                            ##Initialisation de la boucle 1
                            #===================================================================
                            print(compteur_tige) 
                            G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5)
                            print(G.nodes.data())
                            print(G.edges.data())
                            #print(G.nodes.data())
                            #print(G.edges.data())
                            #print(G.nodes.data())
                            #print(positions_ajoutees)
                            #if nom_cle == ('1FJG', 'A') :
                            #    print(G.nodes.data())
                            #    print(G.edges.data())
                             
                            compteur = G.number_of_nodes()+1
                                 
                            #compteur_tige2 = G.node[5]["position"][0] - int_tige
                            compteur_tige = num_ch_2[1]
                             
        #                     print(nom_cle)
        #                     print(occ["num_motif"])
        #                     print(occ["num_occ"])
                            print(compteur_tige) 
                            G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5)  
                            print(G.nodes.data())
                            print(G.edges.data())
                            #print(G.nodes.data())
                            #print(G.edges.data())
                            
                            compteur_tige = num_ch_1[1]
                            int_tige = -1
                            compteur = G.number_of_nodes()+1
#                                 if G.node[num_ch_1[0]]["position"][0] - G.node[num_ch_1[1]]["position"][0] < 0 : ## ordre de haut en bas
#                                     int_tige = 1
#                                 else : ## ordre de bas en haut
#                                     int_tige = -1
                              
                             
                            print(compteur_tige) 
                            G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5) 
                            print(G.nodes.data())
                            print(G.edges.data()) 
         
                            compteur_tige = num_ch_1[0]
                            int_tige = 1
                            compteur = G.number_of_nodes()+1
                             
                            print(compteur_tige) 
                            G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5) 
                             
                            print(G.nodes.data())
                            print(G.edges.data())
                            
                            G = ajout_liaisons_B53_qui_manquent(G)
                             
                            print(G.nodes.data())
                            print(G.edges.data()) 
                            G = ajout_aretes_artificielles(G, "non_regroupement")
                            G = regroupement_liaisons_short_range(G, graphe)
                            
                            G = ajout_aretes_artificielles(G, "regroupement")
    #                         print(G.nodes.data())
    #                         print(G.edges.data())
    #                         print(G.number_of_edges())
                           
                            #print(occ)
        #                     print(G.nodes.data())
        #                     print(G.edges.data())
                            
        #                     compteur = G.number_of_nodes()+1
        #                     compteur_tige = 1
                            #G = tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees) 
    #                         for noeud, data in G.nodes(data=True) :
    #                             if noeud not in G2.nodes or data["type"] != G2.nodes[noeud]["type"]:
    #                                 print(noeud)
                            print(G.nodes.data())
                            print(G.edges.data())
                            
                            for noeud, data in G.nodes(data=True) :
                                print(noeud)
                                print(data)
                                print(G[noeud])
                                
                        compter_graphe += 1
        return G

            
''' new data 02/09/19 '''
def obtenir_extension_new_data(taille_ext):  
    
    motif_6A = nx.MultiDiGraph()
    for i in range(1,12) :
        motif_6A.add_node(i, type=i+10, poids=1, chaine=[], position=())
    motif_6A.add_edge(1,2, label="B53")
    motif_6A.add_edge(1,9, label="TSS")
    motif_6A.add_edge(9,1, label="TSS")
    motif_6A.add_edge(2,3, label="CWW")
    motif_6A.add_edge(3,2, label="CWW")
    motif_6A.add_edge(3,4, label="B53")
    motif_6A.add_edge(4,5, label="B53")
    motif_6A.add_edge(5,9, label="CWW")
    motif_6A.add_edge(9,5, label="CWW")
    motif_6A.add_edge(4,10, label="CWW")
    motif_6A.add_edge(10,4, label="CWW")
    motif_6A.add_edge(8,9, label="B53")
    motif_6A.add_edge(9,10, label="B53")
    motif_6A.add_edge(10,11, label="B53")
    motif_6A.add_edge(8,7, label="THW")
    motif_6A.add_edge(7,8, label="TWH")
    motif_6A.add_edge(6,11, label="CWW")
    motif_6A.add_edge(11,6, label="CWW")
    motif_6A.add_edge(6,7, label="B53")
    
    motif_6A.nodes[3]["chaine"] = [3]
    motif_6A.nodes[4]["chaine"] = [1,3]
    motif_6A.nodes[5]["chaine"] = [1]
    motif_6A.nodes[1]["chaine"] = [5]
    motif_6A.nodes[2]["chaine"] = [6]
    motif_6A.nodes[6]["chaine"] = [8]
    motif_6A.nodes[7]["chaine"] = [7]
    motif_6A.nodes[8]["chaine"] = [2]
    motif_6A.nodes[9]["chaine"] = [2,4]
    motif_6A.nodes[10]["chaine"] = [2,4]
    motif_6A.nodes[11]["chaine"] = [4]
    
    motif_7left = nx.MultiDiGraph()
    for i in range(1,6) :
        motif_7left.add_node(i, type=i+10, poids=1, chaine=[], position=())
    motif_7left.add_edge(1,2, label="B53")
    motif_7left.add_edge(1,3, label="CWW")
    motif_7left.add_edge(3,1, label="CWW")
    motif_7left.add_edge(2,5, label="TSW")
    motif_7left.add_edge(5,2, label="TWS")
    motif_7left.add_edge(3,4, label="CSW")
    motif_7left.add_edge(4,3, label="CWS")
    motif_7left.add_edge(4,5, label="B53")
    
    motif_7left.nodes[1]["chaine"] = [5]
    motif_7left.nodes[2]["chaine"] = [6]
    motif_7left.nodes[3]["chaine"] = [1,3]
    motif_7left.nodes[4]["chaine"] = [2]
    motif_7left.nodes[5]["chaine"] = [4]
    
    
    
    
    #with open("fichier_extension_new_data"+str(taille_ext)+"_avec_regroupement_aretes.txt", "w") as fichier :  
    #with open("Graphs_carnaval/7left.pickle", 'rb') as fichier_pickle :
    with open("all_aminor.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            all_aminor = mon_depickler.load()
            
            liste_pb = []
            for cle in all_aminor.keys() :
                #print(cle)
                if cle == "4y4o" :
                    print("petit rat")
                    with open("Graphs/%s.pickle"%cle, 'rb') as fichier_tout :
                        mon_depickler_graphes = pickle.Unpickler(fichier_tout)
                        graphe = mon_depickler_graphes.load()
                        
                        
                        compteur_nb = 0
                        compteur_nb_2 = 0
                        compter_graphe = 1
                        for graph_motif in all_aminor[cle] :
                            if compter_graphe == 23 :
#                             for noeud, data in graphe.nodes(data=True) :
#                                 print(str(noeud) + ' ' + str(data) + "\n")
                            
                                positions_chaines = [[],[],[],[]]
                                 
                                G = nx.MultiDiGraph()
                                i = 1
                                positions_ajoutees = []
                                print(graph_motif.edges.data())
                                #G, positions_ajoutees = stockage_motif_generique(G, graph_motif, positions_ajoutees, motif_7left)
                                G, positions_ajoutees = stockage_motif(G, graph_motif, positions_ajoutees)

                                
                                print(graph_motif.edges.data())
                                print(G.edges.data())
                                print(positions_ajoutees)
                                
                                #compteur = motif_7left.number_of_nodes()+1
                                compteur = 6
                                
                                print(G.nodes.data())
                                print(G.edges.data())
            #                     print(G.nodes.data())
                                #print(G.edges.data())
                                #fichier.write(str(graphes[('1U9S', 'A')].nodes.data())+"\n")
                                #fichier.write(str(graphes[('1U9S', 'A')].edges.data()))
                                #nom_cle = (occ["num_PDB"], occ["num_ch"])
                                #print(str(nom_cle) + str(occ["num_motif"]) + str(occ["num_occ"])) 
                                 
            #                     if nom_cle == ('4V88', 'A6') and occ["num_motif"] == 17 and occ["num_occ"] == 55 :
            #                         print(nom_cle)
            #                         print(occ["num_motif"])
            #                         print(occ["num_occ"])
                                
                                num_ch_1 = []
                                num_ch_2 = []
                                for noeud, data in G.nodes(data=True) :
                                    print(data)
                                    #if data["chaine"] == [1,3] :#or data["chaine"] == [3] :
                                    if data["chaine"] == [1] or data["chaine"] == [3] :
                                        num_ch_1.insert(0, noeud)
                                    if data["chaine"] == [2] or data["chaine"] == [4] :
                                        num_ch_2.append(noeud)
                                compteur_tige = num_ch_2[0]   
                                
                                print(num_ch_1)   
                                print(num_ch_2)
                                ##Garder toutes les positions de la sequence
                                if len(num_ch_1) == 2 :
                                    if G.node[num_ch_1[0]]["position"][0] - G.node[num_ch_1[1]]["position"][0] < 0 : ## ordre de haut en bas
                                        for i in range(taille_ext) :
                                            if G.node[num_ch_1[0]]["position"][0]-i > 0 and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i) in graphe.nodes():
                                                positions_chaines[0].append((G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i))
                                        for i in range(taille_ext) :
                                            if G.node[num_ch_1[1]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_1[1]]["num_ch"], G.node[num_ch_1[1]]["position"][0]+i) in graphe.nodes() :
                                                positions_chaines[2].append((G.node[num_ch_1[1]]["num_ch"],G.node[num_ch_1[1]]["position"][0]+i))
                                    else : ## ordre de bas en haut
                                        for i in range(taille_ext) :
                                            if G.node[num_ch_1[0]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]+i) in graphe.nodes() :
                                                positions_chaines[0].append((G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]+i))
                                        for i in range(taille_ext) :
                                            if G.node[3]["position"][0]-i > 0 and (G.node[3]["num_ch"], G.node[3]["position"][0]-i) in graphe.nodes():
                                                positions_chaines[2].append((G.node[3]["num_ch"], G.node[3]["position"][0]-i))
                                else :
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_1[0]]["position"][0]-i > 0 and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i) in graphe.nodes():
                                            positions_chaines[0].append((G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]-i))
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_1[0]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_1[0]]["num_ch"], G.node[num_ch_1[0]]["position"][0]+i) in graphe.nodes() :
                                            positions_chaines[2].append((G.node[num_ch_1[0]]["num_ch"],G.node[num_ch_1[0]]["position"][0]+i))
              
                                            
                                if G.node[num_ch_2[0]]["position"][0] - G.node[num_ch_2[1]]["position"][0] > 0 : ## ordre de bas en haut
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_2[1]]["position"][0]-i > 0 and (G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]-i) in graphe.nodes() :
                                            positions_chaines[3].append((G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]-i))
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_2[0]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_2[0]]["num_ch"], G.node[num_ch_2[0]]["position"][0]+i) in graphe.nodes() :
                                            positions_chaines[1].append((G.node[num_ch_2[0]]["num_ch"], G.node[num_ch_2[0]]["position"][0]+i))
                                else : ## ordre de bas en haut
                                    print("gros rat")
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_2[1]]["position"][0]+i < graphe.number_of_nodes() and (G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]+i) in graphe.nodes() :
                                            positions_chaines[3].append((G.node[num_ch_2[1]]["num_ch"], G.node[num_ch_2[1]]["position"][0]+i))
                                    for i in range(taille_ext) :
                                        if G.node[num_ch_2[0]]["position"][0]-i > 0 and  (G.node[num_ch_2[0]]["num_ch"], G.node[2]["position"][0]-i) in graphe.nodes() :
                                            positions_chaines[1].append((G.node[num_ch_2[0]]["num_ch"], G.node[num_ch_2[0]]["position"][0]-i))
                                 
                                print(positions_chaines)  
                                ##Connaitre l'ordre des positions sur la sequence
                                int_tige = 0
                                if G.node[num_ch_2[0]]["position"][0] - G.node[num_ch_2[1]]["position"][0] < 0 : ## ordre de haut en bas
                                    int_tige = 1
                                else : ## ordre de bas en haut
                                    int_tige = -1
                                #compteur_tige2 = G.node[5]["position"][0] + int_tige
                                 
                                ##Initialisation de la boucle 1
                                #===================================================================
                                print(compteur_tige) 
                                G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5)
                                print(G.nodes.data())
                                print(G.edges.data())
                                #print(G.nodes.data())
                                #print(G.edges.data())
                                #print(G.nodes.data())
                                #print(positions_ajoutees)
                                #if nom_cle == ('1FJG', 'A') :
                                #    print(G.nodes.data())
                                #    print(G.edges.data())
                                 
                                compteur = G.number_of_nodes()+1
                                     
                                #compteur_tige2 = G.node[5]["position"][0] - int_tige
                                compteur_tige = num_ch_2[1]
                                 
            #                     print(nom_cle)
            #                     print(occ["num_motif"])
            #                     print(occ["num_occ"])
                                print(compteur_tige) 
                                G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5)  
                                print(G.nodes.data())
                                print(G.edges.data())
                                #print(G.nodes.data())
                                #print(G.edges.data())
                                
                                compteur_tige = num_ch_1[1]
                                int_tige = -1
                                compteur = G.number_of_nodes()+1
#                                 if G.node[num_ch_1[0]]["position"][0] - G.node[num_ch_1[1]]["position"][0] < 0 : ## ordre de haut en bas
#                                     int_tige = 1
#                                 else : ## ordre de bas en haut
#                                     int_tige = -1
                                  
                                 
                                print(compteur_tige) 
                                G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5) 
                                print(G.nodes.data())
                                print(G.edges.data()) 
             
                                compteur_tige = num_ch_1[0]
                                int_tige = 1
                                compteur = G.number_of_nodes()+1
                                 
                                print(compteur_tige) 
                                G = extension_tige_new_data(G, graphe, compteur, compteur_tige, positions_ajoutees, int_tige, taille_ext, positions_chaines, num_ch_1, num_ch_2, 5) 
                                 
                                print(G.nodes.data())
                                print(G.edges.data())
                                
                                G = ajout_liaisons_B53_qui_manquent(G)
                                 
                                print(G.nodes.data())
                                print(G.edges.data()) 
                                G = ajout_aretes_artificielles(G, "non_regroupement")
                                G = regroupement_liaisons_short_range(G, graphe)
                                
                                G = ajout_aretes_artificielles(G, "regroupement")
        #                         print(G.nodes.data())
        #                         print(G.edges.data())
        #                         print(G.number_of_edges())
                               
                                #print(occ)
            #                     print(G.nodes.data())
            #                     print(G.edges.data())
                                
            #                     compteur = G.number_of_nodes()+1
            #                     compteur_tige = 1
                                #G = tige(G, graphes, nom_cle, compteur, compteur_tige, positions_ajoutees) 
        #                         for noeud, data in G.nodes(data=True) :
        #                             if noeud not in G2.nodes or data["type"] != G2.nodes[noeud]["type"]:
        #                                 print(noeud)
                                print(G.nodes.data())
                                print(G.edges.data())
                                
                                for noeud, data in G.nodes(data=True) :
                                    print(noeud)
                                    print(data)
                                    print(G[noeud])
                                     
            
#                                 fichier.write(str(cle)+ " " + str(compter_graphe) +"\n")
#                                 for noeud, data in G.nodes(data=True) :
#                                     fichier.write(str(noeud) + " " + str(data) +"\n")
#                                             
#                                 for u,v,data in G.edges(data=True) :
#                                     fichier.write(str((u,v)) + " " + str(data)+"\n")
                                #break     
#                                 with open(NEW_EXTENSION_PATH_TAILLE+"motif_6A/fichier_{}.pickle".format(str(cle) + "_" + str(compter_graphe)), "rb") as fichier :
#                                     mon_depickler = pickle.Unpickler(fichier)
#                                     ancien_g = mon_depickler.load()
#                                     
#                                     if ancien_g.number_of_edges() != G.number_of_edges() :
#                                         liste_pb.append(str(cle) + "_" + str(compter_graphe))
                                        
                                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_{}_3.pickle".format(str(cle) + "_" + str(compter_graphe)), "wb") as fichier_sortie :
                                        mon_pickler = pickle.Pickler(fichier_sortie)
                                        mon_pickler.dump(G)
                                    
                            compter_graphe += 1
                        #break
    
    print(liste_pb)
if __name__ == '__main__':
    obtenir_extension_new_data(4)
         