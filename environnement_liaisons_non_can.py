'''
Created on 20 mai 2019

@author: coline

Environnement des liaisons non canoniques
Recherche de la proportion des liaisons non canoniques faisant partie d'hélices de différentes tailles
(version CaRNAval, non testee sur version toutes données PDB)

'''

import pickle
import os
import csv
from recup_data.constantes import EXTENSION_PATH_TAILLE, EXTENSION_TOUTES_ARETES


def nb_canoniques_short_range(numero_sommet, graphe_mol): ## normalement, pas possible qu'il y en ait plus qu'une
    nb_can = 0
    voisin_can = -1
#     print("ramous")
    for voisin in graphe_mol[numero_sommet] :
#         print(voisin)
        if graphe_mol.edges[numero_sommet, voisin]["label"] == 'CWW' and \
                            graphe_mol.edges[numero_sommet, voisin]["long_range"] == False and \
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

def nb_non_canoniques_short_range(numero_sommet, graphe_mol):
    nb_non_can = 0
    voisin_non_can = []
    type_non_can = []
#     print(type(graphe_mol))
#     print(type(numero_sommet))
    for voisin in graphe_mol[numero_sommet] :
        if graphe_mol.edges[numero_sommet, voisin]["long_range"] == False and \
        graphe_mol.edges[numero_sommet, voisin]["label"] != 'B53' and \
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
                    type_non_can.append( graphe_mol.edges[numero_sommet, voisin]["label"])
           
    return voisin_non_can, type_non_can


def nb_non_cov(numero_sommet, graphe_mol):
    nb_non_cov = 0
    for voisin in graphe_mol[numero_sommet] :
        if graphe_mol.edges[numero_sommet, voisin]["label"] != 'B53' :
            nb_non_cov += 1
    return nb_non_cov
        
        
def environnement_molecule(graphe, graphe_mol):
    tab_can = []
    tab_non_can = []
    couples_vus = []
    tab_type_non_can = []
    tab_nt = []
    type_nt_a_garder = False
    for u, data in graphe.nodes(data=True) :
        if data["type"] == 3 : ## on cherche les liaisons non canoniques de notre extension qui peuvent appartenir a une helice simple sans autre liaison sur les cotes : le sommet incident est forcement de type 3
            voisin_pos_depart = nb_non_canoniques_short_range(graphe.nodes[u]["position"][0], graphe_mol)[0]
            if len(voisin_pos_depart) == 1 and nb_non_cov(graphe.nodes[u]["position"][0], graphe_mol) == 1  : ## notre liaison non canonique doit etre la seule incidente à son sommet (en dehors des liaisons cov)
                dedans = False
                num_voisin = -1
                for noeud, d in graphe.nodes(data = True) :
                    if voisin_pos_depart[0] in d["position"] and d["type"] != None :
                        dedans = True
                        num_voisin = noeud
                if (num_voisin, u) not in couples_vus and (not dedans or len(nb_non_canoniques_short_range(voisin_pos_depart[0], graphe_mol)[0]) == 1 and nb_non_cov(voisin_pos_depart[0], graphe_mol) == 1 ): ## une seule liaison non covalente pour le sommet courant et aussi pour son voisin
                    ## la liaison ne doit pas deja avoir ete traitee (c est possible si elle appartient a une helice avec une autre liaison non canonique)
                    ## le voisin du sommet de la chaine doit etre de type None ou bien ne pas posseder d'autre liaison non cov
                    couples_vus.append((u,num_voisin))
                    print((u, num_voisin))
                    print(couples_vus)
                    print((6,7) in couples_vus)
                    print("numero liaison : %d\n"%(u))
                    pos_depart = graphe.nodes[u]["position"][0]
                    nb_helice = 0
                    nb_non_can = 0
                    type_non_can = [graphe_mol.edges[pos_depart, voisin_pos_depart[0]]["label"]]
                    type_nt = [graphe_mol.nodes[pos_depart]["nt"], graphe_mol.nodes[voisin_pos_depart[0]]["nt"]]
                    i = 1
                    vieux_i = -1
                    can_short_range = -1
                    non_can_short_range = []
                    fini_non_can = False
                    while vieux_i != i and pos_depart+i < graphe_mol.number_of_nodes():
                        vieux_i = i
    #                     print(pos_depart+i)
   
                        can_short_range = nb_canoniques_short_range(pos_depart+i, graphe_mol)
                        non_can_short_range, type_short_range_non_can = nb_non_canoniques_short_range(pos_depart+i, graphe_mol)
                        non_cov = nb_non_cov(pos_depart+i, graphe_mol)
                        
                        
    #                     print(i)
    #                     print(can_short_range)
    #                     print(non_can_short_range)
                        print(can_short_range)  
                        print(non_can_short_range)    
                        print(non_cov) 
                        print(voisin_pos_depart[0]+i)
                        print(voisin_pos_depart[0]-i)
                        if i == 1 : ## le nt suivant immediat sur la chaine doit posseder une liaison canonique et rien d'autre 
                            if can_short_range != -1 and len(non_can_short_range) == 0 and non_cov == 1 :
                                dedans = False
                                for noeud, d in graphe.nodes(data = True) :
                                    if can_short_range in d["position"] and d["type"] != None :
                                        dedans = True
                                ## toujours pareil son voisin par la liaison canonique doit etre de type None ou ne pas posseder d autre liaison non cov
                                if (not dedans or (nb_canoniques_short_range(can_short_range, graphe_mol) != -1 and len(nb_non_canoniques_short_range(can_short_range, graphe_mol)[0]) == 0 and nb_non_cov(can_short_range, graphe_mol) == 1)) and (abs(voisin_pos_depart[0]+i - can_short_range) < 10 or abs(voisin_pos_depart[0]-i - can_short_range) < 10) :
                                    i = i+1
                                    fini_non_can = False ## le dernier nt vu avait une liaison can, et non, non can
                                    print("ramou")
                        else : ## a partir d une distance 2 de la liaison non canonique de depart, on autorise les liaisons non can
                            #idem qu au dessus : cas d une liaison canonique
                            if can_short_range != -1 and len(non_can_short_range) == 0 and non_cov == 1 :
                                dedans = False
                                for noeud, d in graphe.nodes(data = True) :
                                    if can_short_range in d["position"] and d["type"] != None :
                                        dedans = True
                                
                                if (not dedans or (nb_canoniques_short_range(can_short_range, graphe_mol) != -1 and len(nb_non_canoniques_short_range(can_short_range, graphe_mol)[0]) == 0 and nb_non_cov(can_short_range, graphe_mol) == 1))  and (abs(voisin_pos_depart[0]+i - can_short_range) < 10 or abs(voisin_pos_depart[0]-i - can_short_range) < 10) :
                                    i = i+1
                                    fini_non_can = False
                                    print("ramou")
                                    
                            ## cas d une liaison non canonique
                            elif can_short_range == -1 and len(non_can_short_range) == 1 and non_cov == 1 :
                                dedans = False
                                for noeud, d in graphe.nodes(data = True) :
                                    if non_can_short_range[0] in d["position"] and d["type"] != None:
                                        dedans = True
                                if (not dedans or (nb_canoniques_short_range(non_can_short_range[0], graphe_mol) == -1 and len(nb_non_canoniques_short_range(non_can_short_range[0], graphe_mol)[0]) == 1 and nb_non_cov(non_can_short_range[0], graphe_mol) == 1))  and ( abs(voisin_pos_depart[0]+i - non_can_short_range[0]) < 10 or abs(voisin_pos_depart[0]-i - non_can_short_range[0]) < 10) : 
                                    nb_non_can += 1
                                    for elt in type_short_range_non_can :
                                        type_non_can.append(elt)
                                    i = i+1
                                    fini_non_can = True
                        #print(i)
                    if fini_non_can : ## si la derniere liaison est une non canonique, on ne la compte  
                        nb_helice = i-2
                        nb_non_can = nb_non_can -1 
                        type_non_can.remove(type_non_can[len(type_non_can)-1])
                    else :
                        nb_helice = i-1
    #                 print(nb_can)
    #                 print(nb_non_can)
                    if nb_helice > 0 :
                        i = 1
                        vieux_i = -1
                        nb_helice2 = 0
                        while vieux_i != i and pos_depart-i < graphe_mol.number_of_nodes():
                            vieux_i = i
        #                     print(pos_depart+i)
                            can_short_range = nb_canoniques_short_range(pos_depart-i, graphe_mol)
                            non_can_short_range, type_short_range_non_can = nb_non_canoniques_short_range(pos_depart-i, graphe_mol)
                            non_cov = nb_non_cov(pos_depart-i, graphe_mol)
        #                     print(i)
        #                     print(can_short_range)
        #                     print(non_can_short_range)
                            print(can_short_range)   
                            print(non_can_short_range)  
                            if i == 1 :  
                                if can_short_range != -1 and len(non_can_short_range) == 0 and non_cov == 1 :
                                    dedans = False
                                    for noeud, d in graphe.nodes(data = True) :
                                        if can_short_range in d["position"] and d["type"] != None :
                                            dedans = True
                                    
                                    if (not dedans or (nb_canoniques_short_range(can_short_range, graphe_mol) != -1 and len(nb_non_canoniques_short_range(can_short_range, graphe_mol)[0]) == 0 and nb_non_cov(can_short_range, graphe_mol) == 1)) and (abs(voisin_pos_depart[0]+i - can_short_range) < 10 or abs(voisin_pos_depart[0]-i - can_short_range) < 10 ) :
                                        i = i+1
                                        fini_non_can = False
                                        print("ramou")
                            else :
                                if can_short_range != -1 and len(non_can_short_range) == 0 and non_cov == 1 :
                                    dedans = False
                                    for noeud, d in graphe.nodes(data = True) :
                                        if can_short_range in d["position"] and d["type"] != None :
                                            dedans = True
                                    
                                    if (not dedans or (nb_canoniques_short_range(can_short_range, graphe_mol) != -1 and len(nb_non_canoniques_short_range(can_short_range, graphe_mol)[0]) == 0 and nb_non_cov(can_short_range, graphe_mol) == 1)) and (abs(voisin_pos_depart[0]+i - can_short_range) < 10 or abs(voisin_pos_depart[0]-i - can_short_range) < 10 ) :
                                        i = i+1
                                        print("ramou")
                                elif can_short_range == -1 and len(non_can_short_range) == 1 and non_cov == 1 :
                                    dedans = False
                                    for noeud, d in graphe.nodes(data = True) :
                                        if non_can_short_range[0] in d["position"] and d["type"] != None:
                                            dedans = True
                                
                                    if (not dedans or (nb_canoniques_short_range(non_can_short_range[0], graphe_mol) == -1 and len(nb_non_canoniques_short_range(non_can_short_range[0], graphe_mol)[0]) == 1 and nb_non_cov(non_can_short_range[0], graphe_mol) == 1)) and ( abs(voisin_pos_depart[0]+i - non_can_short_range[0]) < 10 or abs(voisin_pos_depart[0]-i - non_can_short_range[0]) < 10) :
                                        nb_non_can += 1
                                        for elt in type_short_range_non_can :
                                            type_non_can.append(elt)
                                        i = i+1
                                        fini_non_can = True
                    
    #                 print(i)
                        if fini_non_can :   
                            nb_helice2 = i-2
                            nb_non_can = nb_non_can -1 
                            type_non_can.remove(type_non_can[len(type_non_can)-1])
                        else :
                            nb_helice2 = i-1
                    
                    if nb_helice > 0 and nb_helice2 > 0 : ## la liaison non canonique de depart doit etre entouree des deux cotes par des liaisons d helice
                        print("petit rat")
                        print(nb_helice)
                        print(nb_non_can)
                        tab_can.append(nb_helice+nb_helice2)
                        tab_non_can.append(nb_non_can)
                        tab_type_non_can.append(type_non_can)
                        type_nt_a_garder = type_nt
    return tab_can, tab_non_can, tab_type_non_can, type_nt_a_garder
                
def test_tige_contigue(graphe, graphe_mol):   
    garder = []         
    for noeud, data in graphe.nodes(data=True) :
        if data["type"] == 1 :
            vieux_voisin_pos = -1
            voisin_pos = -1
            part_id = []
            for i in range(data["position"][0], data["position"][1]+1) :
                vieux_voisin_pos = voisin_pos
                voisin_pos = nb_canoniques_short_range(i, graphe_mol)
                part_id.append(graphe_mol.nodes[i]["part_id"])
                if vieux_voisin_pos != -1 and abs(vieux_voisin_pos - voisin_pos) != 1 :
                    garder.append(noeud)
            diff_part_id = False        
            if len(part_id) > 0 :
                elt_vieux = part_id[0]
                for j in range(1, len(part_id)) :
                    if elt_vieux != part_id[j] :
                        diff_part_id = True
#                         print(part_id[j])
#                         print(elt_vieux)
                    elt_vieux = part_id[j]
                
                    
    return garder, diff_part_id
if __name__ == '__main__':
    with open("graphs_2.92.pickle", 'rb') as fichier_tout :
        mon_depickler_graphes = pickle.Unpickler(fichier_tout)
        graphes = mon_depickler_graphes.load()
        
        dico_val_can = {}
        dico_val_non_can = {}
        tab_val_can = []
        tab_val_non_can = []
        tab_type_non_can = []
        dico_test = {}
        tab_test_2 = []
        tab_type_nt = []
        compteur = 0
        for fic in os.listdir(EXTENSION_PATH_TAILLE%10) :
            if "pickle" in fic and "couples_possibles" not in fic and len(fic.split("_")) == 6 and "_3.pickle" not in fic :
                print("ramousnif")
                with open(EXTENSION_PATH_TAILLE%10+fic, 'rb') as fichier_graphe :
                    mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
                    graphe = mon_depickler_graphe.load()
                    nom_mol = (fic.split("_")[1], fic.split("_")[2])
                    print(fic)
                    print(nom_mol)
                    nb_can, nb_non_can, type_non_can, type_nt = environnement_molecule(graphe, graphes[nom_mol])
                    tab_val_can.append(nb_can)
                    tab_val_non_can.append(nb_non_can)
                    tab_type_non_can.append(type_non_can)
                    tab_type_nt.append(type_nt)
                    
                    dico_val_can.update({fic : nb_can})
                    dico_val_non_can.update({fic : nb_non_can})
                    
                    compteur += 1
                    
                    tige_contigue, diff = test_tige_contigue(graphe, graphes[nom_mol])
                    pos = []
                    for elt in tige_contigue :
                        pos.append(graphe.nodes[elt]["position"])
                    if len(tige_contigue) != 0 :
                        dico_test.update({fic : [tige_contigue, pos]})
#                     if diff  :
#                         tab_test_2.append(fic)
#                         dico_test.update({fic : tige_contigue})
                    
        print(dico_test)
        print(len(dico_test))
        
        
        for cle in dico_test.keys() :
            vieux_voisin = -1
            voisin = -1
            print(cle)
            print(dico_test[cle][0])
            for elt in dico_test[cle][1] :
                for i in range(elt[0], elt[1]+1) :
                    vieux_voisin = voisin
                    voisin = nb_canoniques_short_range(i, graphes[(cle.split("_")[1], cle.split("_")[2])])
                    if ((vieux_voisin != -1 and abs(vieux_voisin - voisin) < 10) or (vieux_voisin == -1)) and graphes[(cle.split("_")[1], cle.split("_")[2])].nodes[voisin]["part"] in ["Hairpin", "Stem", "Interior Loop"]:
                        print("ok")
                    else :
                        print("a verif")
                        print(vieux_voisin)
                        print(voisin)
                        print(graphes[(cle.split("_")[1], cle.split("_")[2])].nodes[voisin]["part"])
                    
        print(tab_test_2)
        print(len(tab_test_2))
        
        print(compteur)     
        print(tab_val_can)
        print(tab_val_non_can)
        
        print(dico_val_can)
        print(dico_val_non_can)
        
        print(tab_type_non_can)
        
        
        #print(graphes[('5J7L', 'DA')][2543])
#         print(graphes[('4V88', 'A5')][650])
#         print(graphes[('4V88', 'A5')][651])
#         print(graphes[('4V88', 'A5')][652])
#         print(graphes[('4V88', 'A5')][653])
#         print(graphes[('4V88', 'A5')].nodes[649])
#         print(graphes[('4V88', 'A5')].nodes[650])
#         print(graphes[('4V88', 'A5')].nodes[651])
#         print(graphes[('4V88', 'A5')].nodes[652])
#         print(graphes[('4V88', 'A5')].nodes[653])
#         print(graphes[('1GID', 'B')][89])
#         print(graphes[('1GID', 'B')][90])
#         print(graphes[('1GID', 'B')][91])
#         print(graphes[('5FDU', '1A')][2087])
#         print(graphes[('5FDU', '1A')][2088])
#         print(graphes[('5FDU', '1A')][2089])
        #print(graphes[('1FJG', 'A')].nodes[333])
#         print(graphes[('5DM6', 'X')][521])
        print(graphes[('3UCZ', 'R')][54])
        print(graphes[('3UCZ', 'R')][56])
        print(graphes[('4V9F', '0')][627])
        print(graphes[('4V9F', '0')][628])
        print(graphes[('4V9F', '0')][629])
        
        print(graphes[('5FDU', '1A')].nodes[718])
        print(graphes[('5FDU', '1A')].nodes[719])
        print(graphes[('5FDU', '1A')].nodes[720])
#         print(graphes[('1FJG', 'A')][334])
#         print(graphes[('1FJG', 'A')][335])
#         print(graphes[('1FJG', 'A')][336])
#         print(graphes[('1FJG', 'A')][337])
#         print(graphes[('1FJG', 'A')][338])
#         print(graphes[('1FJG', 'A')].nodes[336])
#         print(graphes[('1FJG', 'A')].nodes[1190])
#         print(graphes[('1FJG', 'A')].nodes[1191])
#         print(graphes[('1FJG', 'A')].nodes[1192])
        
        nb = 0
        tab_tot = []
        for elt in tab_val_can :
            for e in elt :
                tab_tot.append(e)
            if len(elt) > 0 :
                nb += 1
        
        print(nb)    
        print(tab_tot)
        
        with open(EXTENSION_TOUTES_ARETES+"fichier_csv_distrib_env_liaisons_non_can.csv", 'w', newline='') as fichier :
            csvwriter = csv.writer(fichier)
            csvwriter.writerow(tab_tot)
            
            
        tab_tot_non_can = []
        for elt in tab_val_non_can :
            for e in elt :
                tab_tot_non_can.append(e)
            
        print(tab_tot_non_can)
        
        intervalles_somme = [0,0,0,0,0,0,0,0,0,0,0]
        intervalles_compteur = [0,0,0,0,0,0,0,0,0,0,0]

        for i in range(len(tab_tot)) :
            if tab_tot[i] < 5 :
                intervalles_somme[tab_tot[i]] += tab_tot_non_can[i]
                intervalles_compteur[tab_tot[i]] += 1
            else :
                intervalles_somme[5] += tab_tot_non_can[i]
                intervalles_compteur[5] += 1
        
        print(intervalles_compteur)
        print(intervalles_somme)
        
        intervalles_tot = [0,0,0,0,0,0,0,0,0,0,0]
        
        for i in range(len(intervalles_somme)) :
            if intervalles_compteur[i] > 0 :
                intervalles_tot[i] = intervalles_somme[i]/intervalles_compteur[i]
            
        print(intervalles_tot)
        
        compteur = 0
        for elt in tab_val_can :
            if len(elt) > 0 :
                compteur += 1
        print(compteur)
        
        with open(EXTENSION_TOUTES_ARETES+"fichier_csv_distrib_env_liaisons_non_can_non_can.csv", 'w', newline='') as fichier :
            csvwriter = csv.writer(fichier)
            csvwriter.writerow(tab_tot_non_can)
        
        dico_type_non_can = {}
        for groupe in tab_type_non_can :
            for elt in groupe :
                for e in elt : 
                    if e not in dico_type_non_can.keys() :
                        dico_type_non_can.update({e : 1})
                    else :
                        dico_type_non_can[e] += 1
        print(dico_type_non_can.keys())
        
        with open(EXTENSION_TOUTES_ARETES+"fichier_csv_distrib_env_type_liaisons_non_can.csv", 'w', newline='') as fichier_type_non_can :
            csvwriter = csv.writer(fichier_type_non_can)
            for cle in dico_type_non_can.keys() :
                csvwriter.writerow([cle, dico_type_non_can[cle]])    
            
        differents = ['fichier_5J7L_DA_25_12_2.pickle', 'fichier_1FJG_A_271_1_2.pickle', 'fichier_5DM6_X_48_28_2.pickle', 'fichier_5FDU_1A_25_68_2.pickle', 'fichier_5J7L_DA_134_1_2.pickle', 'fichier_4V88_A5_25_30_2.pickle', 'fichier_5FDU_1A_290_2_2.pickle', 'fichier_5J7L_DA_272_2_2.pickle', 'fichier_5DM6_X_25_15_2.pickle', 'fichier_4V9F_0_48_16_2.pickle', 'fichier_4V88_A5_287_1_2.pickle', 'fichier_5DM6_X_134_2_2.pickle', 'fichier_5FDU_1A_272_1_2.pickle', 'fichier_4V88_A5_48_6_2.pickle', 'fichier_4V9F_0_48_2_2.pickle', 'fichier_5J7L_DA_48_15_2.pickle', 'fichier_5FDU_1A_48_25_2.pickle', 'fichier_4V9F_0_30_23_2.pickle', 'fichier_4V9F_0_127_6_2.pickle', 'fichier_4V88_A5_30_5_2.pickle', 'fichier_5FDU_1A_134_3_2.pickle', 'fichier_4V9F_0_25_56_2.pickle', 'fichier_4V9F_0_48_13_2.pickle', 'fichier_4V9F_0_207_3_2.pickle', 'fichier_3UCZ_R_62_15_2.pickle', 'fichier_5J7L_DA_48_20_2.pickle', 'fichier_4L81_A_25_77_2.pickle']
     
        for elt in differents :
            if len(dico_val_can[elt]) == 0 :
                print(elt)
        print("ramou")
        for elt in dico_val_can.keys() :
            if len(dico_val_can[elt]) > 0 and elt not in differents :
                print(elt)
        
        compteur = 0
        for elt in tab_type_nt :
            if elt != False :
                print(tab_type_non_can[compteur])
                print(elt)      
            compteur += 1
#         print(len(tab_type_nt))
#         print(tab_type_nt)
        
                
        

        