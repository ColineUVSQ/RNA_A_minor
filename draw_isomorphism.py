'''
Created on 6 déc. 2018

@author: coline
Dessiner les superpositions pour les extensions version matplotlib
'''
from recup_data.constantes import EXTENSION_PATH, EXTENSION_PATH_TAILLE,\
    NEW_EXTENSION_PATH_TAILLE, PATH_MMCIF

import pickle
import networkx as nx
import matplotlib.pyplot as plt
from networkx.classes.function import get_node_attributes
import os
import numpy as np
from collections import OrderedDict
#from recup_data.calcul_sim import calcul_sim_aretes_avec_coeff

liste = ['5J7L_DA_191_3', '5J7L_DA_191_4', '5FDU_1A_301_1', '5J7L_DA_301_2', '5DM6_X_334_1', '5FDU_1A_334_2', '4V9F_0_335_1', '5J7L_DA_335_2', '3JCS_1_137_4', '4V88_A5_290_1', '4V88_A6_314_2', '5J7L_DA_218_3', '4V9F_0_251_2', '1FJG_A_62_8', '5J7L_DA_137_1', '4V9F_0_118_1', '4V9F_0_62_2', '5J7L_DA_271_2', '4V9F_0_224_1', '5DM6_X_197_1', '3GX5_A_138_6', '1FJG_A_317_2', '5J5B_BA_317_1', '1FJG_A_326_1', '5DM6_X_137_3', '5J5B_BA_314_1', '4V9F_0_134_6', '4V9F_0_328_1', '4V9F_0_197_2', '4V9F_0_62_16', '5J7L_DA_282_2', '4V88_A5_137_2', '5FDU_1A_224_3', '5J7L_DA_326_2']
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

#liste_a_tester = [(('5dm6', 9), ('3t1y', 10)),(('4ybb', 54), ('3t1y', 10)),(('4ybb', 12), ('1u9s', 1)),(('4y4o', 12), ('4u3u', 7)),(('4y4o', 28), ('3t1y', 10)),(('2zjr', 3), ('1u9s', 1)),(('1vq8', 3), ('4u3u', 7)),(('1vq8', 3), ('1u9s', 1)),(('1vq8', 15), ('4u3u', 7)),(('3cc2', 6), ('4u3u', 7)),(('3cc2', 21), ('4u3u', 7)),(('3cc2', 21), ('1u9s', 1)),(('3cc2', 22), ('4y4o', 38)),(('3cc2', 22), ('4ybb', 21)),(('5afi', 17), ('3t1y', 10)),(('4v67', 7), ('1u9s', 1)),(('3ccq', 3), ('4ybb', 7)),(('3ccq', 3), ('5ngm', 4)),(('3ccq', 3), ('1u9s', 1)),(('3cpw', 3), ('4u3u', 7)),(('3cpw', 3), ('1u9s', 1)),(('6eri', 11), ('3t1y', 10)),(('6eri', 8), ('4u3u', 7)),(('5wfs', 11), ('4u3u', 7)),(('5wfs', 19), ('1u9s', 1)),(('4u4r', 4), ('3t1y', 10)),(('4u4r', 16), ('4u3u', 7)),(('6ek0', 10), ('1u9s', 1)),(('4ybb', 7), ('4woi', 62)),(('4u27', 54), ('1u9s', 1)),(('4y4o', 38), ('1u9s', 1)),(('4ybb', 21), ('1u9s', 1)),(('6eri', 17), ('4w2f', 40)),(('6eri', 17), ('5mdv', 30)),(('6eri', 17), ('5e81', 27)),(('6eri', 17), ('5f8k', 15)),(('6eri', 17), ('5wdt', 15)),(('6eri', 17), ('5nwy', 12)),(('6eri', 14), ('1u9s', 1)),(('6ek0', 7), ('1u9s', 1)),(('6ek0', 7), ('5nwy', 12)),(('4u3u', 7), ('1u9s', 1))]
liste_a_tester = []
liste_a_tester.extend([(('6eri', 14), ('6eri', 15)),(('5ngm', 1), ('6eri', 14)),(('4u27', 54), ('6eri', 15)),(('4u27', 54), ('5ngm', 1)),(('2zjr', 3), ('3ccq', 3)),(('6eri', 17), ('6ek0', 7)),(('2zjr', 3), ('3cc2', 21)),(('2zjr', 3), ('1vq8', 3)),(('2zjr', 3), ('3cpw', 3)),(('4ybb', 12), ('3cc2', 21)),(('4ybb', 12), ('1vq8', 3)),(('4w2g', 52), ('3cc2', 22)),(('4ybb', 12), ('3cpw', 3)),(('3cc2', 21), ('4v67', 7)),(('1vq8', 3), ('4v67', 7)),(('4v67', 7), ('3cpw', 3)),(('3cc2', 21), ('5wfs', 19)),(('1vq8', 3), ('5wfs', 19)),(('3cpw', 3), ('5wfs', 19)),(('4y4o', 12), ('6ek0', 10)),(('6eri', 8), ('6ek0', 10)),(('4u4r', 16), ('6ek0', 10)),(('4y4o', 12), ('2zjr', 3)),(('2zjr', 3), ('6eri', 8)),(('5wfs', 11), ('6ek0', 10)),(('2zjr', 3), ('5wfs', 11)),(('2zjr', 3), ('4u4r', 16)),(('4y4o', 12), ('4u4r', 18)),(('6eri', 8), ('4u4r', 18)),(('4u4r', 16), ('4u4r', 18)),(('5wfs', 11), ('4u4r', 18))])

liste_a_tester_new = [(('5dm6', 9), ('3t1y', 10)), (('5dm6', 2), ('4y4o', 58)), (('4ybb', 54), ('3t1y', 10)), (('4ybb', 12), ('2qex', 16)), (('4ybb', 12), ('1u9s', 1)), (('4w2g', 52), ('3cc2', 22)), (('4w2g', 52), ('1u9s', 1)), (('4y4o', 28), ('3t1y', 10)), (('2zjr', 1), ('4u4r', 18)), (('2zjr', 3), ('2qex', 16)), (('2zjr', 3), ('1u9s', 1)), (('3cc2', 22), ('4y4o', 38)), (('3cc2', 22), ('4ybb', 21)), (('5afi', 17), ('3t1y', 10)), (('4v7l', 44), ('4u4r', 18)), (('4v67', 7), ('2qex', 16)), (('4v67', 7), ('1u9s', 1)), (('6hma', 3), ('4y4o', 54)), (('3ccq', 3), ('4ybb', 7)), (('3ccq', 3), ('4y4o', 51)), (('3ccq', 3), ('5ngm', 4)), (('6eri', 11), ('3t1y', 10)), (('6eri', 1), ('4ybb', 4)), (('6eri', 1), ('4y4o', 54)), (('6eri', 1), ('6eri', 7)), (('5wfs', 19), ('2qex', 16)), (('5wfs', 19), ('1u9s', 1)), (('6qul', 16), ('4ybb', 4)), (('6qul', 16), ('4y4o', 54)), (('6qul', 16), ('6eri', 7)), (('2qex', 16), ('6ek0', 10)), (('2qex', 16), ('1u9s', 1)), (('4u4r', 4), ('3t1y', 10)), (('6ek0', 10), ('1u9s', 1)), (('4u27', 32), ('4u27', 54)), (('4u27', 54), ('4ybb', 9)), (('4u27', 54), ('5j7l', 10)), (('4u27', 54), ('3t1y', 8)), (('4u27', 54), ('5ngm', 1)), (('4u27', 54), ('6eri', 15)), (('4u27', 54), ('1u9s', 1)), (('4y4o', 38), ('1u9s', 1)), (('4ybb', 21), ('1u9s', 1)), (('6eri', 14), ('1u9s', 1))]
def draw_isomorphisme(taille_ext, grands_graphes):
    #with open("fichier_affichage_isomorphisme_nouvelle_metrique_coeff_all1_cww_non_can_taille_%s.txt"%taille_ext, 'w') as fichier :
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
        with open("dico_graphes_liaisons_near_avec_ou_sans.pickle", 'rb') as fichier_entree :
                mon_depickler_graphes = pickle.Unpickler(fichier_entree)
                dico_fichier_graphes = mon_depickler_graphes.load()
        
        
#         with open(EXTENSION_PATH%taille_ext+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%s.pickle"%taille_ext, 'rb') as fichier_graphe :
#             mon_depickler = pickle.Unpickler(fichier_graphe)
#             dico_graphe = mon_depickler.load()
#             print(len(dico_graphe))
#             with open(EXTENSION_PATH%taille_ext+"sim_extensions_toutes_aretes_coeff_all1_taille_%s.pickle"%taille_ext, 'rb') as fichier_sim : 
#                 mon_depickler_sim = pickle.Unpickler(fichier_sim)
#                 dico_sim = mon_depickler_sim.load()
#                
#                 dico_sim_sorted = OrderedDict(sorted(dico_sim.items(), key= lambda t: t[1], reverse=True))
#                 print(dico_sim_sorted)
                
        #for fic in os.listdir("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext) : 
            #if "graphe_comp"  in fic : 
                #with open("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext +fic, 'rb') as fichier_graphe :
#                 with open(EXTENSION_PATH_TAILLE%taille_ext +"graphe_comp_test_couples_possibles_fichier_5FDU_1A_134_3_fichier_5DM6_X_134_2.pickle", 'rb') as fichier_graphe :
#                     mon_depickler = pickle.Unpickler(fichier_graphe)
#                     dico_graphe = mon_depickler.load()  
                #with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe :    
                #with open("/media/coline/Maxtor/dico_new_120320_sim_par_branche_0.65_avec_liaison_near_plus_infos.pickle", 'rb') as fichier_graphe :    
               # with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_graphe : 
                #with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_graphe : 
                #with open("/media/coline/Maxtor/dico_new_060420_avec_graphes_en_plus.pickle", 'rb') as fichier_graphe : 

                    
                #with open("test_isomorphisme.pickle", 'rb') as fichier_graphe :    
                #with open(NEW_EXTENSION_PATH_TAILLE+"Graphes_communs/Graphes_communs_groupe_sim_par_branche_0.65/dico_graphes.pickle", "rb") as fichier_graphe :

                    
                #with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_612.pickle", 'rb') as fichier_graphe :    
                #with open("dico_algo_heuristique_grands_graphes_taille_4.pickle", 'rb') as fichier_graphe :
                with open("dico_new_060420_avec_graphes_en_plus.pickle", 'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        dico_graphe = mon_depickler.load() 
                        #print(dico_graphe)
                        print(len(dico_graphe))
                        compteur = 0
                    #for couple in liste_a_tester_new :
                        commun = False
                        for elt in dico_graphe.keys() :
                            #print(elt)
                            #if elt == ('4v67', 7) :
                            if elt[0] in [('4ybb', 12), ('3ccq', 3)] and elt[1] in [('4ybb', 12), ('3ccq', 3)]  :
                            #if dico_graphe[elt]["sim"] > 0.7 and dico_graphe[elt]["sim"] < 0.8 :
                            #if compteur < 10 :
                                    if not commun :
                                    
                                        nom1 = "fichier_%s_%s_taille_4.pdb"%(elt[0][0], elt[0][1])
                                        nom2 = "fichier_%s_%s_taille_4.pdb"%(elt[1][0], elt[1][1])
                                        if (nom1, nom2) in rmsd.keys() :
                                            rmsd_couple = rmsd[(nom1, nom2)]
                                        else :
                                            rmsd_couple = rmsd[(nom2, nom1)]
                                        
                                        if not grands_graphes :
                                    
                                            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0][0]+"_"+str(elt[0][1])+"_4.pickle", 'rb') as fichier_graphe1 :
                                                    mon_depickler1 = pickle.Unpickler(fichier_graphe1)
                                                    graphe1 = mon_depickler1.load()
                                                     
                                            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[1][0]+"_"+str(elt[1][1])+"_4.pickle", 'rb') as fichier_graphe2 :
                                                    mon_depickler2 = pickle.Unpickler(fichier_graphe2)
                                                    graphe2 = mon_depickler2.load()
                                        else :
                                            with open("grands_graphes_new_data_taille_4_renommes.pickle", 'rb') as fichier :
                                                mon_depickler = pickle.Unpickler(fichier)
                                                dico_graphe_global = mon_depickler.load()
                                                
                                                graphe1 = dico_graphe_global[elt[0]]
                                                graphe2 = dico_graphe_global[elt[1]]
    #                                     print(type(dico_fichier_graphes[elt[0]]))
    #                                     print(dico_fichier_graphes[elt[0]])
    #                                     print(dico_fichier_graphes[elt[1]])
    #                                     print(dico_graphe[elt]["nums_graphe"])
    #                                     graphe1 = dico_fichier_graphes[elt[0]][dico_graphe[elt]["nums_graphe"][0]]
    #                                     graphe2 = dico_fichier_graphes[elt[1]][dico_graphe[elt]["nums_graphe"][1] - len(dico_fichier_graphes[elt[1]])*dico_graphe[elt]["nums_graphe"][0]]
    #                                     
                                        print(dico_graphe[elt]["graphe"].nodes.data())
                                        print(dico_graphe[elt]["graphe"].edges.data())
                                        #print(graphe1.nodes[27])
                                        #print(graphe2.nodes[27])
    #                                     print(graphe1.edges.data())
    #                                     print(graphe2.edges.data())
    #                                     print(calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[elt]["graphe"], "petit rat", 1, 1, 1))
    #                                     print(dico_graphe[elt]["graphe"].edges.data())
                                        #elt = ('fichier_5J7L_DA_30_15', 'fichier_5FDU_1A_30_17')
                                        #print(dico_graphe[elt].edges.data())
                                    else :
                                        rmsd_couple = -1.0
                                    enlever = False 
                #                         for l in liste : 
                #                                 if l in elt[0] or l in elt[1] :
                #                                     enlever = True
                #                         print(enlever)
                                            
                                    if enlever == False :
                #                                 cle_sim_1 = elt[0].split('_')[1] + "_" + elt[0].split('_')[2] + "_" + elt[0].split('_')[3] + "_" + elt[0].split('_')[4]
                #                                 cle_sim_2 = elt[1].split('_')[1] + "_" + elt[1].split('_')[2] + "_" + elt[1].split('_')[3] + "_" + elt[1].split('_')[4]
                #     
                #                             if (cle_sim_1, cle_sim_2) in dico_sim.keys() :
                #                                     sim = round(dico_sim[(cle_sim_1, cle_sim_2)],2)
                #                             else :
                #                                     sim = round(dico_sim[(cle_sim_2, cle_sim_1)],2)
                                        print(elt)
                                        #sim = round(calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[elt]["graphe"], "petit rat", 1, 1, 1),2)
                                        #print(sim)
                                        sim = dico_graphe[elt]["sim"]
                                        if str(sim)+"__" + str(elt[0]) + "_" + str(elt[1]) + ".png" not in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) :
                                            GC = dico_graphe[elt]["graphe"].copy()
                                            #GC = dico_graphe[elt]["graphe_en_plus"][1].copy()
                    #                         print(GC.nodes.data())
                    #                         
                                            compteur_elt = 0
                    #                         
                                            fig, axs=plt.subplots(figsize=(10, 12), nrows=1, ncols=2)
                    #                         #fig=plt.figure()
                    #             
                                            columns = 2
                                            rows = 1
                    #                         
                                            nom = ""
                                            #elt = ('fichier_5FDU_1A_197_3', 'fichier_5DM6_X_48_9')
                                            print(elt)
                                            
                                            
            #                                 with open("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext+elt[0]+".pickle", "rb") as fichier_1 :
            #                                     mon_depickler_1 = pickle.Unpickler(fichier_1)
            #                                     graphe1 = mon_depickler_1.load()      
            #                                     with open("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext+elt[1]+".pickle", "rb") as fichier_2 :
            #                                         mon_depickler_2 = pickle.Unpickler(fichier_2)
            #                                         graphe2 = mon_depickler_2.load() 
            #                 
            #                                         sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, dico_graphe[elt], elt, 1, 1, 1)
                                            
                                            ''' on veut afficher la similarite et le rang '''                           
            #                                 if (cle_sim_1, cle_sim_2) in dico_sim.keys() :
            #                                     etape = 1
            #                                     rang = 0
            #                                     for elt_sim in dico_sim_sorted :
            #                                         if elt_sim == (cle_sim_1, cle_sim_2) : 
            #                                             rang = etape
            #                                         etape += 1
            #                                      
            #                                 else :
            #                                     etape = 1
            #                                     rang = 0
            #                                     for elt_sim in dico_sim_sorted :
            #                                         if elt_sim == (cle_sim_2, cle_sim_1) : 
            #                                             rang = etape
            #                                         etape += 1
            
                                            
                                            #fig.suptitle('sim = %.2f, rang = %d/4005'%(round(sim,2), rang), fontsize=16)
                                            fig.suptitle('sim = %.2f, rmsd = %s'%(round(sim,2), round(rmsd_couple, 2) if rmsd_couple != None else None), fontsize=16)    
                                            
                                            if commun :
                                                elt = [elt]
                                                with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+elt[0][0]+"_"+str(elt[0][1])+"_4.pickle", 'rb') as fichier_graphe1 :
                                                        mon_depickler1 = pickle.Unpickler(fichier_graphe1)
                                                        graphe1 = mon_depickler1.load()
                                                        
                                            for element in elt :
                                                    print(element)
                                                    print(compteur)
                                                    '''on cree un subplot par graphe de la comparaison'''
                                                    fig.add_subplot(rows, columns, compteur_elt%2+1)
                                                    #fig.axis('off')
                                                    #axs[compteur%2].xaxis.set_visible(False)
                                                    axs[compteur_elt%2].axis("off")
                                                    axs[compteur_elt%2].set_title(element)
                                                
    #                                             with open(NEW_EXTENSION_PATH_TAILLE+"fichier_"+element+".pickle", 'rb') as fichier_entree :
                                                    print(element)
                                                    if compteur_elt == 0 :  
                                                        G = graphe1
                                                    else :
                                                        G = graphe2
    #                                                 mon_depickler = pickle.Unpickler(fichier_entree)
    #                                                 G = mon_depickler.load()
                                                    print(G.nodes.data())
                                                    
                                                    nx.set_node_attributes(G, (33,33), "coordonnees")
                                                    G.nodes[1]["coordonnees"] = (0.0,0.5)
                                                    G.nodes[2]["coordonnees"] = (2.0,0.5)
                                                    G.nodes[3]["coordonnees"] = (0.0,0.0)
                                                    G.nodes[4]["coordonnees"] = (2.0,0.0)
                                                    G.nodes[5]["coordonnees"] = (3.0,0.5)
                                                    
            #                                         fichier.write(str(element)+"\n") 
            #                                         fichier.write(str(G.number_of_nodes())+"\n") 
                                                    #print(G.nodes())
                                                    
                                                    nodes_list = [u for u,d in G.nodes(data=True) if d["type"] != -1] 
                                                    print(nodes_list)
                                                    ordre_noeuds = [1,2,3,4,5]
                                                    
                                                    '''on ordonne les noeuds par chaine'''
                                                    chaines = [[1]]
                                                    for i in range(1,5) :
                                                        compteur = i
                                                        if i != 1 : chaines.append([i])
                                                        liaison_B53 = True
                                                        while liaison_B53 :
                                                            liaison_B53 = False
                                                            temp = compteur
                                    
                                                            for voisin in G.successors(compteur) :
                                                                for arc in G[compteur][voisin] :
                                                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and G[compteur][voisin][arc]["label"] == 'B53' :
                                                                        liaison_B53 = True
                                                                        temp = voisin
                                                                        chaines[len(chaines)-1].append(voisin)
                                                                         
                                                            for voisin in G.predecessors(compteur) :
                                                                for arc in G[voisin][compteur] :
                                                                    if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and G[voisin][compteur][arc]["label"] == 'B53' :
                                                                        liaison_B53 = True
                                                                        temp = voisin
                                                                        chaines[len(chaines)-1].append(voisin)
                                                            compteur = temp
                                                    
                                                    for i in range(4) :
                                                        for elt in chaines[i] :
                                                            if elt not in ordre_noeuds :
                                                                ordre_noeuds.append(elt)
                                                    
                                                    for noeud in ordre_noeuds :
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
            #                                                                 fichier.write("probleme\n")
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
            #                                                                 fichier.write("probleme\n")
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
            #                                                                 fichier.write("probleme\n")
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
            #                                                                 fichier.write("probleme\n")
                                                                            print("probleme") 
                                                                        
                                                    
            #                                         for noeud in G.nodes() :
            #                                             fichier.write(str(noeud) + " " + str(G.nodes[noeud])+"\n")
            #                                         for u,v,edata in G.edges(data=True) :
            #                                             fichier.write(str(u)+ "'" +str(v) + " " + str(edata["label"])+"\n")
                                                    
                                                    #fichier.write(str(G.nodes.data())+"\n")
                                                    #fichier.write(str(G.edges.data())+"\n")                
                                                    
                                                    #print(pos)
                                                    
                                                    #plt.figure(figsize =(5,12))
                                                               
                                                    red_edges = [(1,2),(2,1),(1,3),(3,1),(1,5),(5,1),(2,4),(4,2),(3,4),(4,3),(2,5),(5,2)]
                                                    green_edges = []
                                                    blue_edges = []
                                                    black_edges = []
                                                    grey_edges = []
                                                    weights = []
                                                    #black_edges = [edge for edge in G.edges() if edge not in red_edges]
                                                    edges_list = [(u,v) for u,v,data in G.edges(data=True) if data["label"] != '0']#if data["long_range"] != None]
                                                    #edges_list = [(u,v) for u,v,data in G.edges(data=True) if data["label"] != '0' and contient_arete((u,v), GC, compteur%2)]
                                                    
                                                    print("rapoulou")
                                                    for u,v,edata, in G.edges(data=True) :
                                                        if (u,v) in edges_list :
                                                            if (u,v) not in red_edges :
                                    #                             if contient_arete((u,v), GC, compteur%2) :
                                    #                                 orange_edges.append((u,v))
                                                                if edata["near"] == True :
                                                                    grey_edges.append((u,v))
                                                                elif edata["label"] == "B53" :
                                                                    green_edges.append((u,v))
                                                                elif edata["label"] == "CWW" :
                                                                    blue_edges.append((u,v)) 
                                                                else :
                                                                    black_edges.append((u,v))
                                                            
                                                            ''' on mettra en plus gras les aretes en commun '''
                                                            if contient_arete((u,v), GC, compteur_elt%2) :
                                                                weights.append(3)
                                                                
                                                                print(u,v)
                                                            else :
                                                                weights.append(1)
                                                
                                                    edge_labels=dict([((u,v,),(d["label"]))for u,v,d in G.edges(data=True) if (u,v) in edges_list and d["label"] != 'B53' and d["label"] != 'CWW'])
            
                                                    #print(edge_labels)
                                                   # node_labels=dict([(u,(d["nt"], d["type"]))for u,d in G.nodes(data=True)])## if d["type"] != None])
                                                                                               #print(node_labels)
                                                    
                                                    nodes_list = [u for u,d in G.nodes(data=True) if d["type"] != -1]
                                                    pos = get_node_attributes(G.subgraph(nodes_list), 'coordonnees')
                                                    
                                                    
                                                    '''on mettra en orange les sommets en commun'''
                                                    orange_nodes = []
                                                    pink_nodes = []
                                                    for noeud in nodes_list :
                                                        if contient_sommet(noeud, GC, compteur_elt%2) :
                                                            orange_nodes.append(noeud)
                                                        else :
                                                            pink_nodes.append(noeud)
                                                    node_colors = ['pink' if node in pink_nodes else 'orange' for node in nodes_list]
                                                    
                                                    node_labels=dict([(u, (u,d["type"], d["poids"]))for u,d in G.nodes(data=True) if d["type"] != -1 and d["type"] != None ])#and u in orange_nodes])
        
                                                    #print(nodes_list)
                                                    nx.draw_networkx_nodes(G, pos, nodelist=nodes_list, node_size=300, node_color=node_colors)
                                                                #nx.draw_networkx_nodes(G, pos, cmap=plt.get_cmap('jet'), 
                                                                #           node_color = values, node_size = 500)
                                                    nx.draw_networkx_labels(G, pos, labels = node_labels)
                                                    #nx.draw_networkx_edges(G, pos, edgelist=red_edges, edge_color='r', edge_labels = edge_labels)
                                                    #nx.draw_networkx_edges(G, pos, edgelist=blue_edges, edge_color='b', edge_labels = edge_labels)
                                                    #nx.draw_networkx_edges(G, pos, edgelist=green_edges, edge_color='g', edge_labels = edge_labels)
                                                    #nx.draw_networkx_edges(G, pos, edgelist=black_edges, edge_labels = edge_labels)
                                                    
                                                    
                                                    edge_colors = ['grey' if edge in grey_edges else 'black' if edge in black_edges else 'red' if edge in red_edges else 'blue' if edge in blue_edges else 'green' for edge in edges_list]
                                                    
                                                    
                                                    #nx.draw_networkx_edge_labels(G,pos)
                                                    print(len(edges_list))
                                                    print(len(weights))
                                                    print(weights)
                                                    print(G.edges.data())
                                                    #print(len(edge_labels))
                                                    #nx.draw_networkx_edge_labels(G,pos, edge_labels = edge_labels, font_size=6)
                                                    nx.draw_networkx_edges(G,pos, edgelist=edges_list, edge_color= edge_colors, width=weights)
                                                    #axs[compteur%2].set_axis_off()
                                                    #plt.savefig("graphes_extension/"+element[:len(element)-7]+".png") # save as png
                                                    #plt.savefig("graphes_extension/fichier_1FJG_A_48_8.png") # save as png
                                                    nom = nom + "_" +element[0] +"_"+str(element[1])
                                                    compteur_elt = compteur_elt+1
                                                    
                                                    
                                                    plt.axis('off') 
                                            #plt.savefig("Faux_positifs_png_0.65/"+str(sim)+"_"+nom+".png", format="png")
                     
                                            plt.show()
            #                                 if str(sim)+"_"+nom+".png" not in os.listdir(EXTENSION_PATH_TAILLE%taille_ext) :
            #                                     plt.savefig(EXTENSION_PATH_TAILLE%taille_ext+str(sim)+"_"+nom+".png") # save as png
            
                                            #plt.savefig("/home/coline/Bureau/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_autres_tailles_new_version_avec_non_can/graphes_extension_taille_%d/"%taille_ext+(str(round(sim,2))) + "_"+fic[20:len(fic)-7]+".png") # save as png
                                            plt.clf()
                                            plt.close()
                                        #break
                                        #break
                                #compteur += 1
                            

if __name__ == '__main__':
    for i in range(4, 5) :
        draw_isomorphisme(i, False)             
    