'''
Created on 27 août 2019

@author: coline

Traitement des nouvelles donnees : recherche des types d'ARN, des groupes d'identiques, d'homologues de plusieurs façons

'''

import pickle
import os
import csv
from urllib import request
from gemmi import cif
import multiprocessing
from Bio.Align.Applications import ClustalwCommandline
from recup_data.constantes import PATH_MMCIF, NEW_EXTENSION_PATH_TAILLE
from os import path
import networkx as nx
from urllib import request
from Bio.Blast.Applications import NcbiblastnCommandline
from io import StringIO
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.Blast import NCBIWWW
import matplotlib.pyplot as plt
import urllib
from Bio.Emboss.Applications import NeedleCommandline
from Bio import AlignIO
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.cluster import DBSCAN
import numpy as np
import re
from math import sqrt
from bs4 import BeautifulSoup
from networkx.drawing.layout import spring_layout
import copy
from recup_data.comparaison_clustering_rmsd_sim import calcul_sim_jaccard_pour_clustering_recouvrant,\
    calcul_sim_jaccard
from mpl_toolkits.mplot3d import Axes3D
from recup_data.new_algo_comparaison import creation_graphe_complet
from recup_data.clustering_perez import algo_principal
from recup_data.clustering_mcode import clustering_mcode

liste_pbs = [('6az1', 5), ('6ek0', 9), ('4wqf', 14), ('6ek0', 1), ('6ek0', 8), ('6ek0', 9), ('6az1', 4), ('5wdt', 11), ('6gyv', 1), ('4y4o', 18), ('5dm6', 6), ('4ybb', 27), ('5ibb', 34), ('1vqo', 21), ('5afi', 15), ('6h4n', 23), ('4v9d', 41), ('1k8a', 8), ('4u4r', 28), ('4u3u', 5), ('6ek0', 11), ('2nz4', 3), ('4ena', 1), ('3t1y', 2)]


### Fonctions utilitaires ###

'''04/09/19 
renvoie vrai si les deux graphes d'extension passes en parametres sont identiques d'un point de vue extension et faux sinon
identiques ici : 
noeuds : meme type, meme poids, meme chaine
aretes : meme label, memes noeuds
'''
def identiques(extension1, extension2):
    compteur_idem_n = 0
    for noeud1, data1 in extension1.nodes(data=True) :
        for noeud2, data2 in extension2.nodes(data=True) :
            if noeud1 == noeud2 and data1["type"] == data2["type"] and data1["poids"] == data2["poids"] and data1["chaine"] == data2["chaine"] :
                compteur_idem_n += 1
    
    compteur_idem_a = 0
    for u1, v1, data1 in extension1.edges(data=True) :
        for u2, v2, data2 in extension2.edges(data=True) :
            if u1 == u2 and v1 == v2 and data1["label"] == data2["label"] :
                compteur_idem_a += 1
                
                
    if compteur_idem_n == extension1.number_of_nodes() and compteur_idem_n == extension2.number_of_nodes() :# and compteur_idem_a == extension1.number_of_edges() and compteur_idem_a == extension2.number_of_edges() :
        return True
    else :
        return False


'''04/09/19
renvoie vrai si la difference de positions sur la sequence entre le motif de l'extension elt1 et celui de l'extension elt2 est inferieure a un seuil pour les 3 brins
et faux sinon
(valeurs de position = numeros des noeuds dans le graphe de structure, pas les positions stockees dans la PDB)
dans les deux cas, renvoie aussi les differences de position pour les nucleotides du motif 1, 2 et 5 (representants des 3 brins)'''
def positions_similaires(elt1, elt2, seuil):
#     print(elt1.nodes[1]["position"][0])
#     print(elt2.nodes[1]["position"][0])
#     print(elt1.nodes[2]["position"][0])
#     print(elt2.nodes[2]["position"][0])
    if abs(elt1.nodes[1]["position"][0]- elt2.nodes[1]["position"][0]) < seuil and abs(elt1.nodes[2]["position"][0]- elt2.nodes[2]["position"][0]) < seuil and abs(elt1.nodes[5]["position"][0]- elt2.nodes[5]["position"][0]) < seuil :
        #print(abs(abs(elt1.nodes[1]["position"][0]-elt1.nodes[2]["position"][0])-abs(elt2.nodes[1]["position"][0]-elt2.nodes[2]["position"][0]))) 
        #if abs(abs(elt1.nodes[1]["position"][0]-elt1.nodes[2]["position"][0])-abs(elt2.nodes[1]["position"][0]-elt2.nodes[2]["position"][0])) < 50 :
            return True, abs(elt1.nodes[1]["position"][0]- elt2.nodes[1]["position"][0]), abs(elt1.nodes[2]["position"][0]- elt2.nodes[2]["position"][0]), abs(elt1.nodes[5]["position"][0]- elt2.nodes[5]["position"][0])
        #else :
#             print("rapoutou")
#             print(abs(abs(elt1.nodes[1]["position"][0]-elt1.nodes[2]["position"][0])-abs(elt2.nodes[1]["position"][0]-elt2.nodes[2]["position"][0])))
# 
#             return False, abs(elt1.nodes[1]["position"][0]- elt2.nodes[1]["position"][0]), abs(elt1.nodes[2]["position"][0]- elt2.nodes[2]["position"][0])

    else :
        return False, abs(elt1.nodes[1]["position"][0]- elt2.nodes[1]["position"][0]), abs(elt1.nodes[2]["position"][0]- elt2.nodes[2]["position"][0]), abs(elt1.nodes[5]["position"][0]- elt2.nodes[5]["position"][0])

''' 07/01/20
renvoie vrai si la difference de positions sur la sequence entre le motif de l'extension elt1 et celui de l'extension elt2 est inferieure a un seuil pour les 3 brins
et faux sinon
(valeurs de position = positions stockees dans la PDB)
dans les deux cas, renvoie aussi les differences de position pour les nucleotides du motif 1, 2 et 5 (representants des 3 brins)'''
def pos_similaire_fr3d(ext1, ext2, nom_ext1, nom_ext2, seuil):
    with open("Graphs/%s.pickle"%nom_ext1[0], "rb") as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        graphe1 = mon_depickler.load()
        
    with open("Graphs/%s.pickle"%nom_ext2[0], "rb") as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        graphe2 = mon_depickler.load()
    
    p = re.compile('[-]*[0-9]*')
    
    
    
    num_fr3d = graphe1.nodes[(ext1.nodes[1]["num_ch"],ext1.nodes[1]["position"][0])]["fr3d"] 
    
    if not num_fr3d.isdigit():
        m = p.match(num_fr3d)
        if m :
            num = m.group()
        else :
            print("probleme5")
        pos11 = int(num)
    else :
        pos11 = int(num_fr3d)
    
    num_fr3d = graphe2.nodes[(ext2.nodes[1]["num_ch"],ext2.nodes[1]["position"][0])]["fr3d"]
    if not num_fr3d.isdigit() :
        m = p.match(num_fr3d)
        if m :
            num = m.group()
        else :
            print("probleme5")
        pos21 = int(num)
    else :
        pos21 = int(num_fr3d)
    
    num_fr3d = graphe1.nodes[(ext1.nodes[2]["num_ch"],ext1.nodes[2]["position"][0])]["fr3d"]
    if not num_fr3d.isdigit() :
        m = p.match(num_fr3d)
        if m :
            num = m.group()
        else :
            print("probleme5")
        pos12 = int(num)
    else :
        pos12 = int(num_fr3d)
    
    num_fr3d = graphe2.nodes[(ext2.nodes[2]["num_ch"],ext2.nodes[2]["position"][0])]["fr3d"]
    if not num_fr3d.isdigit() :
        m = p.match(num_fr3d)
        if m :
            num = m.group()
        else :
            print("probleme5")
        pos22 = int(num)
    else :
        pos22 = int(num_fr3d)
        
    
    num_fr3d = graphe1.nodes[(ext1.nodes[5]["num_ch"],ext1.nodes[5]["position"][0])]["fr3d"]
    if not num_fr3d.isdigit() :
        m = p.match(num_fr3d)
        if m :
            num = m.group()
        else :
            print("probleme5")
        pos15 = int(num)
    else :
        pos15 = int(num_fr3d)
    
    num_fr3d = graphe2.nodes[(ext2.nodes[5]["num_ch"],ext2.nodes[5]["position"][0])]["fr3d"]
    if not num_fr3d.isdigit() :
        m = p.match(num_fr3d)
        if m :
            num = m.group()
        else :
            print("probleme5")
        pos25 = int(num)
    else :
        pos25 = int(num_fr3d)

    if abs(pos11- pos21) < seuil and abs(pos12- pos22) < seuil and abs(pos15- pos25) < seuil :
        return True, abs(pos11- pos21), abs(pos12- pos22), abs(pos15- pos25)    
    else :
        return False, abs(pos11- pos21), abs(pos12- pos22), abs(pos15- pos25)  

''' 05/09/19
dans une liste de listes, fusionner deux sous-listes d'indices indice1 et indice2 dans la liste notee groupes '''
def fusion(groupes, indice1, indice2):
    
    for elt in groupes[indice2] :
        if elt not in groupes[indice1] :
            groupes[indice1].append(elt)
            
    groupes.pop(indice2)
    return groupes


''' 09/09/19 
recherche le nombre de groupes d'homologues par type d'ARN et le nombre d'occurrences sans homologues
version 1 de groupes d'homologues (positions similaires du motif a moins de 25 nts de difference)
'''
def nombre_de_groupes_homologues(num_ARN):
    with open(NEW_EXTENSION_PATH_TAILLE+"groupe_%s.pickle"%num_ARN, 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        liste = mon_depickler.load()
    
    nb_seuls = 0
    with open("groupes_%s_homologues.pickle"%num_ARN, 'rb') as fichier_homologues :
        mon_depickler_1 = pickle.Unpickler(fichier_homologues)
        groupes_homologues = mon_depickler_1.load()
        
        nb_homologues = 0
        for groupe in groupes_homologues :
            nb_homologues += 1
            
        for elt in liste :
            existe = False
            for groupe in groupes_homologues :
                for elt_homologue in groupe :
                    if elt == elt_homologue :
                        existe = True
            if not existe :
                nb_seuls += 1
                
    print(nb_seuls)
    print(nb_homologues)
    

''' 09/09/19 
recherche le nombre de groupes d'identiques par type d'ARN et le nombre d'occurrences sans autre occurrence identique
version 1 de groupes d'identiques (au sein d'un groupe d'homologues version 1, meme numeros de sommets du motif (dans le graphe de motif ou dans le fichier PDB (champ fr3d)), memes nts, memes aretes)
'''
def nombre_de_groupes_identiques(num_ARN):
    with open(NEW_EXTENSION_PATH_TAILLE+"groupe_%s.pickle"%num_ARN, 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        liste = mon_depickler.load()

    
    nb_seuls = 0
    with open("groupes_%s_identiques.pickle"%num_ARN, 'rb') as fichier_identiques :
        mon_depickler_1 = pickle.Unpickler(fichier_identiques)
        groupes_identiques = mon_depickler_1.load()
        
        nb_identiques = 0
        for groupes in groupes_identiques :
            for groupe in groupes :
                #for elt in groupe :
                nb_identiques += 1
            
        for elt in liste :
            existe = False
            for groupes in groupes_identiques :
                for groupe in groupes :
                    for elt_identique in groupe :
                        if elt == elt_identique :
                            existe = True
            if not existe :
                nb_seuls += 1
                
    print(nb_seuls)
    print(nb_identiques)


''' 18/12/19 
renvoie la concatenation des sequences des 3 brins de l'occurrence dont on passe le (num PDB, num_occ) en parametre '''
def concatenate_sequences(nom_extension):
    with open("Nouvelles_donnees/fichier_%s_%s.pickle"%(nom_extension[0], nom_extension[1]), "rb") as fichier_extension : 
        mon_depickler = pickle.Unpickler(fichier_extension)
        extension = mon_depickler.load()        
        seq1,seq2,seq3 = recup_sequences_autour_motif(nom_extension[0], extension, 30)
        
        return seq1+seq2+seq3
    
''' renvoie l'ensemble des composantes connexes du graphe passe en parametre sous forme d'une liste de noeuds '''
def recherche_composante_connexe(graphe_sortie):
    deja_vu = []
    composantes_connexes = []
    for noeud in graphe_sortie.nodes() :
        if noeud not in deja_vu :
            composantes_connexes.append([noeud])
            deja_vu.append(noeud)
             
            #parcours en largeur
            file_sommets = [noeud]
            while len(file_sommets) > 0:
                sommet_courant = file_sommets.pop(0)
                enfants_courant = graphe_sortie[sommet_courant]  
                #print(len(file_sommets))
                for enfant in enfants_courant:
                    if enfant not in deja_vu :
                        composantes_connexes[len(composantes_connexes)-1].append(enfant)
                        file_sommets.append(enfant)
                    deja_vu.append(enfant)
#                 print("sommet_courant")
#                 print(sommet_courant)
#                 print(file_sommets)
#     print("composantes connexes")
#     print(composantes_connexes)
#     print(len(composantes_connexes))
    return composantes_connexes    
    
'''04/10/19
renvoie la sequence autour d'une occurrence de motif
nom_extension : num PDB de l'occurrence
extension : graphe d'extension de l'occurrence
nb_nts : nombre de nucleotides autour du motif qu'il faut considerer'''
def recup_sequences_autour_motif(nom_extension, extension, nb_nts):    
    with open("Graphs/%s.pickle"%nom_extension, 'rb') as fichier_graphe :
        mon_depickler = pickle.Unpickler(fichier_graphe)
        graphe = mon_depickler.load()
        
        seq_1 = ""
        seq_2 = ""
        seq_3 = "" 
#         compteur = 0
        for i in range(extension.nodes[3]["position"][0]-nb_nts, extension.nodes[1]["position"][0]+nb_nts+1) :
#             if compteur == 50 : 
#                 print(i)
#                 print((extension.nodes[3]["num_ch"], i))
            if (extension.nodes[3]["num_ch"], i) in graphe.nodes() :
                seq_1 += graphe.nodes[(extension.nodes[3]["num_ch"], i)]["nt"]
#                 print(compteur)
#                 compteur += 1
                
        for i in range(extension.nodes[2]["position"][0]-nb_nts, extension.nodes[4]["position"][0]+nb_nts+1) :
            if (extension.nodes[2]["num_ch"], i) in graphe.nodes() :
                seq_2 += graphe.nodes[(extension.nodes[2]["num_ch"], i)]["nt"]
                
        for i in range(extension.nodes[5]["position"][0]-nb_nts, extension.nodes[5]["position"][0]+nb_nts+2) :
            if (extension.nodes[5]["num_ch"], i) in graphe.nodes() :
                seq_3 += graphe.nodes[(extension.nodes[5]["num_ch"], i)]["nt"]
                
        return seq_1, seq_2, seq_3

### Traitement des donnees un peu a l'aveugle  ###

''' Recherche des redondances dans les occurrences des A-minor
(=ceux qui ont les memes nucleotides en meme position dans la meme molecule)
mais peut etre garder les deux numeros de chaines correspondantes car des differences existent au niveau des aretes des graphes entiers
 non utilisee finalement'''
def recup_aminor():
    with open("fichiers_idem_all_aminor.txt", 'w') as fichier_ecriture :
        with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
            mon_depickler = pickle.Unpickler(fichier_all_aminor)
            all_aminor = mon_depickler.load()
            
            #print(all_aminor)
            print(len(all_aminor))
            
            somme_graphes = 0
            compteur_plusieurs_chaines_identiques = 0
            compteur_plusieurs_chaines_identiques_par_aretes_aussi = 0
            compteur_identiques_noeuds_aretes = 0
            for cle in all_aminor.keys() :
                #if cle == '4v9i' :
                #if cle == '5fdu' :
                    fichier_ecriture.write("####")
                    fichier_ecriture.write(cle+'\n')
                    #print(cle)
                    fichier_ecriture.write("Nombre de graphes : ")
                    fichier_ecriture.write(str(len(all_aminor[cle])))
                    fichier_ecriture.write("\n")
                    #print(len(all_aminor[cle]))
                    #print(all_aminor[cle])
                    idem = 0
                    idem_aretes = 0
                    identiques_noeuds_aretes = 0
                    for i in range(len(all_aminor[cle])) :
                        #print(all_aminor[cle][i].nodes.data())
                        graphe1 = all_aminor[cle][i]
                        for j in range(i+1, len(all_aminor[cle])) :
                            
                            graphe2 = all_aminor[cle][j]
                            
                            liste1 = []
                            liste2 = []
                            compteur_idem = 0
                            for noeud1, data1 in graphe1.nodes(data=True) :
                                liste1.append(noeud1[0])
                                for noeud2, data2 in graphe2.nodes(data=True) :
                                    liste2.append(noeud2[0])
                                    if (noeud1[1] == noeud2[1] or data1["fr3d"] == data2["fr3d"])  and data1["nt"] == data2["nt"] and data1["real_nt"] == data2["real_nt"] :
                                        compteur_idem += 1
                            #print(compteur_idem)
                            if compteur_idem  >= graphe1.number_of_nodes() and compteur_idem >= graphe2.number_of_nodes() :
                                idem += 1
                            #print(idem)
                            compteur_idem_aretes = 0    
                            
                            #print(graphe1.edges.data())
                            #print(graphe2.edges.data())
                            for u1,v1, data1 in graphe1.edges(data=True) :
                                for u2,v2,data2 in graphe2.edges(data=True) :
                                    if ((u1[1] == u2[1] and v1[1] == v2[1]) or (graphe1.nodes[u1]["fr3d"] == graphe2.nodes[u2]["fr3d"] and graphe1.nodes[v1]["fr3d"] == graphe2.nodes[v2]["fr3d"])) and graphe1.nodes[u1]["nt"] == graphe2.nodes[u2]["nt"] and graphe1.nodes[u1]["real_nt"] == graphe2.nodes[u2]["real_nt"] and graphe1.nodes[v1]["nt"] == graphe2.nodes[v2]["nt"] and graphe1.nodes[v1]["real_nt"] == graphe2.nodes[v2]["real_nt"] and data1["label"] == data2["label"] : #and data1["near"] == data2["near"] :
                                        compteur_idem_aretes += 1
                            
                            if compteur_idem_aretes == graphe1.number_of_edges() and compteur_idem_aretes == graphe2.number_of_edges() :           
                                idem_aretes += 1
                            
                            if compteur_idem >= graphe1.number_of_nodes() and compteur_idem >= graphe2.number_of_nodes() and (compteur_idem_aretes == graphe1.number_of_edges() and compteur_idem_aretes == graphe2.number_of_edges()):
                                identiques_noeuds_aretes += 1
                                
                                num_ch = []
                                for noeud in graphe2.nodes() :
                                    if noeud[0] not in num_ch :
                                        num_ch.append(noeud[0])
                                        
                                if len(num_ch) >  1 :
                                    print("petit rat")
                                    print(cle)
                                    print(j)
                                    print(num_ch)
                                
#                             if compteur_idem == graphe1.number_of_nodes() and compteur_idem == graphe2.number_of_nodes() and not (compteur_idem_aretes == graphe1.number_of_edges() and compteur_idem_aretes == graphe2.number_of_edges()):
#                                 print(cle)
#                                 print(graphe1.nodes.data())
#                                 print(graphe1.edges.data())
#                                 print(graphe2.nodes.data())
#                                 print(graphe2.edges.data())
                            
                            if not (compteur_idem >= graphe1.number_of_nodes() and compteur_idem >= graphe2.number_of_nodes()) and (compteur_idem_aretes == graphe1.number_of_edges() and compteur_idem_aretes == graphe2.number_of_edges()):
                                print(cle)
                                print(graphe1.nodes.data())
                                print(graphe1.edges.data())
                                print(graphe2.nodes.data())
                                print(graphe2.edges.data())                                
                                #print(liste1)
                                #print(liste2)
                    #print(idem)
                    fichier_ecriture.write("Nombre de graphes identiques : ")            
                    fichier_ecriture.write(str(idem))
                    fichier_ecriture.write("\n")
                    if idem != 0 :
                        compteur_plusieurs_chaines_identiques += idem
                    
                    if idem_aretes != 0 :
                        compteur_plusieurs_chaines_identiques_par_aretes_aussi += idem_aretes   
                    if identiques_noeuds_aretes != 0 :
                        compteur_identiques_noeuds_aretes += identiques_noeuds_aretes
                    #print(idem)
                    somme_graphes += len(all_aminor[cle])
                    #break
                
            print(somme_graphes)
            print(somme_graphes/len(all_aminor))   
            print(compteur_plusieurs_chaines_identiques)
            print(compteur_plusieurs_chaines_identiques_par_aretes_aussi)
            print(compteur_identiques_noeuds_aretes)

'''
Recherche des chaines identiques dans les graphes correspondant a toute la molecule
identiques = possedent la meme sequence
mais en fait, c'est tres long
non utilisee
'''
def test_chaines_identiques():
    idem = 0
    for fic in os.listdir("Graphs") :
        
        if idem == 0 :
            with open("Graphs/"+fic, 'rb') as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                graphe = mon_depickler.load()
                
                #print(graphe.nodes.data())
                
                ## Rangement par chaine 
                liste_graphes_par_chaine = []
                chaine = -1
                liste_noeuds_temp = []
                for noeud, data in graphe.nodes(data=True) :
                    if chaine != noeud[0] :
                        chaine = noeud[0]
                        if len(liste_noeuds_temp) > 0 :
                            liste_graphes_par_chaine.append(liste_noeuds_temp)
                            liste_noeuds_temp = list([])
                    liste_noeuds_temp.append((noeud, data))
                 
                #print(liste_graphes_par_chaine)   
                for elt in liste_graphes_par_chaine :
                    print(elt)
                print(fic)
                
                ## Recherche chaines identiques
                idem = 0
                for i in range(len(liste_graphes_par_chaine)) :
                    for j in range(i+1,len(liste_graphes_par_chaine)) :
                        
                        liste_noeuds1 = liste_graphes_par_chaine[i]
                        liste_noeuds2 = liste_graphes_par_chaine[j]
                        nom_ch1 = liste_noeuds1[0][0][0]
                        nom_ch2 = liste_noeuds2[0][0][0]
                        compteur_idem = 0
                        for noeud1, data1 in liste_noeuds1 :
                            for noeud2, data2 in liste_noeuds2 :
                                if noeud1[1] == noeud2[1] or (data1["fr3d"] == data2["fr3d"] and data1["nt"] == data2["nt"] and data1["real_nt"] == data2["real_nt"]) :
                                    compteur_idem += 1
                        
                        if compteur_idem == len(liste_noeuds1) and compteur_idem == len(liste_noeuds2) :
                            idem += 1
                            
                            ## Recherche des differences au niveau des aretes entre chaines identiques
                            print(liste_noeuds1)
                            for u,v, data in graphe.edges(data=True) : 
                                #print("gros rat")
                                #print(nom_ch1)
                                #print(nom_ch2)
                                if nom_ch1 == u[0] and nom_ch1 == v[0] :

                                    if ((nom_ch2, u[1]), (nom_ch2, v[1])) not in graphe.edges() :
                                        print("rapoulou")
                                        print(nom_ch2)
                                        print(u,v)
                                        print(data)
                                        
                                elif nom_ch1 == u[0] and nom_ch1 != v[0] :

                                    if ((nom_ch2, u[1]), (v[0], v[1])) not in graphe.edges() :
                                        print("rapoulou")
                                        print(u,v)
                                        print(data)
                                elif nom_ch1 != u[0] and nom_ch1 == v[0] :
                                    
                                    if ((u[0], u[1]), (nom_ch2, v[1])) not in graphe.edges() :
                                        print("rapoulou")
                                        print(u,v)
                                        print(data)
                                
                                if nom_ch2 == u[0] and nom_ch2 == v[0] :
                 
                                    if ((nom_ch1, u[1]), (nom_ch1, v[1])) not in graphe.edges() :
                                        print("rapoulou")
                                        print(nom_ch1)
                                        print(u,v)
                                        print(data)
                                elif nom_ch2 == u[0] and nom_ch2 != v[0] :
                                    
                                    if ((nom_ch1, u[1]), (v[0], v[1])) not in graphe.edges() :
                                        print("rapoulou")
                                        print(u,v)
                                        print(data)
                                elif nom_ch2 != u[0] and nom_ch2 == v[0] :
                                    
                                    if ((u[0], u[1]), (nom_ch1, v[1])) not in graphe.edges() :
                                        print("rapoulou")       
                                        print(u,v)
                                        print(data)
                                
                print(idem)  
                print(fic)
                #print(graphe.edges.data())
#                 for u,v, data in graphe.edges(data=True) : 
#                     print(u)
#                     print(v)
#                     print(data)  
                    
                #break          





'''04/09/19
recherche identiques parmi les graphes d'extensions et stocke leurs references dans un fichier
je ne sais pas si on utilise vraiment cette fonction quelque part '''
def identiques_extensions():
    
    idem = 0
    groupes = []
#     liste_a_faire = []
    liste_fichiers = []
    for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
        if "pickle" in fic : #and ("4tud" in fic or "4v8i" in fic):
            liste_fichiers.append(fic)   

#     with open("groupes_28S_identiques.pickle", 'rb') as fichier :
#         mon_depickler_1 = pickle.Unpickler(fichier)
#         groupes_identiques = mon_depickler_1.load()
        
    
#         for elt in groupes_identiques :
#         for liste_fichiers in elt :
#             print(len(liste_fichiers))
#             print(liste_fichiers)
    for i in range(len(liste_fichiers)) :
                with open(NEW_EXTENSION_PATH_TAILLE + liste_fichiers[i], 'rb') as fichier_1 :
                    mon_depickler_1 = pickle.Unpickler(fichier_1)
                    extension1 = mon_depickler_1.load()
                    for j in range(i+1, len(liste_fichiers)) :
                        print(i)
                        print(j)
                        with open(NEW_EXTENSION_PATH_TAILLE +liste_fichiers[j], 'rb') as fichier_2 :
                            mon_depickler_2 = pickle.Unpickler(fichier_2)
                            extension2 = mon_depickler_2.load()
                            #liste_a_faire.append((extension1, extension2))
                            
                            if identiques(extension1, extension2) :
                                idem += 1
                                groupes.append([liste_fichiers[i], liste_fichiers[j]])
                            else :
                                print("gros rat")
#     p = multiprocessing.Pool(1)
#     result = p.map(identiques, liste_a_faire)
    with open("fichiers_idem_taille_4.pickle", 'wb') as fichier :
        mon_pickler = pickle.Pickler(fichier)
        mon_pickler.dump(groupes)                    
    print(idem)
             

### Recuperation d'informations sur les types d'ARN des occurrences ###

'''28/08/19
Recuperation des infos sur chaque chaine de molecule possedant des occurrences de A-minor
stockees dans fichier csv '''
def recup_type_molecule(rep_entree, rep_sortie):
    #with open("fichier_infos_molecules_new.csv", 'w', newline='') as fichier_csv :
    with open(rep_sortie, 'w', newline='') as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(['Libelle PDB', 'Num PDB', 'Num chaine', "Type d'ARN", 'Organisme', 'Nombre de motifs A-minor', 'Taille de la chaine (nb de nts)'])
        #with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
        with open(rep_entree, 'rb') as fichier_all_aminor :
            mon_depickler = pickle.Unpickler(fichier_all_aminor)
            all_aminor = mon_depickler.load()
            nombre_plusieurs_chaines = 0
            for cle in all_aminor.keys() :
                #if cle == '1xmo' :
                    file = PATH_MMCIF+ cle.upper()+".cif"
                    print(cle)
                    if cle.upper()+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
                        print("petit rat")
                        url = 'http://www.rcsb.org/pdb/files/%s.cif' % cle.upper()
                        request.urlretrieve(url, file)
                    
                    
                    dico_num_chaines = {}
                    for graphe in all_aminor[cle] :
                        liste_nums_chaines_par_graphe = []
                        for noeud in graphe.nodes() :
                            if noeud[0] not in liste_nums_chaines_par_graphe :
                                liste_nums_chaines_par_graphe.append(noeud[0])
                                
                        for elt in liste_nums_chaines_par_graphe :
                            if elt not in dico_num_chaines.keys() :
                                dico_num_chaines.update({elt : 1})
                            else :
                                dico_num_chaines[elt] += 1
                                
                    print(dico_num_chaines)
                    print(len(all_aminor[cle]))
                    
    #                         if noeud[0] not in liste_nums_chaines_par_graphe :
    #                             liste_nums_chaines_par_graphe.append(noeud[0])
                                
    #                     if len(liste_nums_chaines_par_graphe) >  1 :
    #                         print("plusieurs chaines")
    #                         
    #                         print(cle)
    #                         print(graphe.nodes.data())
    #                         print(liste_nums_chaines_par_graphe)
    #                         
    #                         nombre_plusieurs_chaines += 1
                            
                #print(nombre_plusieurs_chaines)
                            
                    
                     
                    try : 
                        doc = cif.read_file(PATH_MMCIF+cle.upper()+".cif")
                        block = doc.sole_block()
                        print(block.name)
                     
                        for elt in dico_num_chaines.keys() : 
                         
                            cat = block.find_mmcif_category("_entity_poly.") 
                            print(list(cat.tags))
                            id = -1
                            for row in cat :
                                #print(tuple(row))
                                num_ch = row[6].split(",")
                                
    #                             print("petit rat")
    #                             print(num_ch)
                                if elt in num_ch :
                                    id = row[0]
                             
                            cat = block.find_mmcif_category("_struct_ref_seq.")
                            
                            for row in cat :
                                #print(tuple(row))
                                if row[1] == id :
                                    taille = int(row[6]) - int(row[4]) +1
                            
                             
                            cat = block.find_mmcif_category("_entity.") 
                              
                            for row in cat :
                                #print(tuple(row))
                                if row[0] == id :
                                    type_ARN = row[3]
                                      
                            cat = block.find_mmcif_category("_entity_src_nat.")         
                                     
                            for row in cat :
                                print(tuple(row))
                                if row[0] == id :
                                    org = row[6]
                                     
                            cat = block.find_mmcif_category("_pdbx_entity_src_syn.")         
                                     
                            for row in cat :
                                print(tuple(row))
                                if row[0] == id :
                                    org = row[5]
                                      
                            print(type_ARN)
                            print(org)
                            print(taille)
                            
                            cat = block.find_mmcif_category("_struct.")
                            print(cat) 
                            compteur = 0
                            for tag in cat.tags :
                                if tag == '_struct.title' :
                                    title = cat[0][compteur]
                                    
                                compteur += 1
                            print(title)
                            
                            csvwriter.writerow([title, cle.upper(), elt, type_ARN, org, dico_num_chaines[elt], taille])
                    except RuntimeError:
                        print("probleme de fichier cif : %s"%cle.upper())
                        

''' modifie le 15/10/19 
Recherche les occurrences A-minor reliant plusieurs types differents d'ARN
et les stocke dans un fichier pickle (occ_multi_chaine)
rep_entree : chemin du fichier de stockage des occurrences du motif 
'''
def entre_plusieurs_chaines(rep_entree):
    groupes_arn = {"ARNt", "Riboswitch", "Ribozyme", "Intron", "SRP", "ARNm", '23S', '28S', '16S', '18S', '25S', '4.5S', '4.8S', 'LSU-alpha', 'LSU-beta'}
    with open(rep_entree, 'rb') as fichier_all_aminor :
            mon_depickler = pickle.Unpickler(fichier_all_aminor)
            all_aminor = mon_depickler.load()
            nombre_plusieurs_chaines = 0
            liste_plusieurs_chaines = []
            
            tous_les_types = {}
            les_arnt_16s = []
            les_arnm_arnt_16S = []
            for cle in all_aminor.keys() :
                #if cle == '5t2a' :
                    #print(cle)
                    file = PATH_MMCIF+ cle.upper()+".cif"
                    compteur = 1
                    
                    
                    
                    for graphe in all_aminor[cle] :
                        liste_nums_chaines_par_graphe = []
                        for noeud in graphe.nodes() :
                            if noeud[0] not in liste_nums_chaines_par_graphe :
                                liste_nums_chaines_par_graphe.append(noeud[0])
                                
                        print(liste_nums_chaines_par_graphe)        
                        if len(liste_nums_chaines_par_graphe) >  1 :
                            print("petit rat")
                            nombre_plusieurs_chaines += 1
                            liste_plusieurs_chaines.append((cle, compteur))
                            
                            doc = cif.read_file(PATH_MMCIF+cle.upper()+".cif")
                            block = doc.sole_block()
                            #print(block.name)
                            
                            
                            types_ARN = []
                            
                            for elt in liste_nums_chaines_par_graphe : 
                             
                                cat = block.find_mmcif_category("_entity_poly.") 
                                #print(list(cat.tags))
                                id = -1
                                for row in cat :
                                        #print(tuple(row))
                                    num_ch = row[6].split(",")
                                        
            #                             print("petit rat")
            #                             print(num_ch)
                                    if elt in num_ch :
                                        id = row[0]
                                
                                cat = block.find_mmcif_category("_entity.") 
                                      
                                for row in cat :
                                        #print(tuple(row))
                                    if row[0] == id :
                                        compter = 0
                                        if "ribosomal" in row[3] or "Ribosomal" in row[3] or "RIBOSOMAL" in row[3] or "rRNA" in row[3] or "ribomosomal" in row[3] or "RRNA" in row[3] or "ribosome" in row[3] or "ribsomal" in row[3] or "ribosmal" in row[3] :
                                            liste_arnr = ['23S', '28S', '16S', '18S', '25S', '4.5S', '4.8S']
                                            for gr in liste_arnr :
                                                if gr in row[3] :
                                                    types_ARN.append(gr)
                                                    compter += 1
                                            if compter == 0 :
                                                
                                                print(block.name)
                                                print(row[3])
                                                types_ARN.append("autre_ARNr")
                                                compter +=1
                                            
                                        elif "tRNA" in row[3] or "TRNA" in row[3] or "transfer RNA" in row[3] or "Transfer RNA" in row[3] or "trNA" in row[3] or "TRANSFER RNA" in row[3] :
                                            types_ARN.append('ARNt')
                                            compter += 1
                                        elif "riboswitch" in row[3] or "Riboswitch" in row[3] or "RIBOSWITCH" in row[3]:
                                            types_ARN.append("Riboswitch")
                                            compter += 1
                                        elif "ribozyme" in row[3] or "Ribozyme" in row[3] or "RIBOZYME" in row[3] or "Ribonuclease" in row[3] or "ribonuclease" in row[3] or "RIBONUCLEASE" in row[3] or "RNase" in row[3] :
                                            types_ARN.append("Ribozyme")
                                            compter += 1
                                        elif "mRNA" in row[3] or "MRNA" in row[3] or "messenger RNA" in row[3] or "Messenger RNA" in row[3] or "MESSENGER RNA" in row[3] :
                                            types_ARN.append("ARNm")
                                            compter += 1
                                        elif "intron" in row[3] or "INTRON" in row[3] :
                                            types_ARN.append("Intron")
                                            compter += 1
                                        elif "SRP" in row[3] :
                                            types_ARN.append("SRP")
                                            compter += 1
                                        else :
                                            types_ARN.append("autre")
                                            compter += 1
                                            print(block.name)
                                            print(row[3])
                                        if compter != 1 :
                                            print("bizarre")
                                            print(compter)
                                            print(row[3])
                                            exit(0)
                            
#                             print(types_ARN)
#                             print(tous_les_types)  

                            
                              
                            est_passe = 0       
                            ajoute = False
                            for tous_types in tous_les_types.keys() :
                                compter = 0
                                deja_utilise = []
                                for elt in types_ARN :
                                    if elt in tous_types and tous_types.index(elt) not in deja_utilise :
                                        deja_utilise
                                        compter += 1
                                        deja_utilise.append(tous_types.index(elt))
                                if compter == len(types_ARN) and compter == len(tous_types) :
                                    tous_les_types[tous_types] += 1
                                    #print("petit rat")
                                    est_passe += 1
                                    ajoute = True
                                    
                                    if tous_types == ('ARNt', '16S') :
                                        les_arnt_16s.append((cle, compteur))
                                        
                                    if tous_types == ('ARNm', '16S', 'ARNt') :
                                        les_arnm_arnt_16S.append((cle, compteur))
                                    
                            if not ajoute :
                                tous_les_types.update({tuple(types_ARN) : 1 })
                                #print("gros rat")
                                est_passe += 1
                            if est_passe != 1 :
                                print("bizarre")
                                print(types_ARN)
                                print(tous_les_types)
                                exit(0)
                        compteur += 1

            print(tous_les_types)
            somme = 0
            for val in tous_les_types.values() :       
                somme += val
            print(somme)
    print(nombre_plusieurs_chaines)
    print(liste_plusieurs_chaines)  
        
#     with open("groupe_arnt_16s.pickle", 'wb') as fichier_multi_chaines_arnt_16S :
#         mon_pickler = pickle.Pickler(fichier_multi_chaines_arnt_16S)
#         mon_pickler.dump(les_arnt_16s) 
#     
#     with open("groupe_arnt_16s_arnm.pickle", 'wb') as fichier_multi_chaines_arnt_16S_arnm :
#         mon_pickler = pickle.Pickler(fichier_multi_chaines_arnt_16S_arnm)
#         mon_pickler.dump(les_arnm_arnt_16S) 
#     
#     with open("occ_multi_chaine.pickle", 'wb') as fichier_multi_chaines :
#         mon_pickler = pickle.Pickler(fichier_multi_chaines)
#         mon_pickler.dump(liste_plusieurs_chaines)

''' renvoie la liste des occurrences dont deux chaines se rejoignent par la sequence
(sous-forme de dictionnaire avec les numeros des deux chaines comme cles et une liste d'occurrences comme valeurs) '''
def chercher_liaison_entre_chaines_par_sequence():
    liste_liens_entre_chaines = {(1,2) : 0, (1,4) : 0, (2,3) : 0, (3,4) : 0}
    liste = []
    for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
        if "fichier" in fic and "pickle" in fic :
            with open(NEW_EXTENSION_PATH_TAILLE+fic, 'rb') as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                graphe = mon_depickler.load()
                
                if abs(graphe.nodes[1]["position"][0] - graphe.nodes[2]["position"][0]) < 10 :
                    liste_liens_entre_chaines[(1,2)] += 1
                    liste.append(fic)
                    print(fic)
                    
#                 if abs(graphe.nodes[1]["position"][0] - graphe.nodes[4]["position"][0]) < 10 :
#                     liste_liens_entre_chaines[(1,4)] += 1
                    
                if abs(graphe.nodes[2]["position"][0] - graphe.nodes[3]["position"][0]) < 10 :
                    liste_liens_entre_chaines[(2,3)] += 1
                    liste.append(fic)
                
#                 if abs(graphe.nodes[3]["position"][0] - graphe.nodes[4]["position"][0]) < 10 :
#                     liste_liens_entre_chaines[(3,4)] += 1
                    
    print(liste_liens_entre_chaines)

'''11/09/19
range les occurrences de typ_motif (A-minor ou autres) par type d'ARN, a partir du fichier csv ou sont stockees les metadonnees des fichiers PDB 
et stocke les references des occurrences (num_PDB+num_occ) dans des fichiers pickle (un par type d'ARN)
(mais ne fait pas la difference pour les occurrences inter-chaines : les stocke dans chaque type d'ARN que l'occurrence contient)
num_ARN : le type d'ARN dont on veut recuperer les occurrences
rep_csv : chemin du fichier csv des metadonnees
rep_pickle : chemin du fichier de stockage des occurrences du motif selectionne
typ_motif : le type de motif (A-minor, motif 6left...)
 ''' 
def groupes_ARN(num_ARN, rep_csv, rep_pickle, typ_motif):
    #with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv :
    with open(rep_csv, 'r', newline='') as fichier_csv :
        csvreader = csv.reader(fichier_csv) 
        groupes_ARNr = {'23S' : [], '28S' : [], '16S' : [], '18S' : [], '25S' : [], '4.5S' : [], '4.8S' : [], 'LSU-alpha' : [], 'LSU-beta' : [], '5S' : []}
        groupes = {"ARNt" : [], "Riboswitch" : [], "Ribozyme" : [], "Intron" : [], "SRP" : [], "ARNm" : []}
        compteur = 0

        for row in csvreader :
            if len(row) == 7 :
                arnr = False
                for cle in groupes_ARNr :
                    if cle in row[3] :
                        arnr = True
                        
                if arnr or "ribosomal" in row[3] or "Ribosomal" in row[3] or "RIBOSOMAL" in row[3] or "rRNA" in row[3] or "ribomosomal" in row[3] or "RRNA" in row[3] or "ribosome" in row[3] or "ribsomal" in row[3] or "ribosmal" in row[3] or "LARGE SUBUNIT" in row[3] :
                    ajoute = False
                    for cle in groupes_ARNr :
                        if cle in row[3] :
                            groupes_ARNr[cle].append((row[1], row[2]))
                            ajoute = True
                    if not ajoute :
                        if "alpha" in row[3] :
                            groupes_ARNr["LSU-alpha"].append((row[1], row[2]))
                        elif "beta" in row[3] :
                            groupes_ARNr["LSU-beta"].append((row[1], row[2]))
                    
                    #groupes["ARNr"].append((row[1], row[2]))
                elif "tRNA" in row[3] or "TRNA" in row[3] or "transfer RNA" in row[3] or "Transfer RNA" in row[3] or "trNA" in row[3] or "TRANSFER RNA" in row[3] :
                    groupes["ARNt"].append((row[1], row[2]))
                elif "riboswitch" in row[3] or "Riboswitch" in row[3] or "RIBOSWITCH" in row[3]:
                    groupes["Riboswitch"].append((row[1], row[2]))
                elif "ribozyme" in row[3] or "Ribozyme" in row[3] or "RIBOZYME" in row[3] or "Ribonuclease" in row[3] or "ribonuclease" in row[3] or "RIBONUCLEASE" in row[3] or "RNase" in row[3] or "Twister" in row[3] or "twister" in row[3] :
                    groupes["Ribozyme"].append((row[1], row[2]))
                elif "mRNA" in row[3] or "MRNA" in row[3] or "messenger RNA" in row[3] or "Messenger RNA" in row[3] or "MESSENGER RNA" in row[3] :
                    groupes["ARNm"].append((row[1], row[2]))
                elif "intron" in row[3] or "INTRON" in row[3] :
                    groupes["Intron"].append((row[1], row[2]))
                elif "SRP" in row[3] :
                    groupes["SRP"].append((row[1], row[2]))
                else :
                    print((row[1], row[2], row[3]))
            compteur += 1      
#         print(groupes)
#         for cle in groupes_ARNr :
#             print(cle)
#             print(len(groupes_ARNr[cle]))
#         print(len(groupes["ARNt"]))
#         print(len(groupes["Riboswitch"]))
#         print(len(groupes["Ribozyme"]))
#         print(len(groupes["ARNm"]))
#         print(len(groupes["Intron"]))
#         print(len(groupes["SRP"]))
#         print(compteur)
                 
    with open(rep_pickle, 'rb') as fichier_all_aminor :
        mon_depickler = pickle.Unpickler(fichier_all_aminor)
        all_aminor = mon_depickler.load()
        
            #print(all_aminor)
        compteur_tot = 0
        print(len(all_aminor))
        liste = []
        for cle in all_aminor.keys() :
            compteur = 1
            #if cle == '4tud' or cle == '4v81' :
            for graphe in all_aminor[cle] :
                #if cle == "1ond" :
                    num_ch = []
                    for noeud in graphe.nodes() :
                        if noeud[0] not in num_ch :
                            num_ch.append(noeud[0])
                    arn = False
                    for elt in num_ch :
                        if (cle.upper(), elt) in groupes[num_ARN] :
                            arn = True
                     
                    if arn : 
                        liste.append((cle, compteur))
                    compteur += 1
                    compteur_tot += 1
    print(liste)
    print(len(liste))
    print(len(all_aminor))
    print(compteur_tot)
    
    with open(NEW_EXTENSION_PATH_TAILLE+"%s/groupe_%s_%s_.pickle"%(typ_motif, typ_motif, num_ARN), 'wb') as fichier_sortie :
        mon_pickler = pickle.Pickler(fichier_sortie)
        mon_pickler.dump(liste)
    


### Recherche identiques  ###   


'''04/09/19
au sein des groupes d'homologues definis par methode sequence-clustalW (voir carnet bleu 08/10/19)
recherche les groupes d'identiques (meme numeros de sommets du motif (dans le graphe de motif ou dans le fichier PDB (champ fr3d)), memes nts, memes aretes)
version 2 de la recherche d'homologues et version 1 de la recherche d'identiques '''                    
def groupes_ARN_creation_homologues_et_identiques(num_ARN):
    with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
        mon_depickler = pickle.Unpickler(fichier_all_aminor)
        all_aminor = mon_depickler.load()
#     with open(NEW_EXTENSION_PATH_TAILLE+"groupe_%s.pickle"%num_ARN, 'rb') as fichier_sortie :
#         mon_depickler = pickle.Unpickler(fichier_sortie)
#         liste = mon_depickler.load()
#         
#         groupes_homologues = []              
#      
#         for i in range(len(liste)) :
# #         if "fichier_"+str(liste[i][0])+"_"+str(liste[i][1])+".pickle" in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
#             with open(NEW_EXTENSION_PATH_TAILLE + "fichier_"+str(liste[i][0])+"_"+str(liste[i][1])+".pickle", 'rb') as fichier_1 :
#                 mon_depickler_1 = pickle.Unpickler(fichier_1)
#                 extension1 = mon_depickler_1.load()
# #         else :
# #             with open(NEW_EXTENSION_PATH_TAILLE + "Doublons/fichier_"+str(liste[i][0])+"_"+str(liste[i][1])+".pickle", 'rb') as fichier_1 :
# #                 mon_depickler_1 = pickle.Unpickler(fichier_1)
# #                 extension1 = mon_depickler_1.load()
#             for j in range(i+1, len(liste)) :
#                     print(i)
#                     print(j)
#                     
# #                     if "fichier_"+str(liste[j][0])+"_"+str(liste[j][1])+".pickle" in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
#                     with open(NEW_EXTENSION_PATH_TAILLE + "fichier_"+str(liste[j][0])+"_"+str(liste[j][1])+".pickle", 'rb') as fichier_2 :
#                             mon_depickler_2 = pickle.Unpickler(fichier_2)
#                             extension2 = mon_depickler_2.load()
# #                     else :
# #                         with open(NEW_EXTENSION_PATH_TAILLE + "Doublons/fichier_"+str(liste[j][0])+"_"+str(liste[j][1])+".pickle", 'rb') as fichier_2 :
# #                             mon_depickler_2 = pickle.Unpickler(fichier_2)
# #                             extension2 = mon_depickler_2.load()
#                         #liste_a_faire.append((extension1, extension2))
# #                     if liste[j] in [('6cae', 26), ('6o97', 1)] and liste[i] in [('6cae', 26), ('6o97', 1)] :
# #                             print(extension1.nodes.data())
# #                             print(extension2.nodes.data())
# #                             print(homologue(extension1, extension2))
# #                             exit()
#                     if homologue(extension1, extension2) :
#                             
#                             ok1 = -1
#                             ok2 = -1
#                             c = 0
#                             for elt in groupes_homologues :
#                                 if liste[i] in elt and liste[j] not in elt :
#                                     ok1 = c
#                                 elif liste[j] in elt and liste[i] not in elt :
#                                     ok2 = c
#                                 if liste[j] in elt and liste[i] in elt :
#                                     ok1 = len(groupes_homologues)+1
#                                     ok2 = len(groupes_homologues)+1
#                                 c += 1
#                             
#                             if ok1 != len(groupes_homologues)+1 and ok2 != len(groupes_homologues)+1 :
#                                 if ok1 != -1 and ok2 != -1 :
#                                     groupes_homologues = fusion(groupes_homologues, ok1, ok2)
#                                 elif ok1 != -1 and ok2 == -1 :
#                                     groupes_homologues[ok1].append(liste[j])
#                                 elif ok2 != -1 and ok1 == -1 :    
#                                     groupes_homologues[ok2].append(liste[i])
#                                 else :
#                                     groupes_homologues.append([liste[i], liste[j]])
#     #                         ok = False
#     #                         for elt in groupes_homologues : 
#     #                             if liste[i] in elt and liste[j]  not in elt :
#     #                                 elt.append(liste[j])
#     #                                 ok = True
#     #                             elif liste[j] in elt and liste[i]  not in elt :
#     #                                 elt.append(liste[i])
#     #                                 ok = True
#     #                             elif liste[i] in elt and liste[j] in elt :
#     #                                 ok = True
#     #                         if not ok :
#     #                             groupes_homologues.append([liste[i], liste[j]])
#     #                         if liste[i] in [('4ujc', 2), ('4uje', 2)] and liste[j] in [('4ujc', 2), ('4uje', 2)] :
#     #                             print(ok1)
#     #                             print(ok2)
#     #                             print(homologue(extension1, extension2))
#     #                             print(groupes_homologues)
#     #                             exit()
#                          
#     print(groupes_homologues)
#     print(len(groupes_homologues))
#     
#     with open("groupes_%s_homologues.pickle"%num_ARN, 'wb') as fichier_homologues :
#         mon_pickler = pickle.Pickler(fichier_homologues)
#         mon_pickler.dump(groupes_homologues)
#     compteur = 1
#     for elt in groupes_28S_homologues :
#         print(len(elt))
#         #print(compteur)
#         compteur += 1

    with open("groupes_%s_homologues_sequences.pickle"%num_ARN, 'rb') as fichier_homologues :
        mon_depickler = pickle.Unpickler(fichier_homologues)
        groupes_homologues = mon_depickler.load()
        
        print(groupes_homologues)
        
        groupes_identiques = []
        for elt in groupes_homologues :
            
            groupes_identiques.append([])
            for i in range(len(elt)) :
                for j in range(i+1, len(elt)) :
                    print("petit rat")
                    e1 = elt[i]
                    e2 = elt[j]
                    print(i)
                    print(j)
                    for cle1 in all_aminor.keys() :
                        if cle1 == e1[0] :
                            compteur1 = 1
                            for graphe1 in all_aminor[cle1] :
                                if compteur1 == e1[1] :
                                    for cle2 in all_aminor.keys() :
                                        if cle2 == e2[0] :
                                            compteur2 = 1
                                            for graphe2 in all_aminor[cle2] :
                                                if compteur2 == e2[1] :
                                                    
                                                    compteur_idem = 0
                                                    for noeud1, data1 in graphe1.nodes(data=True) :
                                                          
                                                        for noeud2, data2 in graphe2.nodes(data=True) :
                                                             
                                                            if (noeud1[1] == noeud2[1] or data1["fr3d"] == data2["fr3d"])  and data1["nt"] == data2["nt"] and data1["real_nt"] == data2["real_nt"] :
                                                                compteur_idem += 1
                                                    #print(compteur_idem)
      
                                                    #print(idem)
                                                    compteur_idem_aretes = 0    
                                                      
                                                    #print(graphe1.edges.data())
                                                    #print(graphe2.edges.data())
                                                    for u1,v1, data1 in graphe1.edges(data=True) :
                                                        for u2,v2,data2 in graphe2.edges(data=True) :
                                                            if ((u1[1] == u2[1] and v1[1] == v2[1]) or (graphe1.nodes[u1]["fr3d"] == graphe2.nodes[u2]["fr3d"] and graphe1.nodes[v1]["fr3d"] == graphe2.nodes[v2]["fr3d"])) and graphe1.nodes[u1]["nt"] == graphe2.nodes[u2]["nt"] and graphe1.nodes[u1]["real_nt"] == graphe2.nodes[u2]["real_nt"] and graphe1.nodes[v1]["nt"] == graphe2.nodes[v2]["nt"] and graphe1.nodes[v1]["real_nt"] == graphe2.nodes[v2]["real_nt"] and data1["label"] == data2["label"] : #and data1["near"] == data2["near"] :
                                                                compteur_idem_aretes += 1
                                                    
                                                    print(compteur_idem)
                                                    print(compteur_idem_aretes)
                                                     
                                                    if compteur_idem  >= graphe1.number_of_nodes() and compteur_idem >= graphe2.number_of_nodes() and compteur_idem_aretes == graphe1.number_of_edges() and compteur_idem_aretes == graphe2.number_of_edges() :           
    #                                                     ok = False
    #                                                     for groupe in groupes_identiques[len(groupes_identiques)-1] :
    # #                                                         print(groupe) 
    #                                                         if e1 in groupe and e2 not in groupe :
    #                                                             groupe.append(e2)
    #                                                             ok = True
    #                                                         elif e2 in groupe and e1 not in groupe :
    #                                                             groupe.append(e1)
    #                                                             ok = True
    #                                                         elif e1 in groupe and e2 in groupe :
    #                                                             ok = True
    #                                                     if not ok :
    #                                                         groupes_identiques[len(groupes_identiques)-1].append([e1, e2])
                                                            
                                                        
                                                        ok1 = -1
                                                        ok2 = -1
                                                        c = 0
                                                        for groupe in groupes_identiques[len(groupes_identiques)-1] :
                                                            if e1 in groupe and e2 not in groupe :
                                                                ok1 = c
                                                            elif e2 in groupe and e1 not in groupe :
                                                                ok2 = c
                                                            if e2 in groupe and e1 in groupe :
                                                                ok1 = len(groupes_identiques[len(groupes_identiques)-1])+1
                                                                ok2 = len(groupes_identiques[len(groupes_identiques)-1])+1
                                                            c += 1
                                                        
                                                        if ok1 != len(groupes_identiques[len(groupes_identiques)-1])+1 and ok2 != len(groupes_identiques[len(groupes_identiques)-1])+1 :
                                                            if ok1 != -1 and ok2 != -1 :
                                                                groupes_identiques[len(groupes_identiques)-1] = fusion(groupes_identiques[len(groupes_identiques)-1], ok1, ok2)
                                                            elif ok1 != -1 and ok2 == -1 :
                                                                groupes_identiques[len(groupes_identiques)-1][ok1].append(e2)
                                                            elif ok2 != -1 and ok1 == -1 :    
                                                                groupes_identiques[len(groupes_identiques)-1][ok2].append(e1)
                                                            else :
                                                                groupes_identiques[len(groupes_identiques)-1].append([e1, e2])
                                                compteur2 += 1
                                      
                                compteur1 += 1
    print("groupes identiques")             
    print(groupes_identiques)
    print(len(groupes_identiques))
    
#     compteur = 0
#     for elt in groupes_homologues :
#         print(elt)
#         print(groupes_identiques[compteur])
#         compteur += 1
        

        
    with open("groupes_%s_identiques_sequences.pickle"%num_ARN, 'wb') as fichier_identiques :
        mon_pickler = pickle.Pickler(fichier_identiques)
        mon_pickler.dump(groupes_identiques)
        
    print("groupes identiques")             
    print(groupes_identiques)
    print(len(groupes_identiques))  
#     for elt in groupes_28S_identiques : 
#         print(len(elt))
                     

''' 15/10/19 
recherche les occurrences du type d'ARN passe en parametre qui ont exactement les memes sequences (30 nts autour du motif pour les 3 brins)
et les stocke dans fichier pickle sous forme d'une liste de listes d'occurrences identiques (num_PDB, num_occ)
version finale de la recherche d'occurrences d'identiques
num_ARN : type d'ARN
type_motif : type de motif (A-minor, motif 6left...)'''
def filtrage_sequences_identiques(numARN, type_motif):
    with open(NEW_EXTENSION_PATH_TAILLE+"%s/groupe_%s_%s_.pickle"%(type_motif, type_motif, numARN), 'rb') as fichier :
        mon_depickler = pickle.Unpickler(fichier)
        groupe = mon_depickler.load()
        print(groupe)
        print(len(groupe))
        deja_range = []
        
        groupes_identiques = []
        for i in range(len(groupe)) :
            #print(groupes_identiques)
            elt1 = groupe[i]
            
#             deja_vu = False
#             for groupe_id in groupes_identiques :
#                 if elt1 in groupe_id :
#                     deja_vu = True
                    
#             if not deja_vu :
            if elt1 not in deja_range :
                print(i)
                groupes_identiques.append([elt1])
                deja_range.append(elt1)
                with open(NEW_EXTENSION_PATH_TAILLE+"%s/fichier_%s_%s.pickle"%(type_motif, elt1[0], str(elt1[1])), "rb") as fichier_extension :
                    mon_depickler_ext = pickle.Unpickler(fichier_extension)
                    extension = mon_depickler_ext.load()
                    seq11, seq12, seq13 = recup_sequences_autour_motif(elt1[0], extension, 30)
                    
                    for j in range(i+1, len(groupe)) :
                        #print(groupes_identiques)
                        elt2 = groupe[j]
                        
#                         deja_vu = False
#                         for groupe_id in groupes_identiques :
#                             if elt2 in groupe_id :
#                                 deja_vu = True
#                                 
#                         if not deja_vu :
                        if elt2 not in deja_range :
                            print(j)
                            
                            with open(NEW_EXTENSION_PATH_TAILLE+"%s/fichier_%s_%s.pickle"%(type_motif, elt2[0], str(elt2[1])), "rb") as fichier_extension :
                                mon_depickler_ext = pickle.Unpickler(fichier_extension)
                                extension = mon_depickler_ext.load()
                                seq21, seq22, seq23 = recup_sequences_autour_motif(elt2[0], extension, 30)
            
                                if seq11 == seq21 and seq12 == seq22 and seq13 == seq23 :
                                    for groupe_id in groupes_identiques :
                                        if elt1 in groupe_id :
                                            groupe_id.append(elt2)
                                            deja_range.append(elt2)
                
        print(groupes_identiques)
        print(len(groupes_identiques)) 
        
        with open(NEW_EXTENSION_PATH_TAILLE+"%s/groupes_identiques_%s_%s.pickle"%(type_motif, type_motif, numARN), 'wb') as fichier_identiques :
            mon_pickler = pickle.Pickler(fichier_identiques)
            mon_pickler.dump(groupes_identiques)   
        
'''
verifie les groupes d'identiques obtenus par recherche de sequence identique autour du motif (30 nts)
en reregardant les sequences
num_ARN : type d'ARN
type_motif : type de motif (A-minor, motif 6left...)
'''        
def verif_filtrage(numARN, type_motif):
    with open(NEW_EXTENSION_PATH_TAILLE+"%s/groupes_identiques_%s_%s.pickle"%(type_motif, type_motif, numARN), 'rb') as fichier_identiques :
        mon_depickler = pickle.Unpickler(fichier_identiques)
        groupes_identiques = mon_depickler.load() 
            
        with open("occ_multi_chaine.pickle", 'rb') as fichier_multi_chaines :
            mon_depickler_2 = pickle.Unpickler(fichier_multi_chaines)
            liste_plusieurs_chaines = mon_depickler_2.load() 
            
            '''supprimer les occurrences entre plusieurs chaines car on les met a part'''
            somme = 0
            compteur = 1
            new_groupes_identiques = []
            for groupe in groupes_identiques :
                new_groupe = []
                for elt in groupe :
                    if elt not in liste_plusieurs_chaines : 
#                         print(compteur)
#                         print("ramou")
                        new_groupe.append(elt)
                        
                new_groupes_identiques.append(new_groupe)
                #print(len(groupe))
                compteur += 1
                somme += len(groupe)
                    
            print(somme)
            
            somme = 0
            for groupe in new_groupes_identiques :
                for i in range(len(groupe)) :
                    with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(groupe[i][0], str(groupe[i][1])), "rb") as fichier_extension :
                        mon_depickler_ext = pickle.Unpickler(fichier_extension)
                        extension1 = mon_depickler_ext.load()
                        seq11,seq12,seq13 = recup_sequences_autour_motif(groupe[i][0], extension1)
                        for j in range(i+1, len(groupe)) :
                            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(groupe[j][0], str(groupe[j][1])), "rb") as fichier_extension :
                                mon_depickler_ext = pickle.Unpickler(fichier_extension)
                                extension2 = mon_depickler_ext.load()
                                seq21,seq22,seq23 = recup_sequences_autour_motif(groupe[j][0], extension2)
                                
                                if seq11 != seq21 or seq12 != seq22 or seq13 != seq23 :
                                    print("gros rat")
                    
                somme += len(groupe)
            print(somme)
        
            
    with open(NEW_EXTENSION_PATH_TAILLE+"groupes_identiques_%s.pickle"%numARN, 'wb') as fichier_identiques :
        mon_depickler = pickle.Pickler(fichier_identiques)
        mon_depickler.dump(new_groupes_identiques)  
        
''' recherche d'un representant dans chaque groupe d'occurrences identiques (version 2 de la recherche d'identiques)
- premier critere : meilleure resolution
- si egalite : le plus grand nombre d'annotation d'aretes
- si egalite : la plus grande sequence
stocke les representants dans un fichier pickle 
num_ARN : type d'ARN
type_motif : type de motif (A-minor, motif 6left...)
 '''
def recherche_representant(numARN, type_motif):
    with open(NEW_EXTENSION_PATH_TAILLE+"%s/groupes_identiques_%s_%s.pickle"%(type_motif, type_motif, numARN), 'rb') as fichier_identiques :
        mon_depickler = pickle.Unpickler(fichier_identiques)
        groupes_identiques = mon_depickler.load()
        
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
            print(numARN)
            print(resolutions)
            
            liste_representant = []
            compter1 = 0
            compter2 = 0
            compter3 = 0
            compter4 = 0
             
            for groupe in groupes_identiques :
                min_res = 50.0
                liste_min_res = []
                for elt in groupe :
                    print(resolutions[elt[0]])
                    if resolutions[elt[0]] != None :
                        if resolutions[elt[0]] < min_res :
                            min_res = resolutions[elt[0]]
                    else :
                        print("petit rat")
                        print(elt)
                print(min_res)
                for elt in groupe :
                    if resolutions[elt[0]] != None :
                        if resolutions[elt[0]] == min_res :
                            liste_min_res.append(elt)    
                     
                print(min_res)   
                print(liste_min_res)    
                
                if len(liste_min_res) > 1 :
                    maxi_nb_aretes = 0
                    liste_maxi_aretes = []
                    for elt in liste_min_res :
                        with open(NEW_EXTENSION_PATH_TAILLE+"%s/fichier_%s_%s.pickle"%(type_motif, elt[0],elt[1]), 'rb') as fichier_extension :
                            mon_depickler = pickle.Unpickler(fichier_extension)
                            extension = mon_depickler.load()
                             
                            nb_aretes = 0
                            for u,v,data in extension.edges(data=True) :
                                if data["near"] == False :
                                    nb_aretes += 1
                             
                            if nb_aretes > maxi_nb_aretes :
                                maxi_nb_aretes = nb_aretes
                    for elt in liste_min_res :
                        with open(NEW_EXTENSION_PATH_TAILLE+"%s/fichier_%s_%s.pickle"%(type_motif, elt[0],elt[1]), 'rb') as fichier_extension :
                            mon_depickler = pickle.Unpickler(fichier_extension)
                            extension = mon_depickler.load()
                             
                            nb_aretes = 0
                            for u,v,data in extension.edges(data=True) :
                                if data["near"] == False :
                                    nb_aretes += 1
                             
                            if nb_aretes == maxi_nb_aretes :
                                liste_maxi_aretes.append(elt)
                 
                    print(maxi_nb_aretes)
                    print(liste_maxi_aretes)
                     
                    if len(liste_maxi_aretes) > 1 :
                        maxi_taille_seq = 0
                        liste_maxi_taille_seq = []
                        for elt in liste_maxi_aretes :
                            with open("Graphs/%s.pickle"%elt[0], 'rb') as fichier_graphe :
                                mon_depickler = pickle.Unpickler(fichier_graphe)
                                graphe = mon_depickler.load()
                                 
                                if graphe.number_of_nodes() > maxi_taille_seq :
                                    maxi_taille_seq = graphe.number_of_nodes()
                         
                        for elt in liste_maxi_aretes :
                            with open("Graphs/%s.pickle"%elt[0], 'rb') as fichier_graphe :
                                mon_depickler = pickle.Unpickler(fichier_graphe)
                                graphe = mon_depickler.load()
                                 
                                if graphe.number_of_nodes() == maxi_taille_seq :
                                    liste_maxi_taille_seq.append(elt)
                                     
                        print(maxi_taille_seq)
                        print(liste_maxi_taille_seq)
                         
                        if len(liste_maxi_taille_seq) > 1 :
                            liste_representant.append(liste_maxi_taille_seq[0])
                            print("tout petit rat")
                            compter4 += 1
                        else :
                            liste_representant.append(liste_maxi_taille_seq[0])
                            compter3 += 1
                         
                    else :
                        
                        liste_representant.append(liste_maxi_aretes[0])
                        
                        compter2 += 1
                else :
                    if len(liste_min_res) != 0 :
                        liste_representant.append(liste_min_res[0])
                    else :
                        liste_representant.append(groupe[0])
                    compter1 += 1
                    
                
        print(compter1)
        print(compter2)
        print(compter3)
        print(compter4)     
        
        print(len(liste_representant))
        
        with open("Nouvelles_donnees/%s/liste_representant_%s.pickle"%(type_motif,numARN), 'wb') as fichier_sortie :
            mon_pickler = pickle.Pickler(fichier_sortie)
            mon_pickler.dump(liste_representant)   


### Recherche homologues ###

## Version ClustalW ##
''' 04/10/19
dans les groupes d'homologues version 1 (positions similaires du motif à 25 nts pres),
- effectue l'alignement par clustalW de toutes les paires d'occurrences au sein de chaque groupe
- cree des sous-groupes d'elts dont la valeur de ressemblance est sup à 90% (les ajoute un par un, donc depend de l'ordre de traitement des paires)
num_ARN : type d'ARN
rep : chemin du repertoire de stockage des fichiers ClustalW
stocke les alignements multiples des sequences de chaque groupe dans des fichiers
 '''
def verif_homologues(num_ARN, rep):
    with open("groupes_%s_homologues.pickle"%num_ARN, 'rb') as fichier_homologues :
        mon_depickler_1 = pickle.Unpickler(fichier_homologues)
        groupes_homologues = mon_depickler_1.load()  
        compteur = 0
        for groupe in groupes_homologues :
            
            #if compteur <= 2 :
                liste_groupe = []
                compteur_sous_groupes = 1
                while len(liste_groupe) < len(groupe) :   
                    new_groupe = [x for x in groupe if x not in liste_groupe]
                    
                    
                    compteur_elt = 1
                    
                    if "fichier_fasta_1_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes) in os.listdir(rep) :
                        os.remove(rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes))
                    fichier_fasta_1 = open(rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), "w") 
                    fichier_fasta_1.close()
                    
                    if "fichier_fasta_2_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes) in os.listdir(rep) :
                        os.remove(rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes))
                    fichier_fasta_2 = open(rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), "w") 
                    fichier_fasta_2.close()
                    for elt in new_groupe :
                        with open(rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), 'r') as fichier_fasta_1 :
                            with open(rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), 'r') as fichier_fasta_2 :
                                with open(rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes), 'w') as fichier_fasta_1_temp :
                                    with open(rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes), 'w') as fichier_fasta_2_temp :
                                        lignes = fichier_fasta_1.readlines()
                                        #print(fichier_fasta_1.readlines())
                                        for ligne in lignes :
                                            fichier_fasta_1_temp.write(ligne)
                                        
                                        lignes = fichier_fasta_2.readlines()
                                        #print(fichier_fasta_1.readlines())
                                        for ligne in lignes :
                                            fichier_fasta_2_temp.write(ligne)
                                            
                                        print(elt)
                                        with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(elt[0], elt[1]), 'rb') as fichier_extension :
                                            mon_depickler = pickle.Unpickler(fichier_extension)
                                            extension = mon_depickler.load()
                                            seq1, seq2 = recup_sequences_autour_motif(elt[0], extension)
                                                    
                                            fichier_fasta_1_temp.write(">elt_%s_%s\n"%elt)
                                            fichier_fasta_1_temp.write(seq1+"\n")
                                                    
                                            fichier_fasta_2_temp.write(">elt_%s_%s\n"%elt)
                                            fichier_fasta_2_temp.write(seq2+"\n")
                                            
                        if compteur_elt > 1 :
                            with open(rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), 'a') as fichier_fasta_1 :
                                with open(rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), 'a') as fichier_fasta_2 :
                                    out_file_1 = rep+"/aligned_1_temp.aln"
                                    out_file_2 = rep+"/aligned_2_temp.aln"       
                                    # Get the command for Clustal Omega
                                    clustalw_cline_1 = ClustalwCommandline(infile=rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes), outfile=out_file_1, output="CLUSTAL")
                                    clustalw_cline_2 = ClustalwCommandline(infile=rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes), outfile=out_file_2, output="CLUSTAL")
            
                                    # Print the executable command
                                    print(clustalw_cline_1)
                                    stdout_1, stderr_1 = clustalw_cline_1()
                                    print(stdout_1)
                                    stdout_2, stderr_2 = clustalw_cline_2()
                                    print(stdout_2)
        #                             with open("stdout_clustalw.txt", "w") as fichier_stdout :
        #                                 fichier_stdout.write(stdout)
        #                             print(type(stdout))
                                            
                                            #nb_seq = len(groupe)
                                    lignes_1 = stdout_1.split("\n")
                                            #print(lignes)
                                    ligne_1 = lignes_1[0]
                                    compte_lignes = 0
                                    pas_bon_1 = False
                                    while compte_lignes < len(lignes_1)  :
                                        if "Aligned." in ligne_1 :
                                            if int(ligne_1.split(" ")[5]) < 90 :
                                                pas_bon_1 = True
                                        ligne_1 = lignes_1[compte_lignes]
                                        compte_lignes += 1
                                        
                                    lignes_2 = stdout_2.split("\n")
                                    ligne_2 = lignes_1[0]
                                    compte_lignes = 0
                                    pas_bon_2 = False
                                    while compte_lignes < len(lignes_2)  :
                                        if "Aligned." in ligne_2 :
                                            if int(ligne_2.split(" ")[5]) < 90 :
                                                pas_bon_2 = True
                                        ligne_2 = lignes_2[compte_lignes]
                                        compte_lignes += 1
                                                
                                    if not pas_bon_1 and not pas_bon_2:
                                        fichier_fasta_1.write(">elt_%s_%s\n"%elt)
                                        fichier_fasta_1.write(seq1+"\n")
                                        fichier_fasta_2.write(">elt_%s_%s\n"%elt)
                                        fichier_fasta_2.write(seq2+"\n")
                                        liste_groupe.append(elt)
                                        
                        
                    #                     if compte_lignes != len(lignes) :   
                    #                         print(ligne) 
                    #                         print(ligne.split(" "))
                    #                         print(ligne.split(" ")[5])
                                            
                                            
                        compteur_elt += 1
                    
                    if len(new_groupe) != 1 :
                        
                        print(liste_groupe)
                        print(len(liste_groupe))
                        
                        if "fichier_fasta_1_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes) in os.listdir(rep) :
                            os.remove(rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes))
                        if "fichier_fasta_2_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes) in os.listdir(rep) :
                            os.remove(rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d_temp.fa"%(compteur, compteur_sous_groupes))
                        if "fichier_fasta_1_groupe_%d_sous_groupe_%d_temp.dnd"%(compteur, compteur_sous_groupes) in os.listdir(rep) :
                            os.remove(rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d_temp.dnd"%(compteur, compteur_sous_groupes))
                        if "fichier_fasta_2_groupe_%d_sous_groupe_%d_temp.dnd"%(compteur, compteur_sous_groupes) in os.listdir(rep) :
                            os.remove(rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d_temp.dnd"%(compteur, compteur_sous_groupes))
                        
                        
                        out_file_1 = rep+"/aligned_1_groupe_%d_sous_groupe_%d.aln"%(compteur,compteur_sous_groupes)
                        out_file_2 = rep+"/aligned_2_groupe_%d_sous_groupe_%d.aln"%(compteur,compteur_sous_groupes)       
                        # Get the command for Clustal Omega
                        clustalw_cline_1 = ClustalwCommandline(infile=rep+"/fichier_fasta_1_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), outfile=out_file_1, output="CLUSTAL")
                        clustalw_cline_2 = ClustalwCommandline(infile=rep+"/fichier_fasta_2_groupe_%d_sous_groupe_%d.fa"%(compteur, compteur_sous_groupes), outfile=out_file_2, output="CLUSTAL")
        
                        # Print the executable command
                        print(clustalw_cline_1)
                        stdout_1, stderr_1 = clustalw_cline_1()
                        print(stdout_1)
                        stdout_2, stderr_2 = clustalw_cline_2()
                        
                        with open(rep+"/stdout_clustalw_1_groupe_%d_sous_groupe_%d.txt"%(compteur,compteur_sous_groupes), "w") as fichier_stdout :
                            fichier_stdout.write(stdout_1)
                        with open(rep+"/stdout_clustalw_2_groupe_%d_sous_groupe_%d.txt"%(compteur,compteur_sous_groupes), "w") as fichier_stdout :
                            fichier_stdout.write(stdout_2)
                    else :
                        liste_groupe.append(new_groupe[0])
                    
                    compteur_sous_groupes += 1
                    
    #                 if compteur_elt == 10 :
    #                     break
                            
                        
                
    #                     nb_seq = len(groupe)
    #                     lignes = stdout.split("\n")
    #                     #print(lignes)
    #                     ligne = lignes[0]
    #                     compte_lignes = 0
    #                     while compte_lignes < len(lignes)  :
    #                         if "Aligned."in ligne :
    #                             if int(ligne.split(" ")[5]) < 90 :
    #                                 print(ligne) 
    #                         ligne = lignes[compte_lignes]
    #                         compte_lignes += 1
    #                     if compte_lignes != len(lignes) :   
    #                         print(ligne) 
    #                         print(ligne.split(" "))
    #                         print(ligne.split(" ")[5])
                            
                            
                compteur += 1    


''' 14/10/19
dans les groupes d'homologues version 1 (positions similaires du motif à 25 nts pres),
- effectue l'alignement par clustalW de toutes les paires d'occurrences au sein de chaque groupe
- cree des sous-groupes d'elts dont la valeur de ressemblance est sup à 90% (fait un graphe de toutes les occurrences et met une arete entre deux occurrences si la ressemblance entre les deux est sup à 90%, puis fait des composantes connexes du graphe les groupes d'homologues)
num_ARN : type d'ARN
rep : chemin du repertoire de stockage des fichiers ClustalW
stocke les alignements multiples des sequences de chaque groupe dans des fichiers
'''
def verif_homologues_version_2(num_ARN, rep):
    with open("groupes_%s_homologues.pickle"%num_ARN, 'rb') as fichier_homologues :
        mon_depickler_1 = pickle.Unpickler(fichier_homologues)
        groupes_homologues = mon_depickler_1.load()  
        compteur = 0
        new_groupes_homologues = []
        for groupe in groupes_homologues :
            #if compteur >= 6 :
                graph = nx.Graph()
                ''' Initialisation des sommets du graphe '''
                for i in range(len(groupe)) :
                    graph.add_node(i)
            
                for i in range(len(groupe)) :
                    with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(groupe[i][0], groupe[i][1]), 'rb') as fichier_extension :
                        mon_depickler = pickle.Unpickler(fichier_extension)
                        extension = mon_depickler.load()
                        seqi1, seqi2 = recup_sequences_autour_motif(groupe[i][0], extension)
                        
                        for j in range(i+1, len(groupe)) :
                            with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(groupe[j][0], groupe[j][1]), 'rb') as fichier_extension :
                                mon_depickler = pickle.Unpickler(fichier_extension)
                                extension = mon_depickler.load()
                                seqj1, seqj2 = recup_sequences_autour_motif(groupe[j][0], extension)
                                
                                
                                with open(rep+"/fichier_fasta_1_groupe_%d_test_graphe.fa"%(compteur), 'w') as fichier_fasta_1 :
                                    with open(rep+"/fichier_fasta_2_groupe_%d_test_graphe.fa"%(compteur), 'w') as fichier_fasta_2 :
                                        
                                        fichier_fasta_1.write(">elt_%s_%s\n"%groupe[i])
                                        fichier_fasta_1.write(seqi1+"\n")
                                        fichier_fasta_2.write(">elt_%s_%s\n"%groupe[i])
                                        fichier_fasta_2.write(seqi2+"\n")
                                        
                                        
                                        fichier_fasta_1.write(">elt_%s_%s\n"%groupe[j])
                                        fichier_fasta_1.write(seqj1+"\n")
                                        fichier_fasta_2.write(">elt_%s_%s\n"%groupe[j])
                                        fichier_fasta_2.write(seqj2+"\n")
                                        
                                out_file_1 = rep+"/aligned_1_temp.aln"
                                out_file_2 = rep+"/aligned_2_temp.aln"       
                                            # Get the command for Clustal Omega
                                            
                                
                                clustalw_cline_1 = ClustalwCommandline(infile=rep+"/fichier_fasta_1_groupe_%d_test_graphe.fa"%(compteur), outfile=out_file_1, output="CLUSTAL")
                                clustalw_cline_2 = ClustalwCommandline(infile=rep+"/fichier_fasta_2_groupe_%d_test_graphe.fa"%(compteur), outfile=out_file_2, output="CLUSTAL")
            
                                # Print the executable command
                                #print(clustalw_cline_1)
                                stdout_1, stderr_1 = clustalw_cline_1()
                                #print(stdout_1)
                                stdout_2, stderr_2 = clustalw_cline_2()
                                #print(stdout_2)
    #                             with open("stdout_clustalw.txt", "w") as fichier_stdout :
    #                                 fichier_stdout.write(stdout)
    #                             print(type(stdout))
                                        
                                        #nb_seq = len(groupe)
                                lignes_1 = stdout_1.split("\n")
                                        #print(lignes)
                                ligne_1 = lignes_1[0]
                                compte_lignes = 0
                                pas_bon_1 = False
                                while compte_lignes < len(lignes_1)  :
                                    if "Aligned." in ligne_1 :
                                        if int(ligne_1.split(" ")[5]) < 90 :
                                            pas_bon_1 = True
                                    ligne_1 = lignes_1[compte_lignes]
                                    compte_lignes += 1
                                    
                                lignes_2 = stdout_2.split("\n")
                                ligne_2 = lignes_1[0]
                                compte_lignes = 0
                                pas_bon_2 = False
                                while compte_lignes < len(lignes_2)  :
                                    if "Aligned." in ligne_2 :
                                        if int(ligne_2.split(" ")[5]) < 90 :
                                            pas_bon_2 = True
                                    ligne_2 = lignes_2[compte_lignes]
                                    compte_lignes += 1
                                            
                                if not pas_bon_1 and not pas_bon_2:
                                    graph.add_edge(i,j)
                
                #print(graph.nodes.data())
                #print(graph.edges.data())
                #print(graph.number_of_edges())
                
                deja_vu = []
                composantes_connexes = []
                for noeud in graph.nodes() :
                    print(groupe[noeud])
                    if noeud not in deja_vu :
                        composantes_connexes.append([noeud])
                        deja_vu.append(noeud)
                         
                        #parcours en largeur
                        file_sommets = [noeud]
                        while len(file_sommets) > 0:
                            sommet_courant = file_sommets.pop(0)
                            enfants_courant = graph[sommet_courant]  
                            #print(len(file_sommets))
                            for enfant in enfants_courant:
                                if enfant not in deja_vu :
                                    composantes_connexes[len(composantes_connexes)-1].append(enfant)
                                    file_sommets.append(enfant)
                                deja_vu.append(enfant)
                
                for composante in composantes_connexes :
                    new_groupes_homologues.append([])
                    for elt in composante :
                        new_groupes_homologues[len(new_groupes_homologues)-1].append(groupe[elt])
                        
                os.remove(rep+"/fichier_fasta_1_groupe_%d_test_graphe.fa"%(compteur))
                os.remove(rep+"/fichier_fasta_2_groupe_%d_test_graphe.fa"%(compteur))
                os.remove(rep+"/fichier_fasta_1_groupe_%d_test_graphe.dnd"%(compteur))
                os.remove(rep+"/fichier_fasta_2_groupe_%d_test_graphe.dnd"%(compteur))
                os.remove(rep+"/aligned_1_temp.aln")
                os.remove(rep+"/aligned_2_temp.aln")
                #print(composantes_connexes)        
                
                #if compteur == 4 :
                #break         
                compteur += 1    
                
        print(new_groupes_homologues)
         
        with open("new_groupes_homologues_%s.pickle"%num_ARN, "wb") as fichier_pickle :
            mon_pickler = pickle.Pickler(fichier_pickle)
            mon_pickler.dump(new_groupes_homologues)

''' 09/10/19
a partir des fichiers d'alignements multiples obtenus par ClustalW
determine les groupes d'homologues avec le critere sur la sequence en plus 
et les stocke dans un fichier pickle sous forme d'une liste de listes
probleme de doublons dans les groupes je crois'''
def regrouper_par_homologues_critere_sequence(numARN):      
    with open("groupes_%s_homologues.pickle"%numARN, 'rb') as fichier_homologues :
        mon_depickler_1 = pickle.Unpickler(fichier_homologues)
        groupes_homologues = mon_depickler_1.load() 
        
        liste_deja_vu = {}
        new_groupes_homologues = []
        compteur = 0
        for groupe in groupes_homologues :
            
            new_groupes_temp = []
            for fic in os.listdir("Nouvelles_donnees/Groupes_homologues/%s"%numARN) :
                
                if "fichier_fasta_1_groupe_%d"%compteur in fic and "dnd" not in fic and "temp" not in fic :
                    with open("Nouvelles_donnees/Groupes_homologues/%s/"%numARN+fic, "r") as fichier_fasta_1 :
                        with open("Nouvelles_donnees/Groupes_homologues/%s/fichier_fasta_2_%s"%(numARN,fic[16:]), "r") as fichier_fasta_2 :
                            print(fic)
                            sous_groupe = []
                            lignes_1 = fichier_fasta_1.readlines()
                            lignes_2 = fichier_fasta_2.readlines() # juste pour verifier
                            
                            compter_lignes = 0
                            for ligne in lignes_1 :
                                if '>' in ligne :
                                    if lignes_2[compter_lignes] != ligne :
                                        print("bizarre")
                                        
                                    if (ligne.split(">")[1].split("_")[1], int(ligne.split(">")[1].split("_")[2])) not in liste_deja_vu.values() :
                                        liste_deja_vu.update({(compteur, fic.split("_")[7][:len(fic.split("_")[7])-2]) : (ligne.split(">")[1].split("_")[1], int(ligne.split(">")[1].split("_")[2]))})
                                        sous_groupe.append((ligne.split(">")[1].split("_")[1], int(ligne.split(">")[1].split("_")[2])))
                                    
                                    else :
                                        print("deux fois")
                                        print(compteur)
                                        print(fic.split("_")[7][:len(fic.split("_")[7])-2])
                                        print(list(liste_deja_vu.keys())[list(liste_deja_vu.values()).index((ligne.split(">")[1].split("_")[1], int(ligne.split(">")[1].split("_")[2])))])
                                
                                compter_lignes += 1              
                            
                            if len(sous_groupe) > 0 :
                                new_groupes_temp.append(sous_groupe)
                                
                            
                            #print(len(sous_groupe))
#             print(len(new_groupes_temp))
#             print(new_groupes_temp)                                       
            for elt in groupe :
                ok = False
                for new_groupe in new_groupes_temp :
                    if elt in new_groupe :
                        ok = True
                if not ok :
                    print(compteur)
                    print(elt)
                    new_groupes_temp.append([elt])
#             print(len(new_groupes_temp))
#             print(new_groupes_temp) 
            new_groupes_homologues.extend(new_groupes_temp)
                    
            compteur += 1
        
        print(len(groupes_homologues))   
        print(len(new_groupes_homologues)) 
        
        with open("groupes_%s_homologues_sequences.pickle"%numARN, 'wb') as fichier_homologues_sequence :
            mon_pickler = pickle.Pickler(fichier_homologues_sequence)
            mon_pickler.dump(new_groupes_homologues)                             

## Version BLAST ##

''' effectue BLAST sur les sequences des 3 brins pour toutes les paires d'occurrences d'un type d'ARN
stocke le resultat dans un fichier pickle '''
def rechercher_homologues_par_blast(numARN):
    with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%numARN, 'rb') as fichier_repr :
            mon_depickler = pickle.Unpickler(fichier_repr)
            liste_representant = mon_depickler.load()  
            
            liste_trouve = {}
            compteur = 0
            for i in range(len(liste_representant)) :
                if liste_representant[i] in [('5ngm', 15),('4ybb', 46),('4y4o', 58)] :
                #if i > 0 :
                    print(i)
                    with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(liste_representant[i][0], str(liste_representant[i][1])), "rb") as fichier_extension :
                        mon_depickler_ext = pickle.Unpickler(fichier_extension)
                        extension1 = mon_depickler_ext.load()
                        seq1 = recup_sequences_autour_motif(liste_representant[i][0], extension1, 30)
                        
                       # print(extension1.nodes.data())
                            
    #                         for noeud,data in graphe.nodes(data=True) :
    #                             if noeud[0] =='2' : 
    #                                 print(noeud)
    #                                 print(data)
                            #print(graphe.nodes.data())
                        
                        for j in range(i+1, len(liste_representant)) :
                                if liste_representant[j] in [('5ngm', 15),('4ybb', 46),('4y4o', 58)] :
                                    print(j)
                                    with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(liste_representant[j][0], str(liste_representant[j][1])), "rb") as fichier_extension :
                                        mon_depickler_ext = pickle.Unpickler(fichier_extension)
                                        extension2 = mon_depickler_ext.load()
                                        seq2 = recup_sequences_autour_motif(liste_representant[j][0], extension2, 30)
                                    
    #                                 print(liste_representant[i])
    #                                 print(liste_representant[j])
                                    
#                                     print(seq11)
#                                     print(seq21)
#                                     print(len(seq11))
#                                     print(len(seq21))
                                    print(liste_representant[i])
                                    print(liste_representant[j])
                                    for k in range(0, 3) :
                                        with open("fichier_seq1.fa", 'w') as fichier1 :
                                            fichier1.write(seq1[k])
                                        with open("fichier_seq2.fa", 'w') as fichier2 :
                                            fichier2.write(seq2[k])
                                        
                                        print(seq1[k])
                                        print(seq2[k])
                                        # Run BLAST and parse the output as XML
                                        output = NcbiblastnCommandline(query="fichier_seq1.fa", subject="fichier_seq2.fa", outfmt=5, task='megablast')()[0]
                                        blast_result_record = NCBIXML.read(StringIO(output))
                                        print(blast_result_record.alignments)
                                        #print(blast_result_record) 
                                        if len(blast_result_record.alignments) > 0  :
                                            #if (liste_representant[i], liste_representant[j]) not in liste_trouve.keys() :
                                            liste_trouve.update({(liste_representant[i], liste_representant[j]) : [{'num' : k, 'length' : 0, 'e-value' : 0, 'score' : 0}]} )
#                                             else :
#                                                 liste_trouve[(liste_representant[i], liste_representant[j])].append({'num' : i, 'length' : 0, 'e-value' : 0, 'score' : 0})
#                                             # Print some information on the result
                                            #compter = 
                                            for alignment in blast_result_record.alignments:
                                                compter  = 0
                                                for hsp in alignment.hsps:
                                                    if compter > 0 :
                                                        liste_trouve[(liste_representant[i], liste_representant[j])].append({'length' : 0, 'e-value' : 0, 'score' : 0})
                                                    print('****Alignment****')
                                                    print('sequence:', alignment.title)
                                                    print('length:', alignment.length)
                                                    print('e value:', hsp.expect)
                                                    print('score :',hsp.score)
                                                    print(hsp.query)
                                                    print(hsp.match)
                                                    print(hsp.sbjct)
                                                    liste_trouve[(liste_representant[i], liste_representant[j])][compter]["length"] = alignment.length
                                                    liste_trouve[(liste_representant[i], liste_representant[j])][compter]["e-value"] = hsp.expect
                                                    liste_trouve[(liste_representant[i], liste_representant[j])][compter]["score"] = hsp.score
                                                    compter += 1
                                                    #exit(0)
                                            compteur += 1
                                        
                                    
                                    
            print(compteur)
            print(liste_trouve)
            
#             with open("/media/coline/Maxtor/liste_homologie_par_blast_%s.pickle"%numARN, 'wb') as fichier_liste_trouve :
#                 mon_pickler = pickle.Pickler(fichier_liste_trouve)
#                 mon_pickler.dump(liste_trouve)

''' cree un graphe des occurrences d'un type d'ARN avec comme ponderation d'aretes les scores de BLAST obtenus
(pas d'aretes si BLAST n'a pas trouve de correspondance) '''
def definit_groupes_homologues_par_blast(numARN):
    with open("/media/coline/Maxtor/liste_homologie_par_blast_%s.pickle"%numARN, 'rb') as fichier_liste_trouve :
        mon_depickler = pickle.Unpickler(fichier_liste_trouve)
        liste_trouve = mon_depickler.load()
        
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%numARN, 'rb') as fichier_repr :
            mon_depickler = pickle.Unpickler(fichier_repr)
            liste_representant = mon_depickler.load()  
            
            graph_blast = nx.Graph()
            
            for elt in liste_representant :
                graph_blast.add_node(elt)
            
            #print(graph_blast.nodes.data())
            for cle in liste_trouve.keys() :
                #print(cle)
                #print(liste_trouve[cle])
                if liste_trouve[cle][0]['length'] == 62 :
                    if cle[0] in [('4y4o', 51), ('4y4p', 20), ('4y4o', 39)] and cle[1] in [('4y4o', 51), ('4y4p', 20), ('4y4o', 39)] :
                        print(cle)
                        print(liste_trouve[cle])
                    graph_blast.add_edge(cle[0], cle[1], poids=liste_trouve[cle][0]['score'])
            node_labels=dict([(u, (u)) for u,d in graph_blast.nodes(data=True)]) #else (u, (u)) for u,d in G.nodes(data=True) ])
            cliques = nx.find_cliques(graph_blast)
            
#             for clique in cliques :
#                 print(clique)
#             
#             for u,v,data in graph_blast.edges(data=True) :
#                 print(data["poids"])
                
                
            for elt in liste_representant :
                deja_vu = False
                for clique in cliques :
                    if elt in clique and not deja_vu : 
                        deja_vu = True
                    elif elt in clique and deja_vu :
                        print("probleme")
            
            pos = spring_layout(graph_blast)
            
            nx.draw_networkx_nodes(graph_blast, pos, node_size=5, node_color="pink")
            nx.draw_networkx_labels(graph_blast, pos, labels = node_labels , fontsize=2)
            nx.draw_networkx_edges(graph_blast,pos)
            plt.plot()
            plt.show()


## Version Needleman-Wunsch (version finale) ##

''' 18/11/19 
effectue l'alignement global par Needleman-Wunsch de toutes les paires d'occurrences d'un type d'ARN passe en parametre 
et stocke les sequences et les resultats d'alignement dans des fichiers texte'''
def alignement_global(num_ARN):
    if isinstance(num_ARN, list) :
        liste_representant = []
        for elt in num_ARN :
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr = mon_depickler_1.load()
                
                liste_representant.extend(liste_repr)
    else :
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
            mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
            liste_representant = mon_depickler_1.load()
        
    os.makedirs("Nouvelles_donnees/alignements_%s"%num_ARN, exist_ok = True)
        
    print(len(liste_representant))
    for elt in liste_representant :
            if "fichier_%s_%s_seq1.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s"%num_ARN) \
            or "fichier_%s_%s_seq2.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s"%num_ARN) \
            or "fichier_%s_%s_seq3.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s"%num_ARN) :
                with open("Nouvelles_donnees/fichier_%s_%s.pickle"%(elt[0], elt[1]), "rb") as fichier_extension : 
                    mon_depickler = pickle.Unpickler(fichier_extension)
                    extension = mon_depickler.load()
                    
                    seq1,seq2,seq3 = recup_sequences_autour_motif(elt[0], extension, 30)
                    
                    with open("Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq1.fa"%(num_ARN, elt[0], elt[1]), 'w') as fichier_seq1 :
                        fichier_seq1.write(">%s_%s_seq1\n"%(elt[0], elt[1]))
                        fichier_seq1.write(seq1)
                        
                    with open("Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq2.fa"%(num_ARN, elt[0], elt[1]), 'w') as fichier_seq2 :
                        fichier_seq2.write(">%s_%s_seq2\n"%(elt[0], elt[1]))
                        fichier_seq2.write(seq2)
                        
                    with open("Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq3.fa"%(num_ARN, elt[0], elt[1]), 'w') as fichier_seq3 :
                        fichier_seq3.write(">%s_%s_seq3\n"%(elt[0], elt[1]))
                        fichier_seq3.write(seq3)
            
    for i in range(len(liste_representant)) :
            for j in range(i+1, len(liste_representant)) :
                if "needle_%s_%s_%s_%s_seq1.txt"%(liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s/"%(num_ARN)) :               
                    needle_cline_1 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq1.fa"%(num_ARN, liste_representant[i][0],liste_representant[i][1]), bsequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq1.fa"%(num_ARN, liste_representant[j][0],liste_representant[j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq1.txt"%(num_ARN,liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]))                
                    print(needle_cline_1)
                    stdout_1, stderr_1 = needle_cline_1()
                    print(stdout_1)
                
                if  "needle_%s_%s_%s_%s_seq2.txt"%(liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s/"%(num_ARN)) :              

                    needle_cline_2 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq2.fa"%(num_ARN, liste_representant[i][0],liste_representant[i][1]), bsequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq2.fa"%(num_ARN, liste_representant[j][0],liste_representant[j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq2.txt"%(num_ARN,liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]))                
                    print(needle_cline_2)
                    stdout_2, stderr_2 = needle_cline_2()
                    print(stdout_2)
                
                if  "needle_%s_%s_%s_%s_seq3.txt"%(liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s/"%(num_ARN)) :              

                    needle_cline_3 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq3.fa"%(num_ARN, liste_representant[i][0],liste_representant[i][1]), bsequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq3.fa"%(num_ARN, liste_representant[j][0],liste_representant[j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq3.txt"%(num_ARN,liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]))                
                    print(needle_cline_3)
                    stdout_3, stderr_3 = needle_cline_3()
                    print(stdout_3)
                                
                
#                 
#                 align = AlignIO.read("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq1.txt"%(num_ARN,liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]), "emboss")
#                 print(align)
                
                
                
                #break
            #break

''' 16/01/20 
effectue l'alignement global par Needleman-Wunsch de toutes les paires d'occurrences entre deux types d'ARN passes en parametre 
et stocke les sequences et les resultats d'alignement dans des fichiers texte'''
def alignement_global_deux_types(num_ARN1, num_ARN2):
    if isinstance(num_ARN1, list) :
        liste_representant_1 = []
        for elt in num_ARN1 :
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr_1 = mon_depickler_1.load()
                
                liste_representant_1.extend(liste_repr_1)
    else :
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN1, 'rb') as fichier_num_arn :
            mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
            liste_representant_1 = mon_depickler_1.load()
    
    if isinstance(num_ARN2, list) :
        liste_representant_2 = []
        for elt in num_ARN2 :
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr_2 = mon_depickler_1.load()
                
                liste_representant_2.extend(liste_repr_2)
    else :
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN2, 'rb') as fichier_num_arn :
            mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
            liste_representant_2 = mon_depickler_1.load()  
    
    os.makedirs("Nouvelles_donnees/alignements_%s_%s"%(num_ARN1, num_ARN2), exist_ok = True)
        
    print(len(liste_representant_1))
    print(len(liste_representant_2))
    liste_tot = list(liste_representant_1)
    liste_tot.extend(liste_representant_2)
    for elt in liste_tot :
            if "fichier_%s_%s_seq1.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s_%s"%(num_ARN1, num_ARN2)) \
            or "fichier_%s_%s_seq2.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s_%s"%(num_ARN1, num_ARN2)) \
            or "fichier_%s_%s_seq3.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s_%s"%(num_ARN1, num_ARN2)) :
                with open("Nouvelles_donnees/fichier_%s_%s.pickle"%(elt[0], elt[1]), "rb") as fichier_extension : 
                    mon_depickler = pickle.Unpickler(fichier_extension)
                    extension = mon_depickler.load()
                    
                    seq1,seq2,seq3 = recup_sequences_autour_motif(elt[0], extension, 30)
                    
                    with open("Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq1.fa"%(num_ARN1, num_ARN2, elt[0], elt[1]), 'w') as fichier_seq1 :
                        fichier_seq1.write(">%s_%s_seq1\n"%(elt[0], elt[1]))
                        fichier_seq1.write(seq1)
                        
                    with open("Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq2.fa"%(num_ARN1, num_ARN2, elt[0], elt[1]), 'w') as fichier_seq2 :
                        fichier_seq2.write(">%s_%s_seq2\n"%(elt[0], elt[1]))
                        fichier_seq2.write(seq2)
                        
                    with open("Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq3.fa"%(num_ARN1, num_ARN2, elt[0], elt[1]), 'w') as fichier_seq3 :
                        fichier_seq3.write(">%s_%s_seq3\n"%(elt[0], elt[1]))
                        fichier_seq3.write(seq3)
            
    for i in range(len(liste_representant_1)) :
            for j in range(len(liste_representant_2)) :
                if "needle_%s_%s_%s_%s_seq1.txt"%(liste_representant_1[i][0],liste_representant_1[i][1], liste_representant_2[j][0],liste_representant_2[j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s_%s/"%(num_ARN1, num_ARN2)) :               
                    needle_cline_1 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq1.fa"%(num_ARN1, num_ARN2, liste_representant_1[i][0],liste_representant_1[i][1]), bsequence="Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq1.fa"%(num_ARN1, num_ARN2, liste_representant_2[j][0],liste_representant_2[j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s_%s/needle_%s_%s_%s_%s_seq1.txt"%(num_ARN1, num_ARN2,liste_representant_1[i][0],liste_representant_1[i][1], liste_representant_2[j][0],liste_representant_2[j][1]))                
                    print(needle_cline_1)
                    stdout_1, stderr_1 = needle_cline_1()
                    print(stdout_1)
                
                if  "needle_%s_%s_%s_%s_seq2.txt"%(liste_representant_1[i][0],liste_representant_1[i][1], liste_representant_2[j][0],liste_representant_2[j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s_%s/"%(num_ARN1, num_ARN2)) :              

                    needle_cline_2 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq2.fa"%(num_ARN1, num_ARN2, liste_representant_1[i][0],liste_representant_1[i][1]), bsequence="Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq2.fa"%(num_ARN1, num_ARN2, liste_representant_2[j][0],liste_representant_2[j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s_%s/needle_%s_%s_%s_%s_seq2.txt"%(num_ARN1, num_ARN2,liste_representant_1[i][0],liste_representant_1[i][1], liste_representant_2[j][0],liste_representant_2[j][1]))                
                    print(needle_cline_2)
                    stdout_2, stderr_2 = needle_cline_2()
                    print(stdout_2)
                
                if  "needle_%s_%s_%s_%s_seq3.txt"%(liste_representant_1[i][0],liste_representant_1[i][1], liste_representant_2[j][0],liste_representant_2[j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s_%s/"%(num_ARN1, num_ARN2)) :              

                    needle_cline_3 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq3.fa"%(num_ARN1, num_ARN2, liste_representant_1[i][0],liste_representant_1[i][1]), bsequence="Nouvelles_donnees/alignements_%s_%s/fichier_%s_%s_seq3.fa"%(num_ARN1, num_ARN2, liste_representant_2[j][0],liste_representant_2[j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s_%s/needle_%s_%s_%s_%s_seq3.txt"%(num_ARN1, num_ARN2,liste_representant_1[i][0],liste_representant_1[i][1], liste_representant_2[j][0],liste_representant_2[j][1]))                
                    print(needle_cline_3)
                    stdout_3, stderr_3 = needle_cline_3()
                    print(stdout_3)
    
                
#                 
#                 align = AlignIO.read("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq1.txt"%(num_ARN,liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1]), "emboss")
#                 print(align)
                
                
                
                #break
            #break



''' 16/01/20 
effectue l'alignement global par Needleman-Wunsch de toutes les paires d'occurrences entre plusieurs types d'ARN passes en parametre 
et stocke les sequences et les resultats d'alignement dans des fichiers texte'''
def alignement_global_plusieurs_types(liste_num_ARN):
    liste_representant_1 = []
    compteur = 0
    for num_ARN in liste_num_ARN :
        if len(liste_representant_1) < compteur + 1 :
            liste_representant_1.append([])
        if isinstance(num_ARN, list) :
            for elt in num_ARN :
                with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_repr_1 = mon_depickler_1.load()
                    
                    liste_representant_1[compteur].extend(liste_repr_1)
        else :
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant_1[compteur] = mon_depickler_1.load()
        compteur += 1
    
    
    
    os.makedirs("Nouvelles_donnees/alignements_%s"%(liste_num_ARN), exist_ok = True)
        
    print(len(liste_representant_1))
    liste_tot = []
    for elt in liste_representant_1 :
        liste_tot.extend(elt)
    print(liste_tot)
    for elt in liste_tot :
            if "fichier_%s_%s_seq1.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s"%(liste_num_ARN)) \
            or "fichier_%s_%s_seq2.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s"%(liste_num_ARN)) \
            or "fichier_%s_%s_seq3.fa"%(elt[0],elt[1]) not in os.listdir("Nouvelles_donnees/alignements_%s"%(liste_num_ARN)) :
                with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(elt[0], elt[1]), "rb") as fichier_extension : 
                    mon_depickler = pickle.Unpickler(fichier_extension)
                    extension = mon_depickler.load()
                    
                    seq1,seq2,seq3 = recup_sequences_autour_motif(elt[0], extension, 30)
                    
                    with open("Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq1.fa"%(liste_num_ARN, elt[0], elt[1]), 'w') as fichier_seq1 :
                        fichier_seq1.write(">%s_%s_seq1\n"%(elt[0], elt[1]))
                        fichier_seq1.write(seq1)
                        
                    with open("Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq2.fa"%(liste_num_ARN, elt[0], elt[1]), 'w') as fichier_seq2 :
                        fichier_seq2.write(">%s_%s_seq2\n"%(elt[0], elt[1]))
                        fichier_seq2.write(seq2)
                        
                    with open("Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq3.fa"%(liste_num_ARN, elt[0], elt[1]), 'w') as fichier_seq3 :
                        fichier_seq3.write(">%s_%s_seq3\n"%(elt[0], elt[1]))
                        fichier_seq3.write(seq3)
    for k in range(len(liste_representant_1)) :
        for l in range(k+1, len(liste_representant_1)) :      
            for i in range(len(liste_representant_1[k])) :
                for j in range(len(liste_representant_1[l])) :
                    if "needle_%s_%s_%s_%s_seq1.txt"%(liste_representant_1[k][i][0],liste_representant_1[k][i][1], liste_representant_1[l][j][0],liste_representant_1[l][j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s/"%(liste_num_ARN)) :               
                        needle_cline_1 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq1.fa"%(liste_num_ARN, liste_representant_1[k][i][0],liste_representant_1[k][i][1]), bsequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq1.fa"%(liste_num_ARN, liste_representant_1[l][j][0],liste_representant_1[l][j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq1.txt"%(liste_num_ARN,liste_representant_1[k][i][0],liste_representant_1[k][i][1], liste_representant_1[l][j][0],liste_representant_1[l][j][1]))                
                        print(needle_cline_1)
                        stdout_1, stderr_1 = needle_cline_1()
                        print(stdout_1)
                    
                    if  "needle_%s_%s_%s_%s_seq2.txt"%(liste_representant_1[k][i][0],liste_representant_1[k][i][1], liste_representant_1[l][j][0],liste_representant_1[l][j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s/"%(liste_num_ARN)) :              
    
                        needle_cline_2 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq2.fa"%(liste_num_ARN, liste_representant_1[k][i][0],liste_representant_1[k][i][1]), bsequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq2.fa"%(liste_num_ARN, liste_representant_1[l][j][0],liste_representant_1[l][j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq2.txt"%(liste_num_ARN,liste_representant_1[k][i][0],liste_representant_1[k][i][1], liste_representant_1[l][j][0],liste_representant_1[l][j][1]))                
                        print(needle_cline_2)
                        stdout_2, stderr_2 = needle_cline_2()
                        print(stdout_2)
                    
                    if  "needle_%s_%s_%s_%s_seq3.txt"%(liste_representant_1[k][i][0],liste_representant_1[k][i][1], liste_representant_1[l][j][0],liste_representant_1[l][j][1]) not in os.listdir("Nouvelles_donnees/alignements_%s/"%(liste_num_ARN)) :              
    
                        needle_cline_3 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq3.fa"%(liste_num_ARN, liste_representant_1[k][i][0],liste_representant_1[k][i][1]), bsequence="Nouvelles_donnees/alignements_%s/fichier_%s_%s_seq3.fa"%(liste_num_ARN, liste_representant_1[l][j][0],liste_representant_1[l][j][1]), gapopen=10, gapextend=0.5, outfile="Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq3.txt"%(liste_num_ARN,liste_representant_1[k][i][0],liste_representant_1[k][i][1], liste_representant_1[l][j][0],liste_representant_1[l][j][1]))                
                        print(needle_cline_3)
                        stdout_3, stderr_3 = needle_cline_3()
                        print(stdout_3)




'''19/11/19
recupere le score d'un alignement global passe en parametre'''
def recup_score(fichier_alignement_global):
    with open(fichier_alignement_global, 'r') as fichier_aln_global :
        ligne = fichier_aln_global.readline()
        while ligne != "" :
            if "Score" in ligne :
                score = float(ligne.split(" ")[2])
                return score
            ligne = fichier_aln_global.readline()


'''19/11/19
construit une matrice de distance a partir des valeurs de score d'alignement stockees en ponderation des aretes du graphe passe en parametre 
dans le but de lancer kmeans'''
def matrice_distance_selon_score_aln(graphe):
    matrice = [[0] *graphe.number_of_nodes() for _ in range(graphe.number_of_nodes())]
    
    for u,v,data in graphe.edges(data=True) :
        matrice[u][v] = 310 - data["aln"]
        matrice[v][u] = 310 - data["aln"]
        
    return matrice


'''19/11/19
traitement d'un alignement global
- distribution des valeurs de score
- distribution des valeurs de RMSD en fonction des valeurs de score avec couleurs differentes pour
paires en positions similaires : jaune
paires homologues verifiees a la main : vert
paires non homologues verifiees a la main : rouge
- creation de fichiers pickle de stockage des valeurs de rmsd, score et score par positions similaires dans des dictionnaires (cle : (paire1, paire2) exemple : (('5dm6', 3), ('4ybb', 21)))
num_ARN : type d'ARN
seuil_hom : seuil de nombres de nts pour les positions similaires
min_ou_max : si on considere le score minimum ou maximum sur les 3 sequences (mais min, c'est mieux, indéniablement)
liste_aretes : une liste d'aretes a mettre d'une certaine couleur
dico_sim, seuil_sim : pas utilise'''
def traitement_aln_global(num_ARN, dico_sim, seuil_hom, seuil_sim, min_ou_max, liste_aretes):
    liste_cles = {}
    if isinstance(num_ARN, list) :
        
        liste_representant = []
        for elt in num_ARN :
            liste_cles.update({elt : []})
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr = mon_depickler_1.load()
                
                liste_representant.extend(liste_repr)
                liste_cles[elt].extend(liste_repr)
        
    else :
        liste_cles.update({num_ARN : []})
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant = mon_depickler_1.load()
                liste_cles[num_ARN].extend(liste_representant)
                
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
                
    with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv :
        csvreader = csv.reader(fichier_csv) 
        
    with open("fichier_type_seq_%s_seuil_%s.csv"%(num_ARN, seuil_hom), 'w', newline='') as fichier_csv :
        csvwriter = csv.writer(fichier_csv) 
        csvwriter.writerow(["seuil pos", "nombre seq1", "nombre seq2", "nombre seq3"])
    
#     with open("fichier_orange_%s.csv"%(num_ARN), 'w', newline='') as fichier_csv_2 :
#         csvwriter2 = csv.writer(fichier_csv_2) 
#         csvwriter2.writerow(["Num paires","Score min","RMSD","Diff des positions"])
#         
    with open("fichier_bleu_frontiere_4_%s.csv"%(num_ARN), 'w', newline='') as fichier_csv_3 :
            csvwriter3 = csv.writer(fichier_csv_3) 
            csvwriter3.writerow(["Num paires","Score min","RMSD","Diff des positions"])
        
#         with open("dico_min_pos_score_%s_seuil_%s.pickle"%(num_ARN, seuil_hom), "rb") as fichier_score_pos :
#             mon_depickler = pickle.Unpickler(fichier_score_pos)
#             dico_min_pos_scores = mon_depickler.load()
           
            liste_verte = []
            # 23S
            #liste_verte = [(('1vqo', 12), ('6nd6', 7)),(('4y4o', 40), ('4u4r', 8)),(('6eri', 9), ('4u4r', 8)), (('6hma', 15), ('1vq8', 15)), (('5dm6', 2), ('6eri', 1)), (('4u27', 5), ('1vqp', 19)),(('4u27', 5), ('3ccr', 4)),(('4u27', 5), ('3cc7', 4)),(('5dm6', 7), ('1vq8', 11)),(('5dm6', 7), ('3ccr', 12)),(('5dm6', 7), ('3ccl', 11)),(('5dm6', 7), ('3ccs', 12)),(('5dm6', 7), ('3cc7', 12)),(('5dm6', 7), ('3ccq', 12)),(('5dm6', 7), ('3cce', 9)),(('5dm6', 7), ('3ccu', 11)),(('4ybb', 8), ('1vq8', 11)),(('4ybb', 8), ('3ccr', 12)),(('4ybb', 8), ('3ccl', 11)),(('4ybb', 8), ('3ccs', 12)),(('4ybb', 8), ('3cc7', 12)),(('4ybb', 8), ('3ccq', 12)),(('4ybb', 8), ('3cce', 9)),(('4ybb', 8), ('3ccu', 11)),(('3cc2', 2), ('6hma', 1)),(('5dm6', 13), ('1vqp', 19)),(('5dm6', 13), ('3ccr', 4)),(('5dm6', 13), ('3cc7', 4)), (('6hma', 9), ('4u4r', 8)), (('4y4o', 25), ('1vq8', 5)),(('4y4o', 25), ('3cc2', 2)),(('6hma', 15), ('4u4r', 16)), (('1vqp', 19), ('6eri', 9)), (('3ccr', 4), ('6eri', 9)),(('6eri', 9), ('3cc7', 4)),(('5afi', 15), ('4u4r', 28)),(('6hma', 9), ('1vqp', 19)),(('6hma', 9), ('3ccr', 4)),(('6hma', 9), ('3cc7', 4)),(('1vq8', 11), ('6hma', 6)),(('3ccr', 12), ('6hma', 6)),(('6hma', 6), ('3ccl', 11)),(('6hma', 6), ('3ccs', 12)),(('6hma', 6), ('3cc7', 12)),(('6hma', 6), ('3ccq', 12)),(('6hma', 6), ('3cce', 9)),(('6hma', 6), ('3ccu', 11)),(('4y4o', 40), ('1vqp', 19)),(('4y4o', 40), ('3ccr', 4)),(('4y4o', 40), ('3cc7', 4)),(('4y4o', 3), ('3ccl', 11)),(('4y4o', 3), ('3cce', 9)),(('4v51', 22), ('3ccl', 11)),(('4v51', 22), ('3cce', 9)),(('5e81', 12), ('3ccl', 11)),(('5e81', 12), ('3cce', 9)),(('3ccl', 11), ('4v90', 13)),(('4v90', 13), ('3cce', 9)),(('1vqo', 12), ('4u4r', 15)),(('4ybb', 1), ('1vq8', 5)),(('4y4o', 3), ('3ccr', 12)),(('4y4o', 3), ('3ccs', 12)),(('4y4o', 3), ('3cc7', 12)),(('4y4o', 3), ('3ccu', 11)),(('4v51', 22), ('3ccr', 12)),(('4v51', 22), ('3ccs', 12)),(('4v51', 22), ('3cc7', 12)),(('4v51', 22), ('3ccu', 11)),(('5e81', 12), ('3ccr', 12)),(('5e81', 12), ('3ccs', 12)),(('5e81', 12), ('3cc7', 12)),(('5e81', 12), ('3ccu', 11)),(('3ccr', 12), ('4v90', 13)),(('4v90', 13), ('3ccs', 12)),(('4v90', 13), ('3cc7', 12)),(('4v90', 13), ('3ccu', 11)),(('4y4o', 3), ('3ccq', 12)),(('4v51', 22), ('3ccq', 12)),(('5e81', 12), ('3ccq', 12)),(('4v90', 13), ('3ccq', 12)),(('4ybb', 1), ('3cc2', 2)),(('4y4o', 3), ('1vq8', 11)),(('4v51', 22), ('1vq8', 11)),(('5e81', 12), ('1vq8', 11)),(('1vq8', 5), ('6hma', 1)),(('1vq8', 11), ('4v90', 13)), (('5dm6', 9), ('1vq8', 21)), (('6hma', 2), ('4u4r', 19)), (('4ybb', 6), ('1yhq', 24)), (('6hma', 15), ('3cc2', 6))]
            # 16S
            #liste_verte.extend([(('5nwy', 17), ('4y4o', 4)), (('6ek0', 2), ('6az1', 1)), (('4u27', 54), ('4y4o', 22)), (('4ybb', 35), ('4u3u', 22))])
            ## homologues cherches a la main avec bonne sim et rmsd
            #liste_verte.extend([(('4y4o', 13), ('4u3u', 22)),(('3cc2', 1), ('4u4r', 36)),(('1vqo', 14), ('4u4r', 36)),(('6eri', 3), ('4u3u', 2)),(('4ybb', 6), ('4u3u', 2)),(('4y4o', 9), ('4u4r', 11)),(('4v67', 23), ('4u4r', 11)),(('4v67', 23), ('6ek0', 6)),(('4y4o', 9), ('6ek0', 6)),(('4y4o', 8), ('3cc2', 1)),(('3t1y', 8), ('4u4r', 11)),(('3t1y', 8), ('6ek0', 6)),(('2qex', 19), ('4u3u', 2)),(('1yhq', 24), ('4u3u', 2)),(('4y4o', 13), ('6az1', 6)),(('4y4o', 8), ('1vqo', 14)),(('5ngm', 1), ('4u4r', 11)),(('6eri', 15), ('4u4r', 11)),(('5ngm', 1), ('6ek0', 6)),(('6eri', 15), ('6ek0', 6)),(('4ybb', 9), ('4u4r', 11)),(('5j7l', 10), ('4u4r', 11)),(('4u27', 32), ('4u4r', 11)),(('5j7l', 10), ('6ek0', 6)),(('4ybb', 9), ('6ek0', 6)),(('4u27', 32), ('6ek0', 6)),(('1vq8', 19), ('4u3u', 2)),(('4ybb', 35), ('6az1', 6)),(('2vqe', 5), ('4u4r', 11)),(('2vqe', 5), ('6ek0', 6))])
            
            ## homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim  a 0.75 et pas avec la rmsd a 2
            #liste_verte.extend([(('4y4o', 25),('1vq8', 5)),(('4y4o', 25),('3cc2', 2)),(('4ybb', 1),('1vq8', 5)),(('4ybb', 1),('3cc2', 2)), (('6hma', 1),('1vq8', 5)),(('6hma', 1),('3cc2', 2))])
            
            liste_rouge = []
            # 23S
            #liste_rouge = [(('1vq8', 18), ('5dm6', 5)),(('5dm6', 5), ('3ccl', 14)),(('5dm6', 5), ('3cce', 12)),(('1vq8', 1), ('6eri', 2)), (('5wfs', 19), ('4u4r', 6)),(('4ybb', 12), ('4y4o', 35)),(('4ybb', 12), ('4wsd', 40)),(('4y4o', 34), ('2qex', 12)),(('4y4o', 28), ('5wfs', 6)),(('4v67', 43), ('2qex', 12)), (('4y4o', 50), ('6ek0', 4)), (('5dm6', 4), ('4ybb', 17)), (('4w2g', 52), ('6ek0', 10)), (('4y4o', 47), ('1vq8', 19)),(('4y4o', 47), ('1yhq', 24)),(('4y4o', 47), ('2qex', 19)),(('4y4o', 47), ('4u4r', 13)),(('4ybb', 13), ('1vqo', 8)),(('4w2f', 23), ('4u4r', 36)),(('1vqo', 8), ('5dm6', 1)),(('1vqo', 8), ('1mms', 1)),(('1vq8', 19), ('6eri', 1)),(('1yhq', 24), ('6eri', 1)),(('6eri', 1), ('2qex', 19)), (('5dm6', 14), ('4y4o', 8)),(('4ybb', 17), ('4y4o', 8)), (('4ybb', 30), ('2zjr', 15)),(('2zjr', 15), ('5wfs', 11)), (('2zjr', 15), ('6eri', 8)), (('4y4o', 12), ('2zjr', 15)), (('4w2g', 52), ('4y4o', 28)), (('4y4o', 28), ('4v67', 7))]
            #liste_rouge.extend([(('4ybb', 12), ('4y4o', 38)),(('4ybb', 12), ('4ybb', 21)),(('4ybb', 12), ('6ek0', 7)),(('4v67', 7), ('4y4o', 38)),(('4v67', 7), ('4ybb', 21)),(('4v67', 7), ('6ek0', 7)),(('5wfs', 19), ('4y4o', 38)),(('5wfs', 19), ('4ybb', 21)),(('5wfs', 19), ('6ek0', 7))])
            ## non homologues meme type rmsd 2-2.5, sim 0.7-0.75
            #liste_rouge.extend([(('4y4o', 23), ('6h4n', 24)),(('2zjr', 3), ('4u4r', 18)),(('4u27', 54), ('4y4o', 38)),(('4u27', 54), ('6ek0', 7)),(('4u27', 54), ('4u3u', 7)), (('4y4o', 23), ('4u27', 3)),(('4y4o', 23), ('5afi', 24))])
            
            ## les non homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim à 0.75 mais pas avec la rmsd à 2
            #liste_rouge.extend([(('4ybb', 12), ('6ek0', 10)),(('4w2g', 52), ('6ek0', 10)),(('2zjr', 3),('6ek0', 10)),(('2zjr', 3),('4u4r', 18)),(('4w2g', 52),('4u4r', 18)),(('4v67', 7), ('4u4r', 18)),(('5wfs', 19), ('4u4r', 18)),(('4v67', 7),('6ek0', 10)),(('5wfs', 19),('6ek0', 10)),(('4y4o', 23),('6h4n', 24)),(('4y4o', 23), ('4u27', 3)),(('4y4o', 23),('5afi', 24))])
            
            liste_clique = []
            #liste_clique_23S_seuil_80_2.5 = [(('5dm6', 10), ('1vq8', 4)), (('5dm6', 10), ('3cc2', 23)), (('5dm6', 10), ('4u4r', 19)), (('5dm6', 10), ('6ek0', 12)), (('5dm6', 9), ('1vq8', 21)), (('5dm6', 9), ('3cc2', 17)), (('5dm6', 9), ('4u4r', 4)), (('5dm6', 9), ('6ek0', 4)), (('6hrm', 22), ('3cma', 14)), (('4ybb', 18), ('1vq8', 4)), (('4ybb', 18), ('3cc2', 23)), (('4ybb', 18), ('4u4r', 19)), (('4ybb', 18), ('6ek0', 12)), (('4ybb', 30), ('4y4o', 12)), (('4ybb', 54), ('1vq8', 21)), (('4ybb', 54), ('3cc2', 17)), (('4ybb', 54), ('6ek0', 4)), (('4y4o', 12), ('6hma', 15)), (('4y4o', 15), ('6eri', 23)), (('4y4o', 3), ('3ccr', 12)), (('4y4o', 3), ('3ccl', 11)), (('4y4o', 3), ('3ccs', 12)), (('4y4o', 3), ('3ccq', 12)), (('4y4o', 28), ('6ek0', 4)), (('4y4o', 2), ('1vqo', 24)), (('4y4o', 2), ('1yij', 19)), (('4y4o', 53), ('4u4r', 19)), (('4y4o', 53), ('6ek0', 12)), (('4ybb', 5), ('1vqo', 24)), (('4ybb', 5), ('1yij', 19)), (('6hma', 2), ('4u4r', 19)), (('6hma', 2), ('6ek0', 12)), (('6hma', 8), ('1vq8', 21)), (('6hma', 8), ('3cc2', 17)), (('6hma', 8), ('6ek0', 4)), (('6hma', 5), ('1vqo', 24)), (('6hma', 5), ('1yij', 19)), (('4v51', 22), ('1vq8', 11)), (('4v51', 22), ('3ccr', 12)), (('4v51', 22), ('3ccl', 11)), (('4v51', 22), ('3ccs', 12)), (('4v51', 22), ('3cc7', 12)), (('4v51', 22), ('3ccq', 12)), (('4v51', 22), ('3cce', 9)), (('4v51', 22), ('3ccu', 11)), (('4v51', 49), ('4u4r', 19)), (('4v51', 49), ('6ek0', 12)), (('2zjr', 3), ('5wfs', 19)), (('5e81', 12), ('1vq8', 11)), (('5e81', 12), ('3ccr', 12)), (('5e81', 12), ('3ccl', 11)), (('5e81', 12), ('3ccs', 12)), (('5e81', 12), ('3cc7', 12)), (('5e81', 12), ('3ccq', 12)), (('5e81', 12), ('3cce', 9)), (('1vq8', 4), ('6ek0', 12)), (('1vq8', 21), ('5afi', 17)), (('1vq8', 21), ('6ek0', 4)), (('1vq8', 22), ('5dm6', 11)), (('1vq8', 22), ('5wfs', 6)), (('1vq8', 22), ('6eri', 24)), (('5dm6', 11), ('3cc2', 12)), (('3cc2', 12), ('5wfs', 6)), (('3cc2', 17), ('5afi', 17)), (('3cc2', 17), ('6eri', 11)), (('3cc2', 17), ('4u4r', 4)), (('3cc2', 17), ('6ek0', 4)), (('3cc2', 23), ('6ek0', 12)), (('5afi', 17), ('6ek0', 4)), (('4v51', 5), ('6eri', 23)), (('3ccr', 12), ('4v90', 13)), (('6hma', 6), ('3ccs', 12)), (('3ccl', 11), ('4v90', 13)), (('6eri', 11), ('4u4r', 4)), (('6eri', 11), ('6ek0', 4)), (('4v90', 13), ('3ccs', 12)), (('5wdt', 12), ('4u4r', 19)), (('5wdt', 12), ('6ek0', 12)), (('6eri', 16), ('6ek0', 3))]
            #liste_clique_16S_seuil_80_2.5 = [(('4ybb', 46), ('5ngm', 15)), (('4ybb', 51), ('5ngm', 9)), (('4y4o', 31), ('4u4r', 10)), (('5ndk', 13), ('4u4r', 10)), (('4v8d', 26), ('4u4r', 10)), (('4v7l', 15), ('4u4r', 10)), (('5f8k', 24), ('4u4r', 10)), (('3t1y', 10), ('4u4r', 10)), (('4v67', 6), ('4u4r', 10)), (('4v90', 11), ('4u4r', 10)), (('5e81', 10), ('4u4r', 10)), (('6eri', 18), ('6az1', 2))]
            #liste_clique_23S_seuil_60_2 = [(('5dm6', 7), ('1vq8', 11)), (('5dm6', 7), ('3ccr', 12)), (('5dm6', 7), ('3ccl', 11)), (('5dm6', 7), ('3ccs', 12)), (('5dm6', 7), ('3cc7', 12)), (('5dm6', 7), ('3ccq', 12)), (('5dm6', 7), ('3cce', 9)), (('5dm6', 7), ('3ccu', 11)), (('5dm6', 10), ('6ek0', 12)), (('4ybb', 6), ('1vq8', 19)), (('4ybb', 6), ('2qex', 19)), (('4ybb', 14), ('1vq8', 13)), (('4ybb', 14), ('2qex', 8)), (('4ybb', 8), ('1vq8', 11)), (('4ybb', 8), ('3ccr', 12)), (('4ybb', 8), ('3ccl', 11)), (('4ybb', 8), ('3ccs', 12)), (('4ybb', 8), ('3cc7', 12)), (('4ybb', 8), ('3ccq', 12)), (('4ybb', 8), ('3cce', 9)), (('4ybb', 8), ('3ccu', 11)), (('4ybb', 18), ('4u4r', 19)), (('4ybb', 18), ('6ek0', 12)), (('4ybb', 30), ('1vq8', 15)), (('4ybb', 30), ('3cc2', 6)), (('4ybb', 30), ('4u4r', 16)), (('4y4o', 34), ('6hma', 13)), (('4y4o', 12), ('1vq8', 15)), (('4y4o', 12), ('3cc2', 6)), (('4y4o', 15), ('2qex', 8)), (('4y4o', 15), ('5tbw', 7)), (('4y4o', 3), ('1vq8', 11)), (('4y4o', 3), ('3ccr', 12)), (('4y4o', 3), ('3ccl', 11)), (('4y4o', 3), ('3ccs', 12)), (('4y4o', 3), ('3cc7', 12)), (('4y4o', 3), ('3ccq', 12)), (('4y4o', 3), ('3cce', 9)), (('4y4o', 3), ('3ccu', 11)), (('4y4o', 53), ('4u4r', 19)), (('4y4o', 53), ('6ek0', 12)), (('4ybb', 5), ('1vqo', 24)), (('4ybb', 5), ('1yij', 19)), (('5mdv', 21), ('1vq8', 13)), (('5mdv', 21), ('2qex', 8)), (('6hma', 15), ('1vq8', 15)), (('6hma', 15), ('4u4r', 16)), (('6hma', 2), ('6ek0', 12)), (('6hma', 11), ('4u4r', 12)), (('4v51', 22), ('1vq8', 11)), (('4v51', 22), ('3ccr', 12)), (('4v51', 22), ('3ccl', 11)), (('4v51', 22), ('3ccs', 12)), (('4v51', 22), ('3cc7', 12)), (('4v51', 22), ('3ccq', 12)), (('4v51', 22), ('3cce', 9)), (('4v51', 22), ('3ccu', 11)), (('4v51', 49), ('4u4r', 19)), (('4v51', 49), ('6ek0', 12)), (('5e81', 12), ('1vq8', 11)), (('5e81', 12), ('3ccr', 12)), (('5e81', 12), ('3ccl', 11)), (('5e81', 12), ('3ccs', 12)), (('5e81', 12), ('3cc7', 12)), (('5e81', 12), ('3ccq', 12)), (('5e81', 12), ('3cce', 9)), (('5e81', 12), ('3ccu', 11)), (('1vq8', 11), ('6hma', 6)), (('1vq8', 11), ('4v90', 13)), (('1vq8', 13), ('6qul', 18)), (('1vq8', 13), ('5tbw', 7)), (('1vq8', 15), ('6eri', 8)), (('1vq8', 15), ('5wfs', 11)), (('1vq8', 15), ('4u4r', 16)), (('1vq8', 19), ('6eri', 3)), (('5dm6', 5), ('6hma', 13)), (('5dm6', 1), ('5d8h', 1)), (('5dm6', 1), ('1mms', 1)), (('5dm6', 1), ('6eri', 4)), (('3cc2', 6), ('6eri', 8)), (('3cc2', 6), ('5wfs', 11)), (('3cc2', 6), ('4u4r', 16)), (('3cc2', 14), ('5tbw', 7)), (('1yhq', 24), ('6eri', 3)), (('5d8h', 1), ('6eri', 4)), (('4v51', 5), ('5tbw', 7)), (('4v9f', 8), ('6qul', 18)), (('4v9f', 8), ('5tbw', 7)), (('3ccr', 12), ('6hma', 6)), (('3ccr', 12), ('4v90', 13)), (('6hma', 6), ('3ccl', 11)), (('6hma', 6), ('3ccs', 12)), (('6hma', 6), ('3cc7', 12)), (('6hma', 6), ('3ccq', 12)), (('6hma', 6), ('3cce', 9)), (('6hma', 6), ('3ccu', 11)), (('3ccl', 11), ('4v90', 13)), (('6eri', 8), ('4u4r', 16)), (('4v90', 13), ('3ccs', 12)), (('4v90', 13), ('3cc7', 12)), (('4v90', 13), ('3ccq', 12)), (('4v90', 13), ('3cce', 9)), (('4v90', 13), ('3ccu', 11)), (('5wdt', 12), ('4u4r', 19)), (('5wdt', 12), ('6ek0', 12)), (('6eri', 3), ('2qex', 19)), (('6qul', 18), ('2qex', 8)), (('6qul', 18), ('5tbw', 7)), (('2qex', 8), ('5tbw', 7)), (('5afi', 15), ('4u4r', 28))]
            #liste_clique_16S_seuil_60_2  = [(('4ybb', 46), ('5ngm', 15)), (('4ybb', 11), ('4u4r', 24)), (('4y4o', 38), ('6ek0', 7)), (('4y4o', 38), ('4u3u', 7)), (('4ybb', 21), ('6ek0', 7)), (('4y4o', 54), ('4u4r', 7))]
            # liste clique seuil 0 2
            #liste_clique = [(('5dm6', 7), ('1vq8', 11)), (('5dm6', 7), ('3ccr', 12)), (('5dm6', 7), ('3ccl', 11)), (('5dm6', 7), ('3ccs', 12)), (('5dm6', 7), ('3cc7', 12)), (('5dm6', 7), ('3ccq', 12)), (('5dm6', 7), ('3cce', 9)), (('5dm6', 7), ('3ccu', 11)), (('4ybb', 8), ('1vq8', 11)), (('4ybb', 8), ('3ccr', 12)), (('4ybb', 8), ('3ccl', 11)), (('4ybb', 8), ('3ccs', 12)), (('4ybb', 8), ('3cc7', 12)), (('4ybb', 8), ('3ccq', 12)), (('4ybb', 8), ('3cce', 9)), (('4ybb', 8), ('3ccu', 11)), (('4ybb', 37), ('4ybb', 1)), (('4ybb', 37), ('4y4o', 25)), (('4ybb', 37), ('6hma', 1)), (('4ybb', 1), ('4y4o', 19)), (('4ybb', 1), ('1vq8', 16)), (('4ybb', 1), ('1yhq', 17)), (('4ybb', 1), ('6h4n', 7)), (('4ybb', 1), ('3ccr', 17)), (('4ybb', 1), ('6hma', 10)), (('4ybb', 1), ('3cd6', 14)), (('4ybb', 1), ('3ccs', 16)), (('4ybb', 1), ('6eri', 16)), (('4ybb', 1), ('3cc7', 17)), (('4ybb', 1), ('6qul', 4)), (('4ybb', 1), ('3ccq', 18)), (('4ybb', 1), ('3ccm', 16)), (('4ybb', 1), ('3ccu', 16)), (('4ybb', 1), ('2qex', 20)), (('4ybb', 1), ('4u4r', 6)), (('4ybb', 1), ('6ek0', 3)), (('4y4o', 23), ('1vq8', 12)), (('4y4o', 23), ('3cc2', 19)), (('4y4o', 23), ('3cpw', 14)), (('4y4o', 34), ('6hma', 13)), (('4y4o', 3), ('1vq8', 11)), (('4y4o', 3), ('3ccr', 12)), (('4y4o', 3), ('3ccl', 11)), (('4y4o', 3), ('3ccs', 12)), (('4y4o', 3), ('3cc7', 12)), (('4y4o', 3), ('3ccq', 12)), (('4y4o', 3), ('3cce', 9)), (('4y4o', 3), ('3ccu', 11)), (('4y4o', 25), ('1vq8', 16)), (('4y4o', 25), ('1yhq', 17)), (('4y4o', 25), ('6h4n', 7)), (('4y4o', 25), ('3ccr', 17)), (('4y4o', 25), ('6hma', 10)), (('4y4o', 25), ('3cd6', 14)), (('4y4o', 25), ('3ccs', 16)), (('4y4o', 25), ('6eri', 16)), (('4y4o', 25), ('3cc7', 17)), (('4y4o', 25), ('6qul', 4)), (('4y4o', 25), ('3ccq', 18)), (('4y4o', 25), ('3ccm', 16)), (('4y4o', 25), ('3ccu', 16)), (('4y4o', 25), ('2qex', 20)), (('4y4o', 25), ('4u4r', 6)), (('4y4o', 25), ('6ek0', 3)), (('4y4o', 19), ('6hma', 1)), (('6hma', 15), ('1vq8', 15)), (('6hma', 15), ('4u4r', 16)), (('4v51', 22), ('1vq8', 11)), (('4v51', 22), ('3ccr', 12)), (('4v51', 22), ('3ccl', 11)), (('4v51', 22), ('3ccs', 12)), (('4v51', 22), ('3cc7', 12)), (('4v51', 22), ('3ccq', 12)), (('4v51', 22), ('3cce', 9)), (('4v51', 22), ('3ccu', 11)), (('5e81', 12), ('1vq8', 11)), (('5e81', 12), ('3ccr', 12)), (('5e81', 12), ('3ccl', 11)), (('5e81', 12), ('3ccs', 12)), (('5e81', 12), ('3cc7', 12)), (('5e81', 12), ('3ccq', 12)), (('5e81', 12), ('3cce', 9)), (('5e81', 12), ('3ccu', 11)), (('1vq8', 11), ('6hma', 6)), (('1vq8', 11), ('4v90', 13)), (('1vq8', 16), ('6hma', 1)), (('5dm6', 5), ('6hma', 13)), (('5dm6', 1), ('5d8h', 1)), (('5dm6', 1), ('1mms', 1)), (('5dm6', 1), ('6eri', 4)), (('1yhq', 17), ('6hma', 1)), (('6h4n', 7), ('6hma', 1)), (('5d8h', 1), ('6eri', 4)), (('3ccr', 12), ('6hma', 6)), (('3ccr', 12), ('4v90', 13)), (('3ccr', 17), ('6hma', 1)), (('6hma', 1), ('6hma', 10)), (('6hma', 1), ('3cd6', 14)), (('6hma', 1), ('3ccs', 16)), (('6hma', 1), ('3cc7', 17)), (('6hma', 1), ('6qul', 4)), (('6hma', 1), ('3ccq', 18)), (('6hma', 1), ('3ccm', 16)), (('6hma', 1), ('3ccu', 16)), (('6hma', 1), ('2qex', 20)), (('6hma', 1), ('4u4r', 6)), (('6hma', 1), ('6ek0', 3)), (('6hma', 6), ('3ccl', 11)), (('6hma', 6), ('3ccs', 12)), (('6hma', 6), ('3cc7', 12)), (('6hma', 6), ('3ccq', 12)), (('6hma', 6), ('3cce', 9)), (('6hma', 6), ('3ccu', 11)), (('3ccl', 11), ('4v90', 13)), (('4v90', 13), ('3ccs', 12)), (('4v90', 13), ('3cc7', 12)), (('4v90', 13), ('3ccq', 12)), (('4v90', 13), ('3cce', 9)), (('4v90', 13), ('3ccu', 11))]

            dico_min_pos_scores = {}
            dico_min_tot = {}
            distrib_min_tot_zoom = []
            dico_min = {}
            distrib_min_tot = []  
            distrib_min = []
            distrib_min_pos_sim = []
            dico_min_pos_sim = {}
            dico_not_min_pos_sim = {}
            distrib_diff_pos = []
            liste_en_dessous_de_100_rmsd = []
            
            liste_seq_min = [0,0,0]
            liste_en_dessous_de_100 = []
            liste_en_dessous_de_100_sim = []
            liste_en_dessous_de_100_pos = []
            dico_min_tot_rmsd = {}
            dico_min_pos_rmsd = {}
            distrib_min_rmsd = []
            distrib_min_pos_rmsd = []
            distrib_score_en_plus = []
            distrib_rmsd_en_plus = []
            
            distrib_min_cluster = []
            distrib_min_cluster_rmsd = []
            compter = 0
            
            dico_type_seq = [0,0,0]
            dico_type_seq_pos = [0,0,0]
            dico_type_seq_bas = [0,0,0]
            dico_type_seq_haut = [0,0,0]
            
            distrib_score_rouge = []
            distrib_score_verte = []
            
            distrib_rmsd_rouge = []
            distrib_rmsd_verte = []
            
            distrib_score_clique = []
            distrib_rmsd_clique = []
            
            clique_rouge = 0
            
            #clusters = traitement_res_cdhit(num_ARN)
            for i in range(len(liste_representant)) :
                for j in range(i+1, len(liste_representant)) :
                    #if liste_representant[i][0] != '6hrm' and  liste_representant[j][0] != '6hrm' :
                        scores = []
                        for k in range(1,4) :
                            scores.append(recup_score("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq%d.txt"%(num_ARN,liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1], k)))
                        #print(scores)
                        if min_ou_max == "min" :
                            distrib_min_tot.append(min(scores))
                            dico_min_tot.update({(liste_representant[i], liste_representant[j]) : min(scores)})
                                
                        elif min_ou_max == "max" :
                            distrib_min_tot.append(max(scores))
                            dico_min_tot.update({(liste_representant[i], liste_representant[j]) : max(scores)})
                        else :
                            print("min ou max pas bon")
                            exit()
                        
                        nom1 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[i][0], liste_representant[i][1])
                        nom2 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[j][0], liste_representant[j][1])
                        
                        if (nom1, nom2) in rmsd.keys() :
#                             if (liste_representant[i], liste_representant[j]) in dico_min_pos_scores.keys() or (liste_representant[j], liste_representant[i]) in dico_min_pos_scores.keys() : 
#                                 dico_min_pos_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom1, nom2)]})
                            distrib_min_rmsd.append(rmsd[(nom1, nom2)])
                            dico_min_tot_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom1, nom2)]})
                            if (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                distrib_score_verte.append(min(scores))
                                distrib_rmsd_verte.append(rmsd[(nom1, nom2)])
                                
                            if (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge :
                                distrib_score_rouge.append(min(scores))
                                distrib_rmsd_rouge.append(rmsd[(nom1, nom2)])
                            
                            if (liste_representant[i], liste_representant[j]) in liste_clique or (liste_representant[j], liste_representant[i]) in liste_clique :
                                distrib_score_clique.append(min(scores))
                                distrib_rmsd_clique.append(rmsd[(nom1, nom2)])
                            
                            if (liste_representant[i], liste_representant[j]) in liste_aretes or (liste_representant[j], liste_representant[i]) in liste_aretes :
                                distrib_score_en_plus.append(min(scores))
                                distrib_rmsd_en_plus.append(rmsd[(nom1, nom2)])

                               
                            if min(scores) <= 60 and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)] <= 3 :
                                csvwriter3.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)]])
                            ## pour clusters cd-hit                      
    #                         for cluster in clusters :
    #                             if liste_representant[i] in cluster and liste_representant[j] in cluster :
    #                                 distrib_min_cluster.append(min(scores))
    #                                 distrib_min_cluster_rmsd.append(rmsd[(nom1, nom2)])
                            
                        else :    
#                             if (liste_representant[i], liste_representant[j]) in dico_min_pos_scores.keys() or (liste_representant[j], liste_representant[i]) in dico_min_pos_scores.keys() : 
#                                 dico_min_pos_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom2, nom1)]})
    
                            distrib_min_rmsd.append(rmsd[(nom2, nom1)])
                            dico_min_tot_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom2, nom1)]})
                            
                            if (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                distrib_score_verte.append(min(scores))
                                distrib_rmsd_verte.append(rmsd[(nom2, nom1)])
                            if (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge :
                                distrib_score_rouge.append(min(scores))
                                distrib_rmsd_rouge.append(rmsd[(nom2, nom1)])
                                
                            if (liste_representant[i], liste_representant[j]) in liste_clique or (liste_representant[j], liste_representant[i]) in liste_clique :
                                distrib_score_clique.append(min(scores))
                                distrib_rmsd_clique.append(rmsd[(nom2, nom1)])
                            
                            if (liste_representant[i], liste_representant[j]) in liste_aretes or (liste_representant[j], liste_representant[i]) in liste_aretes :
                                distrib_score_en_plus.append(min(scores))
                                distrib_rmsd_en_plus.append(rmsd[(nom2, nom1)])
                            
                            if min(scores) <= 60 and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)] <= 3 :
                                csvwriter3.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom2, nom1)]])
                            ## pour clusters cd-hit  
    #                         for cluster in clusters :
    #                             if liste_representant[i] in cluster and liste_representant[j] in cluster :
    #                                 distrib_min_cluster.append(min(scores))
    #                                 distrib_min_cluster_rmsd.append(rmsd[(nom2, nom1)])
    #                         print(liste_representant[i], liste_representant[j])
                        
                        ## pour rechercher les paires en positions similaires sur la sequence 
                        with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0], liste_representant[i][1]), "rb") as fichier_1 :
                            mon_depickler_1 = pickle.Unpickler(fichier_1)
                            graphe1 = mon_depickler_1.load()
                              
                        with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0], liste_representant[j][1]), "rb") as fichier_2 :
                            mon_depickler_2 = pickle.Unpickler(fichier_2)
                            graphe2 = mon_depickler_2.load()
                             
                            print(compter)
                            compter += 1
                             
                            pos_similaire = pos_similaire_fr3d(graphe1, graphe2, liste_representant[i], liste_representant[j], seuil_hom)
                            if liste_representant[i] in [('4y4o', 47), ('6qul',16)] and liste_representant[j] in [('4y4o', 47), ('6qul',16)] :
                                print(pos_similaire)
      
                            #pos_similaire = positions_similaires(graphe1, graphe2, seuil_hom)
                            if pos_similaire[0] :
                                if min(scores) == scores[0] :
                                    dico_type_seq_pos[0] += 1
                                if min(scores) == scores[1] :
                                    dico_type_seq_pos[1] += 1
                                if min(scores) == scores[2] :
                                    dico_type_seq_pos[2] += 1
                                distrib_diff_pos.append(abs(abs(graphe1.nodes[1]["position"][0]-graphe1.nodes[2]["position"][0])-abs(graphe2.nodes[1]["position"][0]-graphe2.nodes[2]["position"][0])))
                                if liste_representant[i][0] != liste_representant[j][0] :
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
                                        #if dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                        if min_ou_max == "min" :
                                            distrib_min_pos_sim.append(min(scores))
                                            dico_min_pos_scores.update({(liste_representant[i], liste_representant[j]) : (min(scores))})
                                            dico_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (min(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                        elif min_ou_max == "max" :
                                            distrib_min_pos_sim.append(max(scores))
                                            dico_min_pos_scores.update({(liste_representant[i], liste_representant[j]) : (max(scores))})
                                            dico_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (max(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                        else :
                                            print("min ou max pas bon")
                                            exit()
                                    else :
                                        #if dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        if min_ou_max == "min" :    
                                            distrib_min_pos_sim.append(min(scores))
                                            dico_min_pos_scores.update({(liste_representant[i], liste_representant[j]) : (min(scores))})
                                            dico_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (min(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                        elif min_ou_max == "max" :
                                            distrib_min_pos_sim.append(max(scores))
                                            dico_min_pos_scores.update({(liste_representant[i], liste_representant[j]) : (max(scores))})
                                            dico_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (max(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                        else :
                                            print("min ou max pas bon")
                                            exit()
                                    if (nom1, nom2) in rmsd.keys() :
                                         
                                        #csvwriter2.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)], pos_similaire])
                                        entree = False
                                        ajoute = False
                                        not_ajoute = False
                                        liste_a_enlever = []
                                        for cle in dico_min_pos_rmsd.keys() :
                                            if not (liste_representant[i] == cle[0] and liste_representant[j] == cle[1]) and not (liste_representant[i] == cle[1] and  liste_representant[j] == cle[0]) and ((liste_representant[i] in cle and (liste_representant[j][0] == cle[0][0] or liste_representant[j][0] == cle[1][0])) or (liste_representant[j] in cle and (liste_representant[i][0] == cle[0][0] or liste_representant[i][0] == cle[1][0])) ):
                                                entree = True
                                                if dico_min_pos_sim[cle][1] > dico_min_pos_sim[(liste_representant[i], liste_representant[j])][1] and dico_min_pos_sim[cle][2] > dico_min_pos_sim[(liste_representant[i], liste_representant[j])][2] and dico_min_pos_sim[cle][3] > dico_min_pos_sim[(liste_representant[i], liste_representant[j])][3] :
                                                    ajoute = True
                                                    liste_a_enlever.append(cle)
                                                    print("roupoulou")
                                                elif dico_min_pos_sim[cle][1] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][1] and dico_min_pos_sim[cle][2] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][2] and dico_min_pos_sim[cle][3] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][3] :
                                                    not_ajoute = True
                                                    del(dico_min_pos_scores[(liste_representant[i], liste_representant[j])])
                                                    del(dico_min_pos_sim[(liste_representant[i], liste_representant[j])])
                                                    print("ripili")
                                                    break
                                                else :
                                                    if dico_min_pos_sim[cle][1] + dico_min_pos_sim[cle][2] + dico_min_pos_sim[cle][3] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][1] + dico_min_pos_sim[(liste_representant[i], liste_representant[j])][2] + dico_min_pos_sim[(liste_representant[i], liste_representant[j])][3] :
                                                        del(dico_min_pos_scores[(liste_representant[i], liste_representant[j])])
                                                        del(dico_min_pos_sim[(liste_representant[i], liste_representant[j])])
                                                        break
                                                    else :
                                                        ajoute = True
                                                        liste_a_enlever.append(cle)
                                                    print("different")
                                        if ajoute :
                                            if not not_ajoute :
                                                distrib_min_pos_rmsd.append(rmsd[(nom1, nom2)])
                                                dico_min_pos_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom1, nom2)]})
                                              
                                            for cle in liste_a_enlever :
                                                del(dico_min_pos_sim[cle])
                                                del(dico_min_pos_scores[cle])
                                                del(dico_min_pos_rmsd[cle])
                                        if not entree :
                                            distrib_min_pos_rmsd.append(rmsd[(nom1, nom2)])
                                            dico_min_pos_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom1, nom2)]})
                                          
                                    else :
                                        #csvwriter2.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)], pos_similaire])
                                        entree = False
                                        ajoute = False
                                        not_ajoute = False
                                        liste_a_enlever = [] 
                                        for cle in dico_min_pos_rmsd.keys() :
                                            if not (liste_representant[i] == cle[0] and liste_representant[j] == cle[1]) and not (liste_representant[i] == cle[1] and  liste_representant[j] == cle[0]) and ((liste_representant[i], liste_representant[j]) != cle and (liste_representant[j], liste_representant[i]) != cle and (liste_representant[i] in cle and (liste_representant[j][0] == cle[0][0] or liste_representant[j][0] == cle[1][0])) or (liste_representant[j] in cle and (liste_representant[i][0] == cle[0][0] or liste_representant[i][0] == cle[1][0]))) :
                                                entree = True
                                                if dico_min_pos_sim[cle][1] > dico_min_pos_sim[(liste_representant[i], liste_representant[j])][1] and dico_min_pos_sim[cle][2] > dico_min_pos_sim[(liste_representant[i], liste_representant[j])][2] and dico_min_pos_sim[cle][3] > dico_min_pos_sim[(liste_representant[i], liste_representant[j])][3] :
                                                    ajoute = True
                                                    liste_a_enlever.append(cle)
                                                    print("roupoulou")
                                                elif dico_min_pos_sim[cle][1] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][1] and dico_min_pos_sim[cle][2] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][2] and dico_min_pos_sim[cle][3] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][3] :
                                                    not_ajoute = True
                                                    del(dico_min_pos_scores[(liste_representant[i], liste_representant[j])])
                                                    del(dico_min_pos_sim[(liste_representant[i], liste_representant[j])])
                                                    print("ripili")
                                                else :
                                                    if dico_min_pos_sim[cle][1] + dico_min_pos_sim[cle][2] + dico_min_pos_sim[cle][3] < dico_min_pos_sim[(liste_representant[i], liste_representant[j])][1] + dico_min_pos_sim[(liste_representant[i], liste_representant[j])][2] + dico_min_pos_sim[(liste_representant[i], liste_representant[j])][3] :
                                                        del(dico_min_pos_scores[(liste_representant[i], liste_representant[j])])
                                                        del(dico_min_pos_sim[(liste_representant[i], liste_representant[j])])
                                                    else :
                                                        ajoute = True
                                                        liste_a_enlever.append(cle)
                                                    print("different")
                                          
                                        if ajoute :
                                            if not not_ajoute :
                                                distrib_min_pos_rmsd.append(rmsd[(nom2, nom1)])
                                                dico_min_pos_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom2, nom1)]})
                                              
                                            for cle in liste_a_enlever :
                                                del(dico_min_pos_sim[cle])
                                                del(dico_min_pos_scores[cle])
                                                del(dico_min_pos_rmsd[cle])
                                                     
                                        if not entree :        
                                                distrib_min_pos_rmsd.append(rmsd[(nom2, nom1)])
                                                dico_min_pos_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom2, nom1)]})
                                      
                            else :
                                  
                                if min(scores) <= 60 and min(scores) >= 20 and (((nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)]  <= 2.5) or ((nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)]  <= 2.5)) :
                                    if (nom1, nom2) in rmsd.keys() :
                                        csvwriter3.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)], pos_similaire])
                                    else :
                                        csvwriter3.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom2, nom1)], pos_similaire])
                                if min_ou_max == "min" :
                                    dico_not_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (min(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                elif min_ou_max == "max" :   
                                    dico_not_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (max(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                else :
                                    print("min ou max pas bon")
                                    exit()
                                                              
                        ## fin de la recherche de paires en positions similaires 
                                
                        if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
                            if dico_sim[(liste_representant[i], liste_representant[j])]["sim"] == 1.0 :
                                if min_ou_max == "min" :
                                    distrib_min.append(min(scores))
                                    dico_min.update({(liste_representant[i], liste_representant[j]) : min(scores) })
                                elif min_ou_max == "max" :
                                    distrib_min.append(max(scores))
                                    dico_min.update({(liste_representant[i], liste_representant[j]) : max(scores) })
                                else :
                                    print("min ou max pas bon")
                                    exit()
                        else :
                            if dico_sim[(liste_representant[j], liste_representant[i])]["sim"] == 1.0 :
                                if min_ou_max == "min" :
                                    distrib_min.append(min(scores))
                                    dico_min.update({(liste_representant[i], liste_representant[j]) : min(scores) })
                                elif min_ou_max == "max" :
                                    distrib_min.append(max(scores))
                                    dico_min.update({(liste_representant[i], liste_representant[j]) : max(scores) })
                                else :
                                    print("min ou max pas bon")
                                    exit()
            
            
            ## distribution des valeurs de score ##
            if len(distrib_min_pos_sim) > 0 :
                min_pos_sim = min(distrib_min_pos_sim)
            else :
                min_pos_sim = 0
            min_tot = min(distrib_min_tot)
            distrib_min_tot.append(-20.0)
            distrib_min_pos_sim.append(-20.0)
            distrib_min_cluster.append(-20.0)
            #print(distrib_min_pos_sim)
            #print(distrib_min_tot)
            #print(liste_seq_min)
            #LogMin, LogMax = np.log10(min(distrib_min_tot)),np.log10(max(distrib_min_tot))
            #newBins = np.logspace(LogMin, LogMax,8)
            #axs = sns.distplot(distrib_min_tot,bins=newBins,kde=False)        
               
            plt.figure(figsize=(9,6))
            axs = sns.distplot(distrib_min_tot, kde=False, bins = len(liste_representant)+1, hist_kws={'log':True})
            sns.distplot(distrib_min_pos_sim, kde=False, bins = len(liste_representant)+1, hist_kws={'log':True})
            sns.distplot(distrib_min_cluster, kde=False, bins = len(liste_representant)+1, hist_kws={'log':True})
            plt.xlim(min_tot-10, 351)
    
            axs.set_xlabel("Score %simum d'alignement global"%min_ou_max)
            axs.set_ylabel("Nombre de paires")
            axs.set_title("Distribution des valeurs de scores d'alignement global \n (minimum des 3 valeurs) pour les %s \n (avec seuil positions similaires : %d min scores : %d)"%(num_ARN, seuil_hom, min_pos_sim))
            plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_%s_seuil_pos_sim_%d.png"%(num_ARN, min_ou_max, seuil_hom))
            plt.show()
            plt.close()
            
            with open("dico_min_tot_score_%s.pickle"%num_ARN, "wb") as fichier_score :
                mon_pickler = pickle.Pickler(fichier_score)
                mon_pickler.dump(dico_min_tot)
                
            with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "wb") as fichier_rmsd :
                mon_pickler = pickle.Pickler(fichier_rmsd)
                mon_pickler.dump(dico_min_tot_rmsd)
            
            with open("dico_min_pos_score_%s_seuil_%s.pickle"%(num_ARN, seuil_hom), "wb") as fichier_score_pos :
                mon_pickler = pickle.Pickler(fichier_score_pos)
                mon_pickler.dump(dico_min_pos_scores) 
                
            #calcul_kmeans_sur_distrib(dico_min_tot, dico_min_tot_rmsd)

#             x = np.arange(55.5,310,5)
#             y = [0.034*xi+0.65 for xi in x]
            
            ## distribution des valeurs de rmsd en fonction des valeurs de score ##
            
            distrib_min_tot.remove(-20.0)
            distrib_min_pos_sim.remove(-20.0)
            distrib_min_cluster.remove(-20.0)
            print(len(dico_min_pos_scores))
            print(dico_min_pos_scores.values())
            print(len(dico_min_pos_rmsd))
            print(dico_min_pos_rmsd.values())
            print(len(distrib_rmsd_en_plus))
            print(len(dico_min_tot.values()))
            plt.figure(figsize=(12,6))
            axs = plt.gca()
            
            plt.scatter(dico_min_tot.values(), dico_min_tot_rmsd.values())
            #plt.scatter(distrib_score_en_plus, distrib_rmsd_en_plus, c='yellow')
            plt.scatter(dico_min_pos_scores.values(), dico_min_pos_rmsd.values())
            plt.scatter(distrib_score_clique, distrib_rmsd_clique, c='purple')
            plt.scatter(distrib_score_rouge, distrib_rmsd_rouge, c='red')
            
            
            plt.scatter(distrib_score_verte, distrib_rmsd_verte, c='green')
            
            plt.scatter(distrib_min_cluster, distrib_min_cluster_rmsd)
            #plt.plot(x,y, c='purple')
            axs.set_xlabel("Score %simum d'alignement global"%min_ou_max)
            axs.set_ylabel("RMSD (en Angstrom)")
            axs.set_xticks(np.arange(0, max(310, max(dico_min_tot.values())+20), 10))
            axs.set_yticks(range(0,16))
            axs.set_title("Distribution des valeurs de RMSD en fonction des valeurs de score d'alignement global \n (seuil positions similaires : %d)"%seuil_hom)
            plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_%s_seuil_pos_rmsd_%d.png"%(num_ARN, min_ou_max, seuil_hom))
            #plt.grid()
            plt.show()
            plt.close()
            
            
            return liste_en_dessous_de_100, liste_en_dessous_de_100_sim, liste_en_dessous_de_100_pos, liste_en_dessous_de_100_rmsd


''' 05/02/20 
-creation de fichiers gephi contenant les sous-graphes induits par les paires de sommets au-dessus/en-dessous des seuils definis en parametre
(pour definition homologues) 
aretes ponderees par score, rmsd, et sim
+ peut contenir des aretes supplementaires de poids plus faible pour obtenir des cliques (ponderees de la meme façon par score, rmsd et sim)
- creation de fichiers gephi contenant les sous-graphes induits par les paires de sommets avec un seuil de sim au-dessus de 0.6 en partant de la liste de sommets definie avec la rmsd et le score du dessus
- renvoie les groupes d'homologues definis ainsi à partir des composantes connexes du graphe de score et de rmsd
et les indique les groupes pour lesquels le clustering perez sur la sim ne met pas tous les elements du groupe dans le même cluster '''
def recup_groupes_homologues(liste_num_ARN, num_ARN, seuil_score, seuil_rmsd, seuil_sim, seuil_eps):
 
    ## creation d'un graphe ponderes par les valeurs de score et de RMSD
    with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
        mon_depickler = pickle.Unpickler(fichier_score)
        dico_min_tot = mon_depickler.load()
                
    with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "rb") as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        dico_min_tot_rmsd =  mon_depickler.load()

        graph = nx.Graph()
        for cle in dico_min_tot.keys() :
            if cle[0] not in liste_pbs and cle[1] not in liste_pbs and cle[0] != '6hrm' and cle[1] != '6hrm' :
                if cle in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[cle] != None :
                    if cle[0] not in graph.nodes() :
                        graph.add_node(cle[0])
                    if cle[1] not in graph.nodes() :
                        graph.add_node(cle[1])
                    graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[cle])
                elif (cle[1], cle[0]) in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[(cle[1], cle[0])] != None :
                    if cle[0] not in graph.nodes() :
                        graph.add_node(cle[0])
                    if cle[1] not in graph.nodes() :
                        graph.add_node(cle[1])
                    graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[(cle[1], cle[0])])
    
        print(graph.number_of_nodes())   

        ## on cree un graphe copy dans lequel on ne garde que les aretes correspondant a notre seuil
        
        a_enlever_aretes = []
        
        for u,v,data in graph.edges(data=True) :
            if data["aln"] < seuil_score or data["rmsd"] > seuil_rmsd :
                a_enlever_aretes.append((u,v))
        
        graph_copy = copy.deepcopy(graph)
        for u,v in a_enlever_aretes :
            graph_copy.remove_edge(u, v)
            
        
        print(graph.number_of_edges())
        print(graph_copy.number_of_edges())
        composantes = recherche_composante_connexe(graph_copy)
        
        ## on supprime les composantes de taille 1
        comp_taille_1 = []
        for composante in composantes :
            if len(composante) == 1 :
                comp_taille_1.append(composante)
                
        for elt in comp_taille_1 :
            composantes.remove(elt)
        
        print(len(composantes))
        
        ## on supprime les noeuds isoles dans le graphe
        a_enlever_noeuds = []
        for noeud in graph.nodes() :
            if len(graph_copy[noeud]) == 0 :
                a_enlever_noeuds.append(noeud)
        
        print(graph.number_of_nodes())
        for noeud in a_enlever_noeuds :
            graph.remove_node(noeud)
            graph_copy.remove_node(noeud)
        print(graph.number_of_nodes())
        print(graph_copy.number_of_nodes())
        print(graph_copy.number_of_edges())
        print(graph.number_of_edges())
        
        ## creation du graphe pondere sur la sim et d'un graphe ne possedant que les aretes au-dessus du seuil de sim
        graphe_complet_sim = creation_graphe_complet(liste_num_ARN)
        graphe_seuil = nx.Graph()
        
        for u,v,data in graphe_complet_sim.edges(data=True) :
            if u in graph.nodes and v in graph.nodes() and u[0] != "6hrm" and v[0] != "6hrm" :
                if data["sim"] >= seuil_sim:
                    graphe_seuil.add_edge(u,v, sim=data["sim"])
                if (u,v) in graph_copy.edges() :
                    graph_copy.edges[u,v]["sim"] = data["sim"]


            
        
        compte_aretes = 0
        for u,v in graph_copy.edges() :
            if (u,v) not in graphe_seuil.edges() :
                compte_aretes += 1
                print(u,v)
        
        print(len(a_enlever_aretes))
        print(graph.number_of_edges())
        print(graph_copy.number_of_edges())
        
        ## clustering dbscan du graphe de rmsd/score
#         matrice = [[2 for _ in range(graph_copy.number_of_nodes())] for _ in range(graph_copy.number_of_nodes())]
#         for u,v, data in graph_copy.edges(data=True) :
#             compte_u = list(graph_copy.nodes()).index(u)
#             compte_v = list(graph_copy.nodes()).index(v)
#             matrice[compte_u][compte_v] = 1-data["sim"]
#             print(compte_u)
#             print(compte_v)
#         #matrice = np.array(matrice)
#         print(matrice)
#         print(graph_copy.nodes.data())
#         print(graph_copy.edges.data())
#         clustering = DBSCAN(min_samples=1, eps=seuil_eps, metric='precomputed').fit(matrice)
#         #clusters, dico_relevance = algo_principal(graph)
#         
#         clusters = [[]]
#         print(len(clusters))
#         print(len(clustering.labels_))
#         print(clustering.labels_)
#         compteur = 0
#         for elt in clustering.labels_ :
#             if elt != -1 :
#                 for _ in range(len(clusters), elt + 1) :
#                         clusters.append([])
#                 clusters[elt].append(list(graph_copy.nodes())[compteur])
#             compteur += 1
#         
#         
#         print(len(clusters))
        
        ## clustering dbscan sur graphe de sim
        matrice_sim = [[2 for _ in range(graphe_seuil.number_of_nodes())] for _ in range(graphe_seuil.number_of_nodes())]
        for u,v, data in graphe_seuil.edges(data=True) :
            compte_u = list(graphe_seuil.nodes()).index(u)
            compte_v = list(graphe_seuil.nodes()).index(v)
            matrice_sim[compte_u][compte_v] = 1-data["sim"] ## 0.0 pour ne pas prendre en compte les valeurs de sim, 1 - data["sim"] si on veut les prendre en compte
        print(graphe_seuil.nodes.data())
        print(graphe_seuil.edges.data())
        clustering_sim = DBSCAN(min_samples=1, eps=seuil_eps, metric='precomputed').fit(matrice_sim)
        #clusters_sim, dico_relevance = algo_principal(graphe_seuil)
        
        clusters_sim = [[]]
        print(len(clusters_sim))
        print(len(clustering_sim.labels_))
        print(clustering_sim.labels_)
        compteur = 0
        for elt in clustering_sim.labels_ :
            if elt != -1 :
                for _ in range(len(clusters_sim), elt + 1) :
                        clusters_sim.append([])
                clusters_sim[elt].append(list(graphe_seuil.nodes())[compteur])
            compteur += 1
        print(graphe_seuil.number_of_nodes())
        print(graph.number_of_nodes())
        
        
        ## calcul Jaccard entre composantes connexes du graphe score/RMSD et clustering dbscan du graphe de sim 
        res = calcul_sim_jaccard(composantes, clusters_sim)
        print(res)
        #print(calcul_sim_jaccard_pour_clustering_recouvrant(clusters, clusters_sim, list(graph.nodes())))
        #return res
        #print(len(clusters))
        
        
        #exit()
        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_sim_%s_seuil_%s.csv"%(num_ARN,num_ARN, str(seuil_score)+"_"+str(seuil_rmsd)+"_"+str(seuil_eps)+"_"+str(seuil_sim)),'w') as fichier_csv:
                    csvwriter = csv.writer(fichier_csv)
                    csvwriter.writerow(["source", "target", "weight", "sim", "dans_rmsd", "aln", "rmsd"])
                    
                    
                    
                    nx.set_edge_attributes(graph, -1, "sim")
                    compteur = 0
                    for u,v,data in graphe_seuil.edges(data=True) :
                        if u in graph.nodes and v in graph.nodes() :
                            #if data["sim"] >= 0.75 :
                            if (u,v) in res[3] or (v,u) in res[3] :
                                compteur += 1
                                csvwriter.writerow([u,v, 2, round(data["sim"], 2), 0, graph.edges[u,v]["aln"], round(graph.edges[u,v]["rmsd"], 2)])
                            else :
                                csvwriter.writerow([u,v, 2, round(data["sim"], 2), 1, graph.edges[u,v]["aln"], round(graph.edges[u,v]["rmsd"], 2)])
                    
                    for u,v in res[3] :
                        if (u,v) not in graphe_seuil.edges() and (v,u) not in graphe_seuil.edges() :
                            csvwriter.writerow([u,v, 1, round(graphe_complet_sim.edges[u,v]["sim"], 2), 0, graph.edges[u,v]["aln"], round(graph.edges[u,v]["rmsd"], 2)])
 
                    
                    print(graph.number_of_nodes())
                    print(graphe_seuil.number_of_nodes())     
                    print(compteur)
                    #exit()
                   
            
            
            
        liste_pour_clique = []          
            #with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_seuil_%s.csv"%(num_ARN,num_ARN, str(seuil3_score)+"_"+str(a_droite)+"_"+str(b_droite)),'w') as fichier_csv:
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_seuil_%s.csv"%(num_ARN,num_ARN, str(seuil_score)+"_"+str(seuil_rmsd)+"_"+str(seuil_eps)+"_"+str(seuil_sim)),'w') as fichier_csv:
                    csvwriter = csv.writer(fichier_csv)
                    csvwriter.writerow(["source", "target", "weight", "aln", "rmsd", "zone_grise", "dans_sim", "val_sim"])
                    compteur = 0
                    for u,v,data in graph.edges(data=True) :
                        
                        
                        if (u,v) in a_enlever_aretes :
                            ok = False
                            for composante in composantes :
                                if u in composante and v in composante :
    #                                 if (u,v) in a_enlever_aretes_perdus :
    #                                     csvwriter.writerow([u,v, 1, data["aln"], round(data["rmsd"],2), 1])
    #                                       
    #                                 else :
                                        #if (u,v) in graphe_seuil.edges() :
                                        if (u,v) in res[2] or (v,u) in res[2] :
                                            csvwriter.writerow([u,v, 1, data["aln"], round(data["rmsd"],2), 0, 0, round(graphe_complet_sim.edges[u,v]["sim"],2)])
                                        else :
                                            csvwriter.writerow([u,v, 1, data["aln"], round(data["rmsd"],2), 0, 1, round(graphe_complet_sim.edges[u,v]["sim"],2)])
                                        liste_pour_clique.append((u,v))
                                        ok = True
                            if ok :
                                compteur += 1
                        else :

                                if (u,v) in res[2] or (v,u) in res[2] :
                                    csvwriter.writerow([u,v, 2, data["aln"], round(data["rmsd"],2), 0, 0, round(graphe_complet_sim.edges[u,v]["sim"],2)])
                                else :
                                    csvwriter.writerow([u,v, 2, data["aln"], round(data["rmsd"],2), 0, 1, round(graphe_complet_sim.edges[u,v]["sim"],2)])
                    
                    for u,v in res[2] :
                        if (u,v) not in graph.edges() and (v,u) not in graph.edges() :
                            csvwriter.writerow([u,v, 0.5, data["aln"], round(data["rmsd"],2), 0, 1, round(graphe_complet_sim.edges[u,v]["sim"],2)])

                    print(compteur)
        
        
        #with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_seuil_%s_noeuds.csv"%(num_ARN, num_ARN, str(seuil3_score)+"_"+str(a_droite)+"_"+str(b_droite)),'w') as fichier_csv:
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_seuil_%s_noeuds.csv"%(num_ARN,num_ARN, str(seuil_score)+"_"+str(seuil_rmsd)+"_"+str(seuil_eps)+"_"+str(seuil_sim)),'w') as fichier_csv:
                    csvwriter = csv.writer(fichier_csv)
                    csvwriter.writerow(["id", "libelle", "num_clusters_dbscan"])
                    
                    compteur = 0
                    for noeud,data in graph.nodes(data=True) :
                        nom_struct = " "
                        #print(noeud)
                        with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                            mon_depickler = pickle.Unpickler(fichier_graphe)
                            graphe = mon_depickler.load()
                            with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                                csvreader = csv.reader(fichier_csv_mol)
                                
                                for row in csvreader :
                                    #print(row)
                                    if len(row) == 7 :
                                        if row[1].lower() == noeud[0] :
                                            if graphe.nodes[1]["num_ch"] == row[2] : 
                                                nom_struct += row[3]
                                            elif graphe.nodes[2]["num_ch"] == row[2] : 
                                                nom_struct += row[3]
                                            elif graphe.nodes[5]["num_ch"] == row[2] :
                                                nom_struct += row[3]
                                                
                        num_clusters = []
                        compteur_clusters = 0
                        for cluster in composantes :
                            if noeud in cluster :
                                num_clusters.append(compteur_clusters)
                            compteur_clusters += 1
                                #print(nom_struct)
                        compteur += 1
                        
                                                
                        
                        csvwriter.writerow([noeud, nom_struct, num_clusters])
                        
                    print(compteur)
        
    
        #exit()
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_sim_%s_seuil_%s_noeuds.csv"%(num_ARN,num_ARN, str(seuil_score)+"_"+str(seuil_rmsd)+"_"+str(seuil_eps)+"_"+str(seuil_sim)),'w') as fichier_csv:
                    csvwriter = csv.writer(fichier_csv)
                    csvwriter.writerow(["id", "libelle", "num_clusters_dbscan", "homologues"])
                    
                    compteur = 0
                    for noeud,data in graph.nodes(data=True) :
                        nom_struct = " "
                        #print(noeud)
                        with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                            mon_depickler = pickle.Unpickler(fichier_graphe)
                            graphe = mon_depickler.load()
                            with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                                csvreader = csv.reader(fichier_csv_mol)
                                
                                for row in csvreader :
                                    #print(row)
                                    if len(row) == 7 :
                                        if row[1].lower() == noeud[0] :
                                            if graphe.nodes[1]["num_ch"] == row[2] : 
                                                nom_struct += row[3]
                                            elif graphe.nodes[2]["num_ch"] == row[2] : 
                                                nom_struct += row[3]
                                            elif graphe.nodes[5]["num_ch"] == row[2] :
                                                nom_struct += row[3]
                                
                                
                        num_clusters = []
                        compteur_clusters = 0
                        for cluster in clusters_sim :
                            if noeud in cluster :
                                num_clusters.append(compteur_clusters)
                            compteur_clusters += 1   
                                #print(nom_struct)
                        compteur += 1
                        
                        num_homologues = -1
                        compteur_homol = 0
                        for cluster in composantes : 
                            if noeud in cluster :
                                num_homologues = compteur_homol
                            compteur_homol += 1   
                                #print(nom_struct)
                        
                                                
                        
                        csvwriter.writerow([noeud, nom_struct, num_clusters, num_homologues])
                        
                    print(compteur)
        return res
    
    
def recup_groupes_non_homologues(liste_num_ARN, num_ARN, seuil_score, seuil_rmsd):  

    with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
        mon_depickler = pickle.Unpickler(fichier_score)
        dico_min_tot = mon_depickler.load()
                
    with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "rb") as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        dico_min_tot_rmsd =  mon_depickler.load()
    #     print(len(dico_min_tot))
    #     print(len(dico_min_tot_rmsd))
    #     print(list(dico_min_tot.keys())[0])
        graph = nx.Graph()
        for cle in dico_min_tot.keys() :
            if cle[0] not in liste_pbs and cle[1] not in liste_pbs :
                if cle in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[cle] != None :
                    if cle[0] not in graph.nodes() :
                        graph.add_node(cle[0])
                    if cle[1] not in graph.nodes() :
                        graph.add_node(cle[1])
                    graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[cle])
                elif (cle[1], cle[0]) in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[(cle[1], cle[0])] != None :
                    if cle[0] not in graph.nodes() :
                        graph.add_node(cle[0])
                    if cle[1] not in graph.nodes() :
                        graph.add_node(cle[1])
                    graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[(cle[1], cle[0])])
    
        print(graph.number_of_nodes())   
            
                
                #if dico_min_tot[cle] >= 71 :
                #nom1 = "fichier_%s_%d_taille_4.pdb"%(cle[0][0], cle[0][1])
                #nom2 = "fichier_%s_%d_taille_4.pdb"%(cle[1][0], cle[1][1])

        a_enlever_aretes = []
        
        for u,v,data in graph.edges(data=True) :           
            if data["aln"] > seuil_score or data["rmsd"] < seuil_rmsd :
                a_enlever_aretes.append((u,v))
        
        
        
        graph_copy = copy.deepcopy(graph)
        for u,v in a_enlever_aretes :
            graph_copy.remove_edge(u, v)
            
        
        print(graph.number_of_edges())
        print(graph_copy.number_of_edges())
        #exit()
           
        #composantes = recherche_composante_connexe(graph_copy)
        #print(len(composantes))
        
        ## noeuds isoles
        a_enlever_noeuds = []
        for noeud in graph.nodes() :
            if len(graph_copy[noeud]) == 0 :
                a_enlever_noeuds.append(noeud)
        print(graph.number_of_nodes())
        for noeud in a_enlever_noeuds :
            graph.remove_node(noeud)
        print(graph.number_of_nodes())
        print(graph_copy.number_of_edges())
        print(graph.number_of_edges())
        #exit()
        
        graphe_complet_sim = creation_graphe_complet(liste_num_ARN)
        graphe_seuil_symetrique = nx.Graph()
        
        for u,v,data in graphe_complet_sim.edges(data=True) :
            if u in graph.nodes() and v in graph.nodes() :
                if data["sim"] < 0.75 :
                    graphe_seuil_symetrique.add_edge(u,v,sim=data["sim"])
        
        compte_aretes = 0    
        liste_aretes_au_dessus = []
        for u,v,data in graph_copy.edges(data=True) :
            if (u,v) not in graphe_seuil_symetrique.edges() :
                compte_aretes += 1
                print(graphe_complet_sim.edges[u,v]["sim"])
                liste_aretes_au_dessus.append((u,v))
                
        print(compte_aretes) 
        print(liste_aretes_au_dessus)
        
        liste_aretes_au_dessus_rmsd = []
        for u,v,data in graphe_seuil_symetrique.edges(data=True) :
            if (u,v) not in graph_copy.edges() :
                compte_aretes += 1
                print(graphe_complet_sim.edges[u,v]["sim"])
                liste_aretes_au_dessus_rmsd.append((u,v))
        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_non_hom_sim_%s_seuil_%s.csv"%(num_ARN,num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["source", "target", "sim", "dans_rmsd", "rmsd", "aln"])
                
                
                compteur = 0
                for u,v,data in graphe_seuil_symetrique.edges(data=True) :
                    if u in graph.nodes() and v in graph.nodes() :
                        if (u,v) in liste_aretes_au_dessus_rmsd or (v,u) in liste_aretes_au_dessus_rmsd :
                            csvwriter.writerow([u,v, round(data["sim"],2), 0, graph.edges[u,v]["rmsd"], graph.edges[u,v]["aln"]])
                        else :
                            csvwriter.writerow([u,v, round(data["sim"],2), 1, graph_copy.edges[u,v]["rmsd"], graph_copy.edges[u,v]["aln"]])
                print(compteur)
                #exit()
               
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_non_hom_sim_%s_seuil_%s_noeuds.csv"%(num_ARN,num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["id", "libelle"])
                
                compteur = 0
                for noeud,data in graphe_seuil_symetrique.nodes(data=True) :
                        nom_struct = -1
                        #print(noeud)
                        with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                            mon_depickler = pickle.Unpickler(fichier_graphe)
                            graphe = mon_depickler.load()
                            with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                                csvreader = csv.reader(fichier_csv_mol)
                                
                                for row in csvreader :
                                    #print(row)
                                    if len(row) == 7 :
                                        if row[1].lower() == noeud[0] :
                                            if graphe.nodes[1]["num_ch"] == row[2] and graphe.nodes[2]["num_ch"] == row[2] and graphe.nodes[5]["num_ch"] == row[2] :
                                                nom_struct = row[3]
                                #print(nom_struct)
                        compteur += 1
                    
                                            
                    
                        csvwriter.writerow([noeud, nom_struct])
                    
                print(compteur)
        
        
        
        
    
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_non_hom_%s_seuil_%s.csv"%(num_ARN,num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["source", "target", "aln", "rmsd", "dans_sim", "sim"])
                
                compteur = 0
                for u,v,data in graph.edges(data=True) :
                    if (u,v) not in a_enlever_aretes :
                        if (u,v) in liste_aretes_au_dessus :
                            csvwriter.writerow([u,v, data["aln"], round(data["rmsd"],2), 0, graphe_complet_sim.edges[u,v]["sim"]])
                        else :
                            csvwriter.writerow([u,v, data["aln"], round(data["rmsd"],2), 1, graphe_seuil_symetrique.edges[u,v]["sim"]])
                print(compteur)
                #exit()
               
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_non_hom_%s_seuil_%s_noeuds.csv"%(num_ARN,num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["id", "libelle"])
                
                compteur = 0
                for noeud,data in graph.nodes(data=True) :
                    nom_struct = -1
                    #print(noeud)
                    with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        graphe = mon_depickler.load()
                        with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                            csvreader = csv.reader(fichier_csv_mol)
                            
                            for row in csvreader :
                                #print(row)
                                if len(row) == 7 :
                                    if row[1].lower() == noeud[0] :
                                        if graphe.nodes[1]["num_ch"] == row[2] and graphe.nodes[2]["num_ch"] == row[2] and graphe.nodes[5]["num_ch"] == row[2] :
                                            nom_struct = row[3]
                            #print(nom_struct)
                    compteur += 1
                    
                                            
                    
                    csvwriter.writerow([noeud, nom_struct])
                    
                print(compteur)
                
        a_enlever_a_enlever_aretes = []
        for (u,v) in a_enlever_aretes :
            if u in graph.nodes() and v in graph.nodes() :
                a_enlever_a_enlever_aretes.append((u,v))
        print(len(a_enlever_a_enlever_aretes))
        print(len(a_enlever_aretes))
        return a_enlever_a_enlever_aretes
                
                


def recup_certaines_aretes_inter_types(liste_num_ARN, num_ARN1, num_ARN2, seuil_score, seuil_rmsd):  

    with open("dico_min_tot_score_%s_%s.pickle"%(num_ARN1, num_ARN2), "rb") as fichier_score :
        mon_depickler = pickle.Unpickler(fichier_score)
        dico_min_tot = mon_depickler.load()
                
    with open("dico_min_tot_rmsd_%s_%s.pickle"%(num_ARN1, num_ARN2), "rb") as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        dico_min_tot_rmsd =  mon_depickler.load()
    #     print(len(dico_min_tot))
    #     print(len(dico_min_tot_rmsd))
    #     print(list(dico_min_tot.keys())[0])
        graph = nx.Graph()
        for cle in dico_min_tot.keys() :
            if cle[0] not in liste_pbs and cle[1] not in liste_pbs :
                if cle in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[cle] != None :
                    if cle[0] not in graph.nodes() :
                        graph.add_node(cle[0])
                    if cle[1] not in graph.nodes() :
                        graph.add_node(cle[1])
                    graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[cle])
                elif (cle[1], cle[0]) in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[(cle[1], cle[0])] != None :
                    if cle[0] not in graph.nodes() :
                        graph.add_node(cle[0])
                    if cle[1] not in graph.nodes() :
                        graph.add_node(cle[1])
                    graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[(cle[1], cle[0])])
        print("gros tas")
        print(graph.number_of_nodes())   
            
                
                #if dico_min_tot[cle] >= 71 :
                #nom1 = "fichier_%s_%d_taille_4.pdb"%(cle[0][0], cle[0][1])
                #nom2 = "fichier_%s_%d_taille_4.pdb"%(cle[1][0], cle[1][1])

        a_enlever_aretes = []
        
        for u,v,data in graph.edges(data=True) :           
            if data["rmsd"] > seuil_rmsd :
                a_enlever_aretes.append((u,v))
                
        graph_copy = copy.deepcopy(graph)
        for u,v in a_enlever_aretes :
            graph_copy.remove_edge(u, v)
            
        
        print(graph.number_of_edges())
        print(graph_copy.number_of_edges())
        #exit()
           
        #composantes = recherche_composante_connexe(graph_copy)
        #print(len(composantes))
        
        ## noeuds isoles
        a_enlever_noeuds = []
        for noeud in graph.nodes() :
            if len(graph_copy[noeud]) == 0 :
                a_enlever_noeuds.append(noeud)
        
        print(graph.number_of_nodes())
        for noeud in a_enlever_noeuds :
            graph.remove_node(noeud)
        print(graph.number_of_nodes())
        print(graph_copy.number_of_edges())
        print(graph.number_of_edges())
        #exit()
        
        graphe_complet_sim = creation_graphe_complet(liste_num_ARN)
        graphe_seuil = nx.Graph()
        print("tadaaa")
        print(graphe_complet_sim.number_of_nodes())
        #print(graphe_complet_sim.nodes[('4y4o', 18)])
        
        for u,v,data in graphe_complet_sim.edges(data=True) :
            if u in graph.nodes() and v in graph.nodes() :
                if data["sim"] >= 0.6 :
                    graphe_seuil.add_edge(u,v,sim=data["sim"])
        print("grrr")
        print(graphe_seuil.number_of_nodes())
        compte_aretes = 0    
        liste_aretes_au_dessus = []
        for u,v,data in graph_copy.edges(data=True) :
            if (u,v) not in graphe_seuil.edges() :
                compte_aretes += 1
                liste_aretes_au_dessus.append((u,v))
                
        print(compte_aretes) 
        print(liste_aretes_au_dessus)
        
        liste_aretes_au_dessus_rmsd = []
        for u,v,data in graphe_seuil.edges(data=True) :
            if (u,v) not in graph_copy.edges() :
                compte_aretes += 1
                print(graphe_complet_sim.edges[u,v]["sim"])
                liste_aretes_au_dessus_rmsd.append((u,v))
        
        liste_noeuds = set()
        with open("Nouvelles_donnees/alignements_%s_%s/fichier_csv_sim_%s_%s_seuil_%s.csv"%(num_ARN1, num_ARN2, num_ARN1,num_ARN2, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["source", "target", "sim", "rmsd", "dans_rmsd"])
                
                
                compteur = 0
                for u,v,data in graphe_seuil.edges(data=True) :
                    if u in graph.nodes() and v in graph.nodes() and ((u,v) in dico_min_tot_rmsd.keys() or (v,u) in dico_min_tot_rmsd.keys() ) :
                        if (u,v) in liste_aretes_au_dessus_rmsd or (v,u) in liste_aretes_au_dessus_rmsd :
                            print(u,v)
                            #print(dico_min_tot_rmsd[("fichier_%s_%s_taille_4.pdb"%(u[0], u[1]), "fichier_%s_%s_taille_4.pdb"%(v[0], v[1]))] if ("fichier_%s_%s_taille_4.pdb"%(u[0], u[1]), "fichier_%s_%s_taille_4.pdb"%(v[0], v[1])) in dico_min_tot_rmsd.keys() else dico_min_tot_rmsd[("fichier_%s_%s_taille_4.pdb"%(v[0], v[1]), "fichier_%s_%s_taille_4.pdb"%(u[0], u[1]))] )
                            csvwriter.writerow([u,v, round(data["sim"],2), round(graph.edges[u,v]["rmsd"],2), 0])
                        else :
                            csvwriter.writerow([u,v, round(data["sim"],2), round(graph.edges[u,v]["rmsd"],2), 1])
                        liste_noeuds.add(u)
                        liste_noeuds.add(v)
                print(compteur)
                #exit()
               
        with open("Nouvelles_donnees/alignements_%s_%s/fichier_csv_sim_%s_%s_seuil_%s_noeuds.csv"%(num_ARN1, num_ARN2, num_ARN1,num_ARN2, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["id", "libelle"])
                
                compteur = 0
                for noeud,data in graphe_seuil.nodes(data=True) :
                    if noeud in liste_noeuds :
                        nom_struct = -1
                        #print(noeud)
                        with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                            mon_depickler = pickle.Unpickler(fichier_graphe)
                            graphe = mon_depickler.load()
                            with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                                csvreader = csv.reader(fichier_csv_mol)
                                
                                for row in csvreader :
                                    #print(row)
                                    if len(row) == 7 :
                                        if row[1].lower() == noeud[0] :
                                            if graphe.nodes[1]["num_ch"] == row[2] and graphe.nodes[2]["num_ch"] == row[2] and graphe.nodes[5]["num_ch"] == row[2] :
                                                nom_struct = row[3]
                                #print(nom_struct)
                        compteur += 1
                    
                                            
                    
                        csvwriter.writerow([noeud, nom_struct])
                    
                print(compteur)
        

        print(compte_aretes) 
        print(liste_aretes_au_dessus)
        
    
        with open("Nouvelles_donnees/alignements_%s_%s/fichier_csv_%s_%s_seuil_%s.csv"%(num_ARN1, num_ARN2, num_ARN1,num_ARN2, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["source", "target", "aln", "rmsd", "dans_sim", "val_sim"])
                
                compteur = 0
                for u,v,data in graph.edges(data=True) :
                    if (u,v) not in a_enlever_aretes :
                        if (u,v) in liste_aretes_au_dessus :
                            print(u,v)
                            csvwriter.writerow([u,v, data["aln"], round(data["rmsd"],2), 0, round(graphe_complet_sim.edges[u,v]["sim"],2)])
                        else :
                            csvwriter.writerow([u,v, data["aln"], round(data["rmsd"],2), 1, round(graphe_complet_sim.edges[u,v]["sim"],2)])
                print(compteur)
                #exit()
               
        with open("Nouvelles_donnees/alignements_%s_%s/fichier_csv_%s_%s_seuil_%s_noeuds.csv"%(num_ARN1, num_ARN2, num_ARN1,num_ARN2, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["id", "libelle"])
                
                compteur = 0
                for noeud,data in graph.nodes(data=True) :
                    nom_struct = -1
                    #print(noeud)
                    with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        graphe = mon_depickler.load()
                        with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                            csvreader = csv.reader(fichier_csv_mol)
                            
                            for row in csvreader :
                                #print(row)
                                if len(row) == 7 :
                                    if row[1].lower() == noeud[0] :
                                        if graphe.nodes[1]["num_ch"] == row[2] and graphe.nodes[2]["num_ch"] == row[2] and graphe.nodes[5]["num_ch"] == row[2] :
                                            nom_struct = row[3]
                            #print(nom_struct)
                    compteur += 1
                    
                                            
                    
                    csvwriter.writerow([noeud, nom_struct])
                    
                print(compteur)
      
                
def recup_certaines_aretes_plusieurs_types(liste_num_ARN, num_ARN, seuil_score, seuil_rmsd):  
    
    if isinstance(num_ARN, list) :
        print(num_ARN)
        dico_pos = {}
        dico_min_tot_rmsd = {}
        dico_min_tot = {}
        if num_ARN in [['23S', '25S', '28S'], ['16S', '18S']] :
            print("rapala")
            with open("dico_min_pos_score_%s_seuil_30.pickle"%(num_ARN), "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos.update(mon_depickler.load())
                
        with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "rb") as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            dico_min_tot_rmsd.update(mon_depickler.load())   
                
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
            mon_depickler = pickle.Unpickler(fichier_score)
            dico_min_tot.update(mon_depickler.load())
            
        for elt in num_ARN :
            
            #liste_cles.update({elt : []})
            if isinstance(elt, list) :
                print(elt)
                if elt in [['23S', '25S', '28S'], ['16S', '18S']] :
                    print("rapala")
                    with open("dico_min_pos_score_%s_seuil_30.pickle"%(elt), "rb") as fichier_score_pos :
                        mon_depickler = pickle.Unpickler(fichier_score_pos)
                        dico_pos.update(mon_depickler.load())
                        
                with open("dico_min_tot_rmsd_%s.pickle"%elt, "rb") as fichier_rmsd :
                    mon_depickler = pickle.Unpickler(fichier_rmsd)
                    dico_min_tot_rmsd.update(mon_depickler.load())   
                        
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_min_tot.update(mon_depickler.load())

            else :
                print(elt)
                with open("dico_min_tot_rmsd_%s.pickle"%elt, "rb") as fichier_rmsd :
                    mon_depickler = pickle.Unpickler(fichier_rmsd)
                    dico_min_tot_rmsd.update(mon_depickler.load())   
                        
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_min_tot.update(mon_depickler.load())

                #liste_cles[elt].extend(liste_repr)
        
    else :
                
        with open("dico_min_pos_score_%s_seuil_30.pickle"%num_ARN, "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos = mon_depickler.load()
                
        with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "rb") as fichier_rmsd :
                mon_depickler = pickle.Unpickler(fichier_rmsd)
                dico_min_tot_rmsd = mon_depickler.load()  
                
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
                mon_depickler = pickle.Unpickler(fichier_score)
                dico_min_tot = mon_depickler.load()
    
    
    #     print(len(dico_min_tot))
    #     print(len(dico_min_tot_rmsd))
    #     print(list(dico_min_tot.keys())[0])
    graph = nx.Graph()
    for cle in dico_min_tot.keys() :
        if cle in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[cle] != None :
            if cle[0] not in graph.nodes() :
                graph.add_node(cle[0])
            if cle[1] not in graph.nodes() :
                graph.add_node(cle[1])
            graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[cle])
        elif (cle[1], cle[0]) in dico_min_tot_rmsd.keys() and dico_min_tot_rmsd[(cle[1], cle[0])] != None :
            if cle[0] not in graph.nodes() :
                graph.add_node(cle[0])
            if cle[1] not in graph.nodes() :
                graph.add_node(cle[1])
            graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=dico_min_tot_rmsd[(cle[1], cle[0])])

    print(graph.number_of_nodes())   
        
            
            #if dico_min_tot[cle] >= 71 :
            #nom1 = "fichier_%s_%d_taille_4.pdb"%(cle[0][0], cle[0][1])
            #nom2 = "fichier_%s_%d_taille_4.pdb"%(cle[1][0], cle[1][1])

    a_enlever_aretes = []
    
    for u,v,data in graph.edges(data=True) :           
        if data["rmsd"] > seuil_rmsd :
            a_enlever_aretes.append((u,v))
            
    graph_copy = copy.deepcopy(graph)
    for u,v in a_enlever_aretes :
        graph_copy.remove_edge(u, v)
        
    
    print(graph.number_of_edges())
    print(graph_copy.number_of_edges())
    #exit()
       
    #composantes = recherche_composante_connexe(graph_copy)
    #print(len(composantes))
    
    ## noeuds isoles
    a_enlever_noeuds = []
    for noeud in graph.nodes() :
        if len(graph_copy[noeud]) == 0 :
            a_enlever_noeuds.append(noeud)
    
    print(graph.number_of_nodes())
    for noeud in a_enlever_noeuds :
        graph.remove_node(noeud)
    print(graph.number_of_nodes())
    print(graph_copy.number_of_edges())
    print(graph.number_of_edges())
    #exit()
    
    graphe_complet_sim = creation_graphe_complet(liste_num_ARN)
    graphe_seuil = nx.Graph()
    liste_noeuds = set()
    with open("Nouvelles_donnees/alignements_%s/fichier_csv_sim_%s_seuil_%s.csv"%(num_ARN, num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["source", "target", "sim"])
            
            
            compteur = 0
            for u,v,data in graphe_complet_sim.edges(data=True) :
                if u in graph.nodes() and v in graph.nodes() :
                    if data["sim"] >= 0.6 :
                        graphe_seuil.add_edge(u,v,sim=data["sim"])
                        liste_noeuds.add(u)
                        liste_noeuds.add(v)
                        csvwriter.writerow([u,v, round(data["sim"],2)])
            print(compteur)
            #exit()
           
    with open("Nouvelles_donnees/alignements_%s/fichier_csv_sim_%s_seuil_%s_noeuds.csv"%(num_ARN, num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["id", "libelle"])
            
            compteur = 0
            for noeud,data in graphe_complet_sim.nodes(data=True) :
                if noeud in liste_noeuds :
                    nom_struct = -1
                    #print(noeud)
                    with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        graphe = mon_depickler.load()
                        with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                            csvreader = csv.reader(fichier_csv_mol)
                            
                            for row in csvreader :
                                #print(row)
                                if len(row) == 7 :
                                    if row[1].lower() == noeud[0] :
                                        if graphe.nodes[1]["num_ch"] == row[2] and graphe.nodes[2]["num_ch"] == row[2] and graphe.nodes[5]["num_ch"] == row[2] :
                                            nom_struct = row[3]
                            #print(nom_struct)
                    compteur += 1
                
                                        
                
                    csvwriter.writerow([noeud, nom_struct])
                
            print(compteur)
    
    
    compte_aretes = 0    
    liste_aretes_au_dessus = []
    for u,v,data in graph_copy.edges(data=True) :
        if (u,v) not in graphe_seuil.edges() :
            compte_aretes += 1
            print(graphe_complet_sim.edges[u,v]["sim"])
            liste_aretes_au_dessus.append((u,v))
            
    print(compte_aretes) 
    print(liste_aretes_au_dessus)
    

    with open("Nouvelles_donnees/alignements_%s/fichier_csv_%s_seuil_%s.csv"%(num_ARN, num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["source", "target", "aln", "rmsd", "dans_sim", "val_sim"])
            
            compteur = 0
            for u,v,data in graph.edges(data=True) :
                if (u,v) not in a_enlever_aretes :
                    if (u,v) in liste_aretes_au_dessus :
                        csvwriter.writerow([u,v, data["aln"], round(data["rmsd"],2), 1, round(graphe_complet_sim.edges[u,v]["sim"],2)])
                    else :
                        csvwriter.writerow([u,v, data["aln"], round(data["rmsd"],2), 0, round(graphe_complet_sim.edges[u,v]["sim"],2)])
            print(compteur)
            #exit()
           
    with open("Nouvelles_donnees/alignements_%s/fichier_csv_%s_seuil_%s_noeuds.csv"%(num_ARN, num_ARN, str(seuil_score)+'_'+str(seuil_rmsd)),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["id", "libelle"])
            
            compteur = 0
            for noeud,data in graph.nodes(data=True) :
                nom_struct = -1
                #print(noeud)
                with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                    mon_depickler = pickle.Unpickler(fichier_graphe)
                    graphe = mon_depickler.load()
                    with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                        csvreader = csv.reader(fichier_csv_mol)
                        
                        for row in csvreader :
                            #print(row)
                            if len(row) == 7 :
                                if row[1].lower() == noeud[0] :
                                    if graphe.nodes[1]["num_ch"] == row[2] and graphe.nodes[2]["num_ch"] == row[2] and graphe.nodes[5]["num_ch"] == row[2] :
                                        nom_struct = row[3]
                        #print(nom_struct)
                compteur += 1
                
                                        
                
                csvwriter.writerow([noeud, nom_struct])
                
            print(compteur)              


                

    
def traitement_aln_global_inter_types(liste_num_ARN, min_ou_max, seuil_rmsd):
    liste_representant = []
    compteur = 0
    for num_ARN in liste_num_ARN :
        if len(liste_representant) < compteur+1 :
            liste_representant.append([])
        if isinstance(num_ARN, list) :

            for elt in num_ARN :
                with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_repr = mon_depickler_1.load()
                    
                    liste_representant[compteur].extend(liste_repr)
            
        else :
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_representant[compteur] = mon_depickler_1.load()           
        compteur += 1
                
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
        
    distrib_min_tot = []
    dico_min_tot = {}
    dico_type_seq_haut = [0,0,0]
    dico_type_seq = [0,0,0]
    dico_type_seq_bas = [0,0,0]
    distrib_min_rmsd = []
    dico_min_tot_rmsd = {}
    for k in range(len(liste_representant)) :
        for l in range(k+1, len(liste_representant)) :
            for i in range(len(liste_representant[k])) :
                        for j in range(len(liste_representant[l])) :
                            if liste_representant[k][i][0] != '6hrm' and  liste_representant[l][j][0] != '6hrm' and not (liste_representant[k][i] == ('2r8s', 1) and liste_representant[l][j] == ('6d8o', 2)) and not (liste_representant[k][i] == ('1gid', 1) and liste_representant[l][j] == ('6d8o', 2)) and not (liste_representant[k][i] == ('1l8v', 1) and liste_representant[l][j] == ('6d8o', 2)):
                                scores = []
                                for c in range(1,4) :
                                    scores.append(recup_score("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq%d.txt"%(liste_num_ARN,liste_representant[k][i][0],liste_representant[k][i][1], liste_representant[l][j][0],liste_representant[l][j][1], c)))
                                nom1 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[k][i][0], liste_representant[k][i][1])
                                nom2 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[l][j][0], liste_representant[l][j][1])
                                if min_ou_max == "min" :
                                    if (nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)] <= seuil_rmsd :
                                        
                                        distrib_min_tot.append(min(scores))
                                        dico_min_tot.update({(liste_representant[k][i], liste_representant[l][j]) : min(scores)})
                                    elif (nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)] <= seuil_rmsd :
                                        distrib_min_tot.append(min(scores))
                                        dico_min_tot.update({(liste_representant[k][i], liste_representant[l][j]) : min(scores)})
                                    if min(scores) == scores[0] :
                                        if scores[0] >= 100 :
                                            dico_type_seq_haut[0] += 1
                                        else :
                                            dico_type_seq_bas[0] += 1
                                        dico_type_seq[0] += 1
                                    if min(scores) == scores[1] :
                                        if scores[1] >= 100 :
                                            dico_type_seq_haut[1] += 1
                                        else :
                                            dico_type_seq_bas[1] += 1
                                        dico_type_seq[1] += 1
                                    if min(scores) == scores[2] :
                                        if scores[2] >= 100 :
                                            dico_type_seq_haut[2] += 1
                                        else :
                                            dico_type_seq_bas[2] += 1
                                        dico_type_seq[2] += 1
                                        
                                elif min_ou_max == "max" :
                                    distrib_min_tot.append(max(scores))
                                    dico_min_tot.update({(liste_representant[k][i], liste_representant[l][j]) : max(scores)})
                                else :
                                    print("min ou max pas bon")
                                    exit()
                                
                                
                                
                                if (nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)] <= seuil_rmsd :
            
                                    distrib_min_rmsd.append(rmsd[(nom1, nom2)])
                                    dico_min_tot_rmsd.update({(liste_representant[k][i], liste_representant[l][j]) : rmsd[(nom1, nom2)]})
            
                                    
                                elif (nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)] <= seuil_rmsd :
            
                                    distrib_min_rmsd.append(rmsd[(nom2, nom1)])
                                    dico_min_tot_rmsd.update({(liste_representant[k][i], liste_representant[l][j]) : rmsd[(nom2, nom1)]})

                                if min(scores) > 100 : 
                                    print(liste_representant[k][i], liste_representant[l][j])
                                    print(min(scores))
                                    print(rmsd[(nom1, nom2)])
        
    plt.figure(figsize=(9,6))
    axs = sns.distplot(distrib_min_tot, kde=False, bins = len(liste_representant[0][0])+1, hist_kws={'log':True})
    plt.xlim(0, 351)

    axs.set_xlabel("Score %simum d'alignement global"%min_ou_max)
    axs.set_ylabel("Nombre de paires")
    axs.set_title("Distribution des valeurs de scores d'alignement global \n (minimum des 3 valeurs) pour les inter-types \n %s"%(liste_num_ARN))
    #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_%s.png"%(liste_num_ARN, min_ou_max))
    
#     with open("dico_min_tot_score_%s.pickle"%num_ARN, "wb") as fichier_score :
#         mon_pickler = pickle.Pickler(fichier_score)
#         mon_pickler.dump(dico_min_tot)
#         
#     with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "wb") as fichier_rmsd :
#         mon_pickler = pickle.Pickler(fichier_rmsd)
#         mon_pickler.dump(dico_min_tot_rmsd)
        
    #calcul_kmeans_sur_distrib(dico_min_tot, dico_min_tot_rmsd)
    print(dico_type_seq)
    print(dico_type_seq_bas)
    print(dico_type_seq_haut)
    print(len(dico_min_tot))
    
    plt.show()
    plt.close()

    plt.figure(figsize=(12,6))
    axs = plt.gca()
    plt.scatter(dico_min_tot.values(), dico_min_tot_rmsd.values())
    axs.set_xlabel("Score %simum d'alignement global"%min_ou_max)
    axs.set_ylabel("RMSD (en Angstrom)")
    axs.set_xticks(np.arange(0, max(310, max(dico_min_tot.values())+20), 10))
    axs.set_title("Distribution des valeurs de RMSD en fonction des valeurs de score d'alignement global")
    #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_%s_rmsd.png"%(liste_num_ARN, min_ou_max))
    plt.show()
    plt.close()
    
#     with open("dico_min_tot_score_%s.pickle"%(liste_num_ARN), "wb") as fichier_score :
#         mon_pickler = pickle.Pickler(fichier_score)
#         mon_pickler.dump(dico_min_tot)
#                 
#     with open("dico_min_tot_rmsd_%s.pickle"%(liste_num_ARN), "wb") as fichier_rmsd :
#         mon_pickler = pickle.Pickler(fichier_rmsd)
#         mon_pickler.dump(dico_min_tot_rmsd)
    
    
    print(dico_min_tot.values())
    print(max(dico_min_tot.values()))
    
    
    

''' 03/02/20 '''
def distrib_sim_rmsd(num_ARN, liste_homologues, seuil_rmsd, seuil_sim):
    for groupe in liste_homologues :
        if ('2r8s', 1) in groupe :
            groupe.append(('6d8o', 2))
    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%0, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
        
    with open("/media/coline/Maxtor/dico_algo_heuristique_grands_graphes_taille_4.pickle", 'rb') as fichier_sim :
        mon_depickler = pickle.Unpickler(fichier_sim)
        dico_sim = mon_depickler.load()
        
    #liste_cles = {}
    if isinstance(num_ARN, list) :
        
        liste_representant = []
        dico_pos = {}
        dico_score = {}
        
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
            mon_depickler = pickle.Unpickler(fichier_score)
            dico_score.update(mon_depickler.load())
        
        if num_ARN in [['23S', '25S', '28S'], ['16S', '18S']] :
            print("rapala")
            with open("dico_min_pos_score_%s_seuil_30.pickle"%(num_ARN), "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos.update(mon_depickler.load())
                
            
        for elt in num_ARN :
            
            #liste_cles.update({elt : []})
            if isinstance(elt, list) :
                
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_score.update(mon_depickler.load()) 
                
                if elt in [['23S', '25S', '28S'], ['16S', '18S']] :
                    print("rapala")
                    with open("dico_min_pos_score_%s_seuil_30.pickle"%(elt), "rb") as fichier_score_pos :
                        mon_depickler = pickle.Unpickler(fichier_score_pos)
                        dico_pos.update(mon_depickler.load())
                        
                    
                for e in elt :
                    with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%e, 'rb') as fichier_num_arn :
                        mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                        liste_repr = mon_depickler_1.load()
                        
                        liste_representant.extend(liste_repr)
            else :
                with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_representant.extend(mon_depickler_1.load())
                    
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_score.update(mon_depickler.load()) 
                #liste_cles[elt].extend(liste_repr)
        
    else :
        #liste_cles.update({num_ARN : []})
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant = mon_depickler_1.load()
                #liste_cles[num_ARN].extend(liste_representant)
                
        with open("dico_min_pos_score_%s_seuil_30.pickle"%num_ARN, "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos = mon_depickler.load()
                
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_score = mon_depickler.load()

    with open("fichier_forte_sim_faible_rmsd_non_homol_rapoulou.pickle", "w", newline="") as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(["Paire","Type ARN 1", "Type ARN 2", "Sim", "RMSD", "Score d'alignement"])
    
        # 23S
        liste_verte = [(('1vqo', 12), ('6nd6', 7)),(('4y4o', 40), ('4u4r', 8)),(('6eri', 9), ('4u4r', 8)), (('6hma', 15), ('1vq8', 15)), (('5dm6', 2), ('6eri', 1)), (('4u27', 5), ('1vqp', 19)),(('4u27', 5), ('3ccr', 4)),(('4u27', 5), ('3cc7', 4)),(('5dm6', 7), ('1vq8', 11)),(('5dm6', 7), ('3ccr', 12)),(('5dm6', 7), ('3ccl', 11)),(('5dm6', 7), ('3ccs', 12)),(('5dm6', 7), ('3cc7', 12)),(('5dm6', 7), ('3ccq', 12)),(('5dm6', 7), ('3cce', 9)),(('5dm6', 7), ('3ccu', 11)),(('4ybb', 8), ('1vq8', 11)),(('4ybb', 8), ('3ccr', 12)),(('4ybb', 8), ('3ccl', 11)),(('4ybb', 8), ('3ccs', 12)),(('4ybb', 8), ('3cc7', 12)),(('4ybb', 8), ('3ccq', 12)),(('4ybb', 8), ('3cce', 9)),(('4ybb', 8), ('3ccu', 11)),(('3cc2', 2), ('6hma', 1)),(('5dm6', 13), ('1vqp', 19)),(('5dm6', 13), ('3ccr', 4)),(('5dm6', 13), ('3cc7', 4)), (('6hma', 9), ('4u4r', 8)), (('4y4o', 25), ('1vq8', 5)),(('4y4o', 25), ('3cc2', 2)),(('6hma', 15), ('4u4r', 16)), (('1vqp', 19), ('6eri', 9)), (('3ccr', 4), ('6eri', 9)),(('6eri', 9), ('3cc7', 4)),(('5afi', 15), ('4u4r', 28)),(('6hma', 9), ('1vqp', 19)),(('6hma', 9), ('3ccr', 4)),(('6hma', 9), ('3cc7', 4)),(('1vq8', 11), ('6hma', 6)),(('3ccr', 12), ('6hma', 6)),(('6hma', 6), ('3ccl', 11)),(('6hma', 6), ('3ccs', 12)),(('6hma', 6), ('3cc7', 12)),(('6hma', 6), ('3ccq', 12)),(('6hma', 6), ('3cce', 9)),(('6hma', 6), ('3ccu', 11)),(('4y4o', 40), ('1vqp', 19)),(('4y4o', 40), ('3ccr', 4)),(('4y4o', 40), ('3cc7', 4)),(('4y4o', 3), ('3ccl', 11)),(('4y4o', 3), ('3cce', 9)),(('4v51', 22), ('3ccl', 11)),(('4v51', 22), ('3cce', 9)),(('5e81', 12), ('3ccl', 11)),(('5e81', 12), ('3cce', 9)),(('3ccl', 11), ('4v90', 13)),(('4v90', 13), ('3cce', 9)),(('1vqo', 12), ('4u4r', 15)),(('4ybb', 1), ('1vq8', 5)),(('4y4o', 3), ('3ccr', 12)),(('4y4o', 3), ('3ccs', 12)),(('4y4o', 3), ('3cc7', 12)),(('4y4o', 3), ('3ccu', 11)),(('4v51', 22), ('3ccr', 12)),(('4v51', 22), ('3ccs', 12)),(('4v51', 22), ('3cc7', 12)),(('4v51', 22), ('3ccu', 11)),(('5e81', 12), ('3ccr', 12)),(('5e81', 12), ('3ccs', 12)),(('5e81', 12), ('3cc7', 12)),(('5e81', 12), ('3ccu', 11)),(('3ccr', 12), ('4v90', 13)),(('4v90', 13), ('3ccs', 12)),(('4v90', 13), ('3cc7', 12)),(('4v90', 13), ('3ccu', 11)),(('4y4o', 3), ('3ccq', 12)),(('4v51', 22), ('3ccq', 12)),(('5e81', 12), ('3ccq', 12)),(('4v90', 13), ('3ccq', 12)),(('4ybb', 1), ('3cc2', 2)),(('4y4o', 3), ('1vq8', 11)),(('4v51', 22), ('1vq8', 11)),(('5e81', 12), ('1vq8', 11)),(('1vq8', 5), ('6hma', 1)),(('1vq8', 11), ('4v90', 13)), (('5dm6', 9), ('1vq8', 21)), (('6hma', 2), ('4u4r', 19)), (('4ybb', 6), ('1yhq', 24)), (('6hma', 15), ('3cc2', 6))]
        # 16S
        liste_verte.extend([(('5nwy', 17), ('4y4o', 4)), (('6ek0', 2), ('6az1', 1)), (('4u27', 54), ('4y4o', 22)), (('4ybb', 35), ('4u3u', 22))])
        ## homologues cherches a la main avec bonne sim et rmsd
        liste_verte.extend([(('4y4o', 13), ('4u3u', 22)),(('3cc2', 1), ('4u4r', 36)),(('1vqo', 14), ('4u4r', 36)),(('6eri', 3), ('4u3u', 2)),(('4ybb', 6), ('4u3u', 2)),(('4y4o', 9), ('4u4r', 11)),(('4v67', 23), ('4u4r', 11)),(('4v67', 23), ('6ek0', 6)),(('4y4o', 9), ('6ek0', 6)),(('4y4o', 8), ('3cc2', 1)),(('3t1y', 8), ('4u4r', 11)),(('3t1y', 8), ('6ek0', 6)),(('2qex', 19), ('4u3u', 2)),(('1yhq', 24), ('4u3u', 2)),(('4y4o', 13), ('6az1', 6)),(('4y4o', 8), ('1vqo', 14)),(('5ngm', 1), ('4u4r', 11)),(('6eri', 15), ('4u4r', 11)),(('5ngm', 1), ('6ek0', 6)),(('6eri', 15), ('6ek0', 6)),(('4ybb', 9), ('4u4r', 11)),(('5j7l', 10), ('4u4r', 11)),(('4u27', 32), ('4u4r', 11)),(('5j7l', 10), ('6ek0', 6)),(('4ybb', 9), ('6ek0', 6)),(('4u27', 32), ('6ek0', 6)),(('1vq8', 19), ('4u3u', 2)),(('4ybb', 35), ('6az1', 6)),(('2vqe', 5), ('4u4r', 11)),(('2vqe', 5), ('6ek0', 6))])
        
        ## homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim  a 0.75 et pas avec la rmsd a 2
        liste_verte.extend([(('4y4o', 25),('1vq8', 5)),(('4y4o', 25),('3cc2', 2)),(('4ybb', 1),('1vq8', 5)),(('4ybb', 1),('3cc2', 2)), (('6hma', 1),('1vq8', 5)),(('6hma', 1),('3cc2', 2))])
        
        # 23S
        liste_rouge = [(('1vq8', 18), ('5dm6', 5)),(('5dm6', 5), ('3ccl', 14)),(('5dm6', 5), ('3cce', 12)),(('1vq8', 1), ('6eri', 2)), (('5wfs', 19), ('4u4r', 6)),(('4ybb', 12), ('4y4o', 35)),(('4ybb', 12), ('4wsd', 40)),(('4y4o', 34), ('2qex', 12)),(('4y4o', 28), ('5wfs', 6)),(('4v67', 43), ('2qex', 12)), (('4y4o', 50), ('6ek0', 4)), (('5dm6', 4), ('4ybb', 17)), (('4w2g', 52), ('6ek0', 10)), (('4y4o', 47), ('1vq8', 19)),(('4y4o', 47), ('1yhq', 24)),(('4y4o', 47), ('2qex', 19)),(('4y4o', 47), ('4u4r', 13)),(('4ybb', 13), ('1vqo', 8)),(('4w2f', 23), ('4u4r', 36)),(('1vqo', 8), ('5dm6', 1)),(('1vqo', 8), ('1mms', 1)),(('1vq8', 19), ('6eri', 1)),(('1yhq', 24), ('6eri', 1)),(('6eri', 1), ('2qex', 19)), (('5dm6', 14), ('4y4o', 8)),(('4ybb', 17), ('4y4o', 8)), (('4ybb', 30), ('2zjr', 15)),(('2zjr', 15), ('5wfs', 11)), (('2zjr', 15), ('6eri', 8)), (('4y4o', 12), ('2zjr', 15)), (('4w2g', 52), ('4y4o', 28)), (('4y4o', 28), ('4v67', 7))]
        #liste_rouge.extend([(('4ybb', 12), ('4y4o', 38)),(('4ybb', 12), ('4ybb', 21)),(('4ybb', 12), ('6ek0', 7)),(('4v67', 7), ('4y4o', 38)),(('4v67', 7), ('4ybb', 21)),(('4v67', 7), ('6ek0', 7)),(('5wfs', 19), ('4y4o', 38)),(('5wfs', 19), ('4ybb', 21)),(('5wfs', 19), ('6ek0', 7))])
        ## non homologues meme type rmsd 2-2.5, sim 0.7-0.75
        liste_rouge.extend([(('4y4o', 23), ('6h4n', 24)),(('2zjr', 3), ('4u4r', 18)),(('4u27', 54), ('4y4o', 38)),(('4u27', 54), ('6ek0', 7)),(('4u27', 54), ('4u3u', 7)), (('4y4o', 23), ('4u27', 3)),(('4y4o', 23), ('5afi', 24))])
        
        ## les non homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim à 0.75 mais pas avec la rmsd à 2
        liste_rouge.extend([(('4ybb', 12), ('6ek0', 10)),(('4w2g', 52), ('6ek0', 10)),(('2zjr', 3),('6ek0', 10)),(('2zjr', 3),('4u4r', 18)),(('4w2g', 52),('4u4r', 18)),(('4v67', 7), ('4u4r', 18)),(('5wfs', 19), ('4u4r', 18)),(('4v67', 7),('6ek0', 10)),(('5wfs', 19),('6ek0', 10)),(('4y4o', 23),('6h4n', 24)),(('4y4o', 23), ('4u27', 3)),(('4y4o', 23),('5afi', 24))])
        #liste_rouge = []
    
        distrib_sim = []
        distrib_rmsd = []
        distrib_sim_rouge = []
        distrib_sim_verte = []
        distrib_rmsd_rouge = []
        distrib_rmsd_verte = []
        distrib_sim_pos = []
        distrib_rmsd_pos = []
        
        compteur = 0
        compter = 0
        compteer = 0
        dico_type = {}
        for elt in liste_representant :
            dico_type.update({elt : renvoi_num_ARN(elt)})
        
        print("combien d'elts ?")
        print(len(liste_representant))
        for i in range(len(liste_representant)) :
            for j in range(i+1, len(liste_representant)) :
                if liste_representant[i][0] != '6hrm' and  liste_representant[j][0] != '6hrm' :
                    ok = True
                    for groupe in liste_homologues :
                        if liste_representant[i] in groupe and liste_representant[j] in groupe :
                            ok = False
                    if ok :
#                         type1 = renvoi_num_ARN(liste_representant[i])
#                         type2 = renvoi_num_ARN(liste_representant[j])
                        type1 = dico_type[liste_representant[i]]
                        type2 = dico_type[liste_representant[j]]
                        diff_type = False
                        if (type1 != type2 and type1 not in ["23S","28S", "25S", "16S", "18S"] and type2 not in ["23S","28S", "25S", "16S", "18S"]) or ((type1 in ["23S","28S", "25S"] and type2 not in ["23S","28S","25S"]) or (type2 in ["23S","28S", "25S"] and type1 not in ["23S","28S","25S"]) ) or ((type1 in ["16S","18S"] and type2 not in  ["16S","18S"]) or (type2 in ["16S", "18S"] and type1 not in ["16S", "18S"])) :
                            diff_type = True
                            print(type1)
                            print(type2)
                        nom1 = "fichier_%s_%d_taille_0.pdb"%(liste_representant[i][0], liste_representant[i][1])
                        nom2 = "fichier_%s_%d_taille_0.pdb"%(liste_representant[j][0], liste_representant[j][1])
                        compteur += 1
                        
                        if (nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)] <= seuil_rmsd:
                            print("raaaa")
                            #exit()
                            if (liste_representant[i], liste_representant[j]) in dico_pos.keys() or (liste_representant[j], liste_representant[i]) in dico_pos.keys() :
                                
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                    
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_sim_pos.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    distrib_rmsd_pos.append(rmsd[(nom1,nom2)])
                                elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_sim_pos.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd_pos.append(rmsd[(nom1,nom2)])
                            elif (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_sim_verte.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    distrib_rmsd_verte.append(rmsd[(nom1,nom2)])
                                elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_sim_verte.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd_verte.append(rmsd[(nom1,nom2)])
                            
                            elif (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge  or  diff_type:
                                
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                    if dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] > 0.7 : 
                                        print("gros tas")
                                        print((liste_representant[i], liste_representant[j]))
                                    if (liste_representant[i], liste_representant[j]) in dico_score.keys() :
                                        csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])
                                    else :
                                        csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[j], liste_representant[i])] ])
                                        
                                    distrib_sim_rouge.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    distrib_rmsd_rouge.append(rmsd[(nom1,nom2)])
                                elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                    if dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] > 0.7 : 
                                        print("gros tas")
                                        print((liste_representant[i], liste_representant[j]))
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_sim_rouge.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd_rouge.append(rmsd[(nom1,nom2)])
                            

                            else :
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
#                                     if rmsd[(nom1, nom2)] <= 4 and dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)] ])
                                    distrib_sim.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_rmsd.append(rmsd[(nom1,nom2)])
                                elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                    #if rmsd[(nom1, nom2)] <= 4 and dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                    #    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)] ])
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_sim.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd.append(rmsd[(nom1,nom2)])
                        elif (nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)] <= seuil_rmsd :
                            print("riiii")
                            #exit()
                            if (liste_representant[i], liste_representant[j]) in dico_pos.keys() or (liste_representant[j], liste_representant[i]) in dico_pos.keys() :
                                
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                    distrib_sim_pos.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    distrib_rmsd_pos.append(rmsd[(nom2,nom1)])
                                elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim:
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])
                                    distrib_sim_pos.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd_pos.append(rmsd[(nom2,nom1)])
                            elif (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])
                                    distrib_sim_verte.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    distrib_rmsd_verte.append(rmsd[(nom2,nom1)])
                                elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])
                                    distrib_sim_verte.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd_verte.append(rmsd[(nom2,nom1)])
                            
                            elif (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge or diff_type :
                                
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])                                   
                                    distrib_sim_rouge.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    distrib_rmsd_rouge.append(rmsd[(nom2,nom1)])
                                elif  (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])                                    
                                    distrib_sim_rouge.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd_rouge.append(rmsd[(nom2,nom1)])
                            
                            else :                            
                                
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
#                                     if rmsd[(nom2, nom1)] <= 4 and dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)] ])
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])
                                    distrib_sim.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                    distrib_rmsd.append(rmsd[(nom2,nom1)])
                                elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
#                                     if rmsd[(nom2, nom1)] <= 4 and dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)] ])
                                    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])
                                    distrib_sim.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                    distrib_rmsd.append(rmsd[(nom2,nom1)])
        #                 else : 
        #                     print(nom1, nom2)
                    else :
                        compter += 1
                    print(compteur)
                else :
                    compteer += 1
                    
    print(compter)
    print(compteur)
    print(compteer)
    print(len(distrib_sim))
    print(len(distrib_rmsd))
    
    for i in range(len(distrib_sim)) :
        if distrib_sim[i] == 0.6 : 
            print(distrib_sim[i])
            print(distrib_rmsd[i])
    plt.figure(figsize=(12,6))
    axs = plt.gca()
    
    plt.scatter(distrib_sim_rouge, distrib_rmsd_rouge, c='red')
    plt.scatter(distrib_sim, distrib_rmsd)
    plt.scatter(distrib_sim_pos, distrib_rmsd_pos)
    plt.scatter(distrib_sim_verte, distrib_rmsd_verte, c='green')
    
    axs.set_xlabel("Similarite")
    axs.set_ylabel("RMSD (en Angstrom)")
    #axs.set_xticks(np.arange(0,1.1,0.1))
    #axs.set_xticks(np.arange(0, max(310, max(dico_min_tot.values())+20), 10))
    axs.set_title("Distribution des valeurs de RMSD en fonction des valeurs de similarite \n")
    #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_%s_seuil_pos_rmsd_%d.png"%(num_ARN, min_ou_max, seuil_hom))
    plt.show()
    plt.close()


''' 03/02/20 '''
def distrib_sim_rmsd_inter_types(num_ARN1, num_ARN2):
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
        
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_sim :
        mon_depickler = pickle.Unpickler(fichier_sim)
        dico_sim = mon_depickler.load()
        
    liste_cles_1 = {}
    if isinstance(num_ARN1, list) :
        
        liste_representant_1 = []
        for elt in num_ARN1 :
            liste_cles_1.update({elt : []})
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr = mon_depickler_1.load()
                
                liste_representant_1.extend(liste_repr)
                liste_cles_1[elt].extend(liste_repr)
        
    else :
        liste_cles_1.update({num_ARN1 : []})
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN1, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant_1 = mon_depickler_1.load()
                liste_cles_1[num_ARN1].extend(liste_representant_1)
                
    liste_cles_2 = {}
    if isinstance(num_ARN2, list) :
        
        liste_representant_2 = []
        for elt in num_ARN2 :
            liste_cles_2.update({elt : []})
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr = mon_depickler_1.load()
                
                liste_representant_2.extend(liste_repr)
                liste_cles_2[elt].extend(liste_repr)
        
    else :
        liste_cles_2.update({num_ARN2 : []})
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN2, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant_2 = mon_depickler_1.load()
                liste_cles_2[num_ARN2].extend(liste_representant_2)
    
    
    distrib_sim = []
    distrib_rmsd = []
    
    compteur = 0
    for i in range(len(liste_representant_1)) :
        for j in range(len(liste_representant_2)) :
            if '6hrm' not in liste_representant_1[i][0] and '6hrm' not in liste_representant_2[j][0] :
                nom1 = "fichier_%s_%d_taille_4.pdb"%(liste_representant_1[i][0], liste_representant_1[i][1])
                nom2 = "fichier_%s_%d_taille_4.pdb"%(liste_representant_2[j][0], liste_representant_2[j][1])
                compteur += 1
                if (nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None :
                    distrib_rmsd.append(rmsd[(nom1,nom2)])
                    if (liste_representant_1[i], liste_representant_2[j]) in dico_sim.keys() :
                        distrib_sim.append(dico_sim[ (liste_representant_1[i], liste_representant_2[j])]["sim"])
                    else :
                        distrib_sim.append(dico_sim[ (liste_representant_2[j], liste_representant_1[i])]["sim"])
                    
                elif (nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None :
                    
                    distrib_rmsd.append(rmsd[(nom2,nom1)])
                    if (liste_representant_1[i], liste_representant_2[j]) in dico_sim.keys() :
                        distrib_sim.append(dico_sim[ (liste_representant_1[i], liste_representant_2[j])]["sim"])
                    else :
                        distrib_sim.append(dico_sim[ (liste_representant_2[j], liste_representant_1[i])]["sim"])
                else : 
                    print(nom1, nom2)
    print(compteur)
    print(len(distrib_sim))
    print(len(distrib_rmsd))
    
    for i in range(len(distrib_sim)) :
        if distrib_sim[i] == 0.6 : 
            print(distrib_sim[i])
            print(distrib_rmsd[i])
    plt.figure(figsize=(12,6))
    axs = plt.gca()
    plt.scatter(distrib_sim, distrib_rmsd)
    axs.set_xlabel("Similarite")
    axs.set_ylabel("RMSD (en Angstrom)")
    #axs.set_xticks(np.arange(0, max(310, max(dico_min_tot.values())+20), 10))
    axs.set_title("Distribution des valeurs de RMSD en fonction des valeurs de similarite \n")
    axs.set_xticks(np.arange(0,1.1,0.1))
    #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_%s_seuil_pos_rmsd_%d.png"%(num_ARN, min_ou_max, seuil_hom))
    plt.show()
    plt.close()
    
''' 03/02/20 '''
def distrib_sim_rmsd_plusieurs_types(liste_num_ARN):
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
        
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_sim :
        mon_depickler = pickle.Unpickler(fichier_sim)
        dico_sim = mon_depickler.load()
    
    liste_representant_1 = []
    compteur = 0
    for num_ARN in liste_num_ARN :
        if len(liste_representant_1) < compteur+1 :
            liste_representant_1.append([])
        if isinstance(num_ARN, list) :
            for elt in num_ARN :
                with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_repr = mon_depickler_1.load()
                    
                    liste_representant_1[compteur].extend(liste_repr)
            
        else :
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_representant_1[compteur] = mon_depickler_1.load()               
        compteur += 1   
    
    
    distrib_sim = []
    distrib_rmsd = []
    
    compteur = 0
    for k in range(len(liste_representant_1)) :
        for l in range(k+1, len(liste_representant_1)) :
            for i in range(len(liste_representant_1[k])) :
                    for j in range(len(liste_representant_1[l])) :
                            if liste_representant_1[k][i][0] != '6hrm' and  liste_representant_1[l][j][0] != '6hrm'  and not (liste_representant_1[k][i] == ('2r8s', 1) and liste_representant_1[l][j] == ('6d8o', 2)) and not (liste_representant_1[k][i] == ('1gid', 1) and liste_representant_1[l][j] == ('6d8o', 2)) and not (liste_representant_1[k][i] == ('1l8v', 1) and liste_representant_1[l][j] == ('6d8o', 2)):
                                print(liste_representant_1[k][i])
                                nom1 = "fichier_%s_%d_taille_4.pdb"%(liste_representant_1[k][i][0], liste_representant_1[k][i][1])
                                nom2 = "fichier_%s_%d_taille_4.pdb"%(liste_representant_1[l][j][0], liste_representant_1[l][j][1])
                                compteur += 1
                                if (nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None :
                                    distrib_rmsd.append(rmsd[(nom1,nom2)])
                                    if (liste_representant_1[k][i], liste_representant_1[l][j]) in dico_sim.keys() :
                                        distrib_sim.append(dico_sim[ (liste_representant_1[k][i], liste_representant_1[l][j])]["sim"])
                                    else :
                                        distrib_sim.append(dico_sim[ (liste_representant_1[l][j], liste_representant_1[k][i])]["sim"])
                                    
                                elif (nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None :
                                    
                                    distrib_rmsd.append(rmsd[(nom2,nom1)])
                                    if (liste_representant_1[k][i], liste_representant_1[l][j]) in dico_sim.keys() :
                                        distrib_sim.append(dico_sim[ (liste_representant_1[k][i], liste_representant_1[l][j])]["sim"])
                                    else :
                                        distrib_sim.append(dico_sim[ (liste_representant_1[l][j], liste_representant_1[k][i])]["sim"])
                                else : 
                                    print(nom1, nom2)
    print(compteur)
    print(len(distrib_sim))
    print(len(distrib_rmsd))
    
    for i in range(len(distrib_sim)) :
        if distrib_sim[i] == 0.6 : 
            print(distrib_sim[i])
            print(distrib_rmsd[i])
    plt.figure(figsize=(12,6))
    axs = plt.gca()
    plt.scatter(distrib_sim, distrib_rmsd)
    axs.set_xlabel("Similarite")
    axs.set_ylabel("RMSD (en Angstrom)")
    #axs.set_xticks(np.arange(0, max(310, max(dico_min_tot.values())+20), 10))
    axs.set_title("Distribution des valeurs de RMSD en fonction des valeurs de similarite \n")
    axs.set_xticks(np.arange(0,1.1,0.1))
    #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_%s_seuil_pos_rmsd_%d.png"%(num_ARN, min_ou_max, seuil_hom))
    plt.show()
    plt.close()


''' 06/02/20 '''
def distrib_3D_sim_rmsd_score(num_ARN, liste_homologues):
    
    for groupe in liste_homologues :
        if ('2r8s', 1) in groupe :
            groupe.append(('6d8o', 2))
    
    liste_verte = [(('1vqo', 12), ('6nd6', 7)),(('4y4o', 40), ('4u4r', 8)),(('6eri', 9), ('4u4r', 8)), (('6hma', 15), ('1vq8', 15)), (('5dm6', 2), ('6eri', 1)), (('4u27', 5), ('1vqp', 19)),(('4u27', 5), ('3ccr', 4)),(('4u27', 5), ('3cc7', 4)),(('5dm6', 7), ('1vq8', 11)),(('5dm6', 7), ('3ccr', 12)),(('5dm6', 7), ('3ccl', 11)),(('5dm6', 7), ('3ccs', 12)),(('5dm6', 7), ('3cc7', 12)),(('5dm6', 7), ('3ccq', 12)),(('5dm6', 7), ('3cce', 9)),(('5dm6', 7), ('3ccu', 11)),(('4ybb', 8), ('1vq8', 11)),(('4ybb', 8), ('3ccr', 12)),(('4ybb', 8), ('3ccl', 11)),(('4ybb', 8), ('3ccs', 12)),(('4ybb', 8), ('3cc7', 12)),(('4ybb', 8), ('3ccq', 12)),(('4ybb', 8), ('3cce', 9)),(('4ybb', 8), ('3ccu', 11)),(('3cc2', 2), ('6hma', 1)),(('5dm6', 13), ('1vqp', 19)),(('5dm6', 13), ('3ccr', 4)),(('5dm6', 13), ('3cc7', 4)), (('6hma', 9), ('4u4r', 8)), (('4y4o', 25), ('1vq8', 5)),(('4y4o', 25), ('3cc2', 2)),(('6hma', 15), ('4u4r', 16)), (('1vqp', 19), ('6eri', 9)), (('3ccr', 4), ('6eri', 9)),(('6eri', 9), ('3cc7', 4)),(('5afi', 15), ('4u4r', 28)),(('6hma', 9), ('1vqp', 19)),(('6hma', 9), ('3ccr', 4)),(('6hma', 9), ('3cc7', 4)),(('1vq8', 11), ('6hma', 6)),(('3ccr', 12), ('6hma', 6)),(('6hma', 6), ('3ccl', 11)),(('6hma', 6), ('3ccs', 12)),(('6hma', 6), ('3cc7', 12)),(('6hma', 6), ('3ccq', 12)),(('6hma', 6), ('3cce', 9)),(('6hma', 6), ('3ccu', 11)),(('4y4o', 40), ('1vqp', 19)),(('4y4o', 40), ('3ccr', 4)),(('4y4o', 40), ('3cc7', 4)),(('4y4o', 3), ('3ccl', 11)),(('4y4o', 3), ('3cce', 9)),(('4v51', 22), ('3ccl', 11)),(('4v51', 22), ('3cce', 9)),(('5e81', 12), ('3ccl', 11)),(('5e81', 12), ('3cce', 9)),(('3ccl', 11), ('4v90', 13)),(('4v90', 13), ('3cce', 9)),(('1vqo', 12), ('4u4r', 15)),(('4ybb', 1), ('1vq8', 5)),(('4y4o', 3), ('3ccr', 12)),(('4y4o', 3), ('3ccs', 12)),(('4y4o', 3), ('3cc7', 12)),(('4y4o', 3), ('3ccu', 11)),(('4v51', 22), ('3ccr', 12)),(('4v51', 22), ('3ccs', 12)),(('4v51', 22), ('3cc7', 12)),(('4v51', 22), ('3ccu', 11)),(('5e81', 12), ('3ccr', 12)),(('5e81', 12), ('3ccs', 12)),(('5e81', 12), ('3cc7', 12)),(('5e81', 12), ('3ccu', 11)),(('3ccr', 12), ('4v90', 13)),(('4v90', 13), ('3ccs', 12)),(('4v90', 13), ('3cc7', 12)),(('4v90', 13), ('3ccu', 11)),(('4y4o', 3), ('3ccq', 12)),(('4v51', 22), ('3ccq', 12)),(('5e81', 12), ('3ccq', 12)),(('4v90', 13), ('3ccq', 12)),(('4ybb', 1), ('3cc2', 2)),(('4y4o', 3), ('1vq8', 11)),(('4v51', 22), ('1vq8', 11)),(('5e81', 12), ('1vq8', 11)),(('1vq8', 5), ('6hma', 1)),(('1vq8', 11), ('4v90', 13)), (('5dm6', 9), ('1vq8', 21)), (('6hma', 2), ('4u4r', 19)), (('4ybb', 6), ('1yhq', 24)), (('6hma', 15), ('3cc2', 6))]
    liste_verte.extend([(('5nwy', 17), ('4y4o', 4)), (('6ek0', 2), ('6az1', 1)), (('4u27', 54), ('4y4o', 22)), (('4ybb', 35), ('4u3u', 22))])

    
    #liste_verte = []
    liste_rouge = [(('1vq8', 18), ('5dm6', 5)),(('5dm6', 5), ('3ccl', 14)),(('5dm6', 5), ('3cce', 12)),(('1vq8', 1), ('6eri', 2)), (('5wfs', 19), ('4u4r', 6)),(('4ybb', 12), ('4y4o', 35)),(('4ybb', 12), ('4wsd', 40)),(('4y4o', 34), ('2qex', 12)),(('4y4o', 28), ('5wfs', 6)),(('4v67', 43), ('2qex', 12)), (('4y4o', 50), ('6ek0', 4)), (('5dm6', 4), ('4ybb', 17)), (('4w2g', 52), ('6ek0', 10)), (('4y4o', 47), ('1vq8', 19)),(('4y4o', 47), ('1yhq', 24)),(('4y4o', 47), ('2qex', 19)),(('4y4o', 47), ('4u4r', 13)),(('4ybb', 13), ('1vqo', 8)),(('4w2f', 23), ('4u4r', 36)),(('1vqo', 8), ('5dm6', 1)),(('1vqo', 8), ('1mms', 1)),(('1vq8', 19), ('6eri', 1)),(('1yhq', 24), ('6eri', 1)),(('6eri', 1), ('2qex', 19)), (('5dm6', 14), ('4y4o', 8)),(('4ybb', 17), ('4y4o', 8)), (('4ybb', 30), ('2zjr', 15)),(('2zjr', 15), ('5wfs', 11)), (('2zjr', 15), ('6eri', 8)), (('4y4o', 12), ('2zjr', 15)), (('4w2g', 52), ('4y4o', 28)), (('4y4o', 28), ('4v67', 7)), (('4y4o', 29), ('4y4o', 8)),(('1vqo', 8), ('1vqo', 14)),(('1vqo', 8), ('3cc2', 1)),(('4u4r', 12), ('4u4r', 36))]
    #liste_rouge = []
    
    
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_sim :
        mon_depickler = pickle.Unpickler(fichier_sim)
        dico_sim = mon_depickler.load()
        
    if isinstance(num_ARN, list) :
        print(num_ARN)
        liste_representant = []
        dico_pos = {}
        dico_min_tot_rmsd = {}
        dico_min_tot = {}
        if num_ARN in [['23S', '25S', '28S'], ['16S', '18S']] :
            print("rapala")
            with open("dico_min_pos_score_%s_seuil_30.pickle"%(num_ARN), "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos.update(mon_depickler.load())
                
        with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "rb") as fichier_rmsd :
            mon_depickler = pickle.Unpickler(fichier_rmsd)
            dico_min_tot_rmsd.update(mon_depickler.load())   
                
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
            mon_depickler = pickle.Unpickler(fichier_score)
            dico_min_tot.update(mon_depickler.load())
            
        for elt in num_ARN :
            
            #liste_cles.update({elt : []})
            if isinstance(elt, list) :
                print(elt)
                if elt in [['23S', '25S', '28S'], ['16S', '18S']] :
                    print("rapala")
                    with open("dico_min_pos_score_%s_seuil_30.pickle"%(elt), "rb") as fichier_score_pos :
                        mon_depickler = pickle.Unpickler(fichier_score_pos)
                        dico_pos.update(mon_depickler.load())
                        
                with open("dico_min_tot_rmsd_%s.pickle"%elt, "rb") as fichier_rmsd :
                    mon_depickler = pickle.Unpickler(fichier_rmsd)
                    dico_min_tot_rmsd.update(mon_depickler.load())   
                        
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_min_tot.update(mon_depickler.load())
                
                
                for e in elt :
                    print(e)
                    with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%e, 'rb') as fichier_num_arn :
                        mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                        liste_repr = mon_depickler_1.load()
                        
                        liste_representant.extend(liste_repr)
            else :
                print(elt)
                with open("dico_min_tot_rmsd_%s.pickle"%elt, "rb") as fichier_rmsd :
                    mon_depickler = pickle.Unpickler(fichier_rmsd)
                    dico_min_tot_rmsd.update(mon_depickler.load())   
                        
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_min_tot.update(mon_depickler.load())
                
                with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_representant.extend(mon_depickler_1.load())
                #liste_cles[elt].extend(liste_repr)
        
    else :
        #liste_cles.update({num_ARN : []})
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant = mon_depickler_1.load()
                #liste_cles[num_ARN].extend(liste_representant)
                
        with open("dico_min_pos_score_%s_seuil_30.pickle"%num_ARN, "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos = mon_depickler.load()
                
        with open("dico_min_tot_rmsd_%s.pickle"%num_ARN, "rb") as fichier_rmsd :
                mon_depickler = pickle.Unpickler(fichier_rmsd)
                dico_min_tot_rmsd = mon_depickler.load()  
                
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
                mon_depickler = pickle.Unpickler(fichier_score)
                dico_min_tot = mon_depickler.load()
    
    with open("fichier_forte_sim_faible_rmsd_non_homol.pickle", "w", newline="") as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(["Paire","Type ARN 1", "Type ARN 2", "Sim", "RMSD", "Score d'alignement"])
    

        print(len(dico_min_tot_rmsd))
    
            
        distrib_score = []
        distrib_rmsd = []
        distrib_sim = []
            
        distrib_score_pos = []
        distrib_rmsd_pos = []
        distrib_sim_pos = []
            
        distrib_score_vert = []
        distrib_rmsd_vert = []
        distrib_sim_vert = []
            
        distrib_score_rouge = []
        distrib_rmsd_rouge = []
        distrib_sim_rouge = []
        
        dico_type = {}
        for elt in liste_representant :
            dico_type.update({elt : renvoi_num_ARN(elt)})
            
            
        for i in range(len(liste_representant)) :
                for j in range(i+1, len(liste_representant)) :
                    ok = True
                    for groupe in liste_homologues :
                        if liste_representant[i] in groupe and liste_representant[j] in groupe :
                            ok = False
                    if ok and liste_representant[i][0] != '6hrm' and liste_representant[j][0] != '6hrm' :
                        diff_type = False
#                         type1 = renvoi_num_ARN(liste_representant[i])
#                         type2 = renvoi_num_ARN(liste_representant[j])
                        type1 = dico_type[liste_representant[i]]
                        type2 = dico_type[liste_representant[j]]
                        diff_type = False
                        if (type1 != type2 and type1 not in ["23S","28S", "25S", "16S", "18S"] and type2 not in ["23S","28S", "25S", "16S", "18S"]) or ((type1 in ["23S","28S", "25S"] and type2 not in ["23S","28S","25S"]) or (type2 in ["23S","28S", "25S"] and type1 not in ["23S","28S","25S"]) ) or ((type1 in ["16S","18S"] and type2 not in  ["16S","18S"]) or (type2 in ["16S", "18S"] and type1 not in ["16S", "18S"])) :
                            diff_type = True
                            print(type1)
                            print(type2)
                        
                        if dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])] != None :
                            if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
                                if (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                    distrib_sim_vert.append(1-dico_sim[(liste_representant[i], liste_representant[j])]["sim"])
                                elif (liste_representant[i], liste_representant[j]) in dico_pos.keys() or (liste_representant[j], liste_representant[i]) in dico_pos.keys() :
                                    distrib_sim_pos.append(1-dico_sim[(liste_representant[i], liste_representant[j])]["sim"])
                                elif (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge or diff_type :
                                    if dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])] <= 4 :
                                        csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[(liste_representant[i], liste_representant[j])]["sim"], dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])], dico_min_tot[(liste_representant[i], liste_representant[j])] ])
                                    distrib_sim_rouge.append(1-dico_sim[(liste_representant[i], liste_representant[j])]["sim"])
                                else :
                                    distrib_sim.append(1-dico_sim[(liste_representant[i], liste_representant[j])]["sim"])
                            else :
                                if (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                    distrib_sim_vert.append(1-dico_sim[(liste_representant[j], liste_representant[i])]["sim"])
                                elif (liste_representant[i], liste_representant[j]) in dico_pos.keys() or (liste_representant[j], liste_representant[i]) in dico_pos.keys() :
                                    distrib_sim_pos.append(1-dico_sim[(liste_representant[j], liste_representant[i])]["sim"])
                                elif (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge or diff_type:
                                    if dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])] <= 4 :
                                        csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[(liste_representant[j], liste_representant[i])]["sim"], dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])], dico_min_tot[(liste_representant[i], liste_representant[j])]])
                                    distrib_sim_rouge.append(1-dico_sim[(liste_representant[j], liste_representant[i])]["sim"])
                                else :
                                    distrib_sim.append(1-dico_sim[(liste_representant[j], liste_representant[i])]["sim"])
        
                                
                            if dico_min_tot[(liste_representant[i], liste_representant[j])] >= 50 and dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])] <= 4 :
                                print(liste_representant[i], liste_representant[j])
                                print(dico_min_tot[(liste_representant[i], liste_representant[j])])
                                print(dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])])
                            
                            if (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                distrib_rmsd_vert.append(dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])])
                                distrib_score_vert.append(-dico_min_tot[(liste_representant[i], liste_representant[j])])
                            elif (liste_representant[i], liste_representant[j]) in dico_pos.keys() or (liste_representant[j], liste_representant[i]) in dico_pos.keys() :
                                distrib_rmsd_pos.append(dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])])
                                distrib_score_pos.append(-dico_min_tot[(liste_representant[i], liste_representant[j])])
                            elif (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge or diff_type:
                                distrib_rmsd_rouge.append(dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])])
                                distrib_score_rouge.append(-dico_min_tot[(liste_representant[i], liste_representant[j])])
                            else :
                                distrib_rmsd.append(dico_min_tot_rmsd[(liste_representant[i], liste_representant[j])])
                                distrib_score.append(-dico_min_tot[(liste_representant[i], liste_representant[j])])
                        
                        
        #print(distrib_rmsd)
        #print(distrib_score)
        #print(distrib_sim)
        print(distrib_rmsd_vert)
        print(distrib_score_vert)
        print(distrib_sim_vert)
               
        fig = plt.figure()
        ax = Axes3D(fig)
        
        ax.scatter(distrib_sim_rouge, distrib_score_rouge, distrib_rmsd_rouge, c='red')
        ax.scatter(distrib_sim, distrib_score, distrib_rmsd)
        ax.scatter(distrib_sim_pos, distrib_score_pos, distrib_rmsd_pos, c='orange')
        ax.scatter(distrib_sim_vert, distrib_score_vert, distrib_rmsd_vert, c='green')
        
        #ax.set_yticklabels(np.arange(175, -5, -25))
        #ax.set_xticklabels([1.0,0.8,0.6,0.4,0.2,0])
        ax.set_zlabel("RMSD")
        ax.set_ylabel("Score d'alignement global")
        ax.set_xlabel("Similarité")
        plt.show()
        

''' 18/12/19
cree les fichiers de sequences concatenees de toutes les occurrences d'A-minor du type d'ARN passe en parametre
pour etre ensuite utilisees par le programme CDHIT '''
def creation_fichier_input_cdhit(num_ARN):
    if isinstance(num_ARN, list) :
        liste_representant = []
        for elt in num_ARN :
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr = mon_depickler_1.load()
                
                liste_representant.extend(liste_repr)
    else :
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
            mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
            liste_representant = mon_depickler_1.load()
    
    with open("Nouvelles_donnees/fichier_fasta_seq_concatenees_%s.fa"%num_ARN, 'w') as fichier_fasta :
        for elt in liste_representant :
            seqs = concatenate_sequences(elt)
            
            fichier_fasta.write(">%s_%d\n"%(elt[0], elt[1]))
            fichier_fasta.write(seqs+"\n")
    
''' 18/12/19 
recupere les clusters obtenus par l'algo CDHIT (version a choisir)
pour un type d'ARN passe en parametre
'''
def traitement_res_cdhit(num_ARN):
    with open("/home/coline/cdhit/res_%s_80_gap.clstr"%num_ARN, 'r') as fichier_fasta :
        clusters = []
        ligne = fichier_fasta.readline()
        while ligne != "" :
            if "Cluster" in ligne : 
                clusters.append([])
                ligne = fichier_fasta.readline()
                while ligne != "" and "Cluster" not in ligne :
                    #print(ligne)
                    elt = ligne.split(" ")[1][1:len(ligne.split(" ")[1])-3]
                    print(elt)
                    clusters[len(clusters)-1].append((elt.split("_")[0], int(elt.split("_")[1])))
                    ligne = fichier_fasta.readline()
                
            else :
                print("bizarre")
                exit(0)
    print(clusters)
    return clusters

    
''' 07/01/20 
calculer un clustering en deux groupes des paires d'occurrences 
en fonction du score d'alignement global des sequences et de la RMSD (calcul de la distance euclidienne entre toutes les paires 
mises dans un graphique de la RMSD en fonction du score d'alignement)
en parametres : les valeurs de RMSD et de score d'alignement global de toutes les paires 
tres long a executer
'''
def calcul_kmeans_sur_distrib(dico_min_tot_scores, dico_min_tot_rmsd):
    matrice_distance_euclidienne = [[0] * len(dico_min_tot_scores) for _ in range(len(dico_min_tot_scores))]
    
    for i in range(len(dico_min_tot_scores)) :
        for j in range(i, len(dico_min_tot_scores)) :
            print(i)
            print(j)
            if i == j :
                matrice_distance_euclidienne[i][j] = 0
            else :
                if list(dico_min_tot_scores.keys())[i] == list(dico_min_tot_rmsd.keys())[i] and  list(dico_min_tot_scores.keys())[j] == list(dico_min_tot_rmsd.keys())[j] :
                    cle1 = list(dico_min_tot_scores.keys())[i]
                    cle2 = list(dico_min_tot_scores.keys())[j]
                    if dico_min_tot_rmsd[cle1] != None and dico_min_tot_rmsd[cle2] != None :
                        matrice_distance_euclidienne[i][j] = sqrt((dico_min_tot_rmsd[cle1] - dico_min_tot_rmsd[cle2])**2 + (dico_min_tot_scores[cle1] - dico_min_tot_scores[cle2])**2)
                        matrice_distance_euclidienne[j][i] = sqrt((dico_min_tot_rmsd[cle1] - dico_min_tot_rmsd[cle2])**2 + (dico_min_tot_scores[cle1] - dico_min_tot_scores[cle2])**2)
                        #print(matrice_distance_euclidienne[i][j])
                else :
                    print("pas pratique")
    #print(matrice_distance_euclidienne)
    
    with open("matrice_distance_euclidienne_score_rmsd_23S_25S_28S.pickle", 'wb') as fichier_pickle :
        mon_pickler = pickle.Pickler(fichier_pickle)
        mon_pickler.dump(matrice_distance_euclidienne)
    
    clustering = KMeans(n_clusters=2, n_init=50).fit(matrice_distance_euclidienne)  # init=init, n_init=1).fit(matrice)
    tab_clustering = [['']]
    print(len(tab_clustering))
    print(len(clustering.labels_))
    print(clustering.labels_)
    compteur = 0
    dico_groupe_1_scores = {}
    dico_groupe_1_rmsd = {}
    dico_groupe_2_scores = {}
    dico_groupe_2_rmsd = {}
    for elt in clustering.labels_ :
        if elt == 1 :
            dico_groupe_1_scores.update({list(dico_min_tot_scores.keys())[compteur] : dico_min_tot_scores[list(dico_min_tot_scores.keys())[compteur]] })
            dico_groupe_1_rmsd.update({list(dico_min_tot_scores.keys())[compteur] : dico_min_tot_rmsd[list(dico_min_tot_scores.keys())[compteur]] })
        elif elt == 2 :
            dico_groupe_2_scores.update({list(dico_min_tot_scores.keys())[compteur] : dico_min_tot_scores[list(dico_min_tot_scores.keys())[compteur]] })
            dico_groupe_2_rmsd.update({list(dico_min_tot_scores.keys())[compteur] : dico_min_tot_rmsd[list(dico_min_tot_scores.keys())[compteur]] })
        compteur += 1
        
    plt.scatter(dico_groupe_1_scores.values(), dico_groupe_1_rmsd.values())
    plt.scatter(dico_groupe_2_scores.values(), dico_groupe_2_rmsd.values())
    plt.show()
    

### Filtrage resolution ###
''' pour chaque structure PDB possedant un ou plusieurs motif A-minor 
recupere la resolution dans le fichier PDB
stocke toutes les resolutions dans un dictionnaire '''
def stocker_resolution():
    with open("resolutions.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        resolutions = mon_depickler.load()
        
        compteur = 0
        for fic in os.listdir("Graphs") :
            if fic.split(".")[0] not in resolutions.keys() :
                print(fic)
                print(fic.split(".")[0])
                print(compteur)
                res = get_resolution(fic.split(".")[0])
                resolutions.update({fic.split(".")[0] : res})
                compteur += 1
            
    with open("resolutions.pickle", 'wb') as fichier_pickle :
            mon_pickler = pickle.Pickler(fichier_pickle)
            mon_pickler.dump(resolutions)



''' pour les structures PDB ou la premiere methode de recuperation de la resolution (en utilisant le champ exp_header_0_diffraction_resolution du fichier PDB)
n'a pas marche, on utilise une deuxieme methode (champ _em_3d_reconstruction. de la structure PDB)
on stocke toutes les resolutions dans un dictionnaire  '''
def stocker_resolution_qui_manque():
    with open("resolutions.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        resolutions = mon_depickler.load()
    
        for cle in resolutions.keys() :
            if resolutions[cle] == None :
                if cle.upper()+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
                    print("petit rat")
                    url = 'http://www.rcsb.org/pdb/files/%s.cif' % cle.upper()
                    urllib.request.urlretrieve(url, PATH_MMCIF+cle.upper()+".cif")
                    
                try:
                    doc = cif.read_file(PATH_MMCIF+cle.upper()+".cif")
                    
                    block = doc.sole_block()
                    print(block.name)
        
                     
                    cat = block.find("_em_3d_reconstruction.",["resolution"])
                    print(cat) 
                    for row in cat :
                        print(row[0])
                        resolutions[cle] = float(row[0])
                except:
                    print("Error in RNA %s" % PATH_MMCIF+cle.upper()+".cif")
                    #exit()
                
                
                    
                    
                    #print(row[9])
                
    with open("resolutions.pickle", 'wb') as fichier_pickle :
        mon_pickler = pickle.Pickler(fichier_pickle)
        mon_pickler.dump(resolutions)
    
'''16/10/19 code de Vladimir pour obtenir la resolution d'une structure PDB '''
def get_resolution(pdb_id):
    url = f'https://www.rcsb.org/structure/{pdb_id}'

    soup = BeautifulSoup(request.urlopen(url), 'html.parser')

    el = soup.find('li', id='exp_header_0_diffraction_resolution')
    print(el)
    if el != None :
        res = float(el.text.split()[-2])
    else :
        res = None
    return res 



''' calcule la proportion d'occurrences ayant une resolution inferieure a 3A sur le nombre total d'occurrences d'un type d'ARN
et stocke la liste de ces occurrences dans un fichier pickle '''
def distrib_resolution(num_ARN, type_motif):
    with open("resolutions.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        resolutions = mon_depickler.load()
        
        proportion_inf_3 = 0
        total = 0
        
        
        if isinstance(num_ARN, list) :
            for typ in num_ARN :
                with open("Nouvelles_donnees/%s/liste_representant_%s.pickle"%(type_motif,typ), 'rb') as fichier_repr :
                    mon_depickler = pickle.Unpickler(fichier_repr)
                    liste_representant = mon_depickler.load() 
                    
                    #print(resolutions)
                    for elt in liste_representant :
                        if resolutions[elt[0]] <= 3.0 :
                            proportion_inf_3 += 1
                        total += 1
        else :
            with open("Nouvelles_donnees/%s/liste_representant_%s.pickle"%(type_motif,num_ARN), 'rb') as fichier_repr :
                    mon_depickler = pickle.Unpickler(fichier_repr)
                    liste_representant = mon_depickler.load() 
                    
                    liste_representant_3a = []
                    #print(resolutions)
                    for elt in liste_representant :
                        #print(resolutions[elt[0]])
                        if resolutions[elt[0]] != None :
                            if resolutions[elt[0]] <= 3.0 :
                                proportion_inf_3 += 1
                                liste_representant_3a.append(elt)
                        else :
                            print("gros rat")
                        total += 1
        print(num_ARN)
        print(proportion_inf_3)
        print(total)         
        
        with open("Nouvelles_donnees/%s/liste_representant_%s_res_3a.pickle"%(type_motif,num_ARN), 'wb') as fichier_repr :
            mon_pickler = pickle.Pickler(fichier_repr)
            mon_pickler.dump(liste_representant_3a) 

''' 03/02/20 '''
def chaines_incompletes(num_ARN):
    liste_cles = {}
    if isinstance(num_ARN, list) :
        
        liste_representant = []
        for elt in num_ARN :
            liste_cles.update({elt : []})
            with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_repr = mon_depickler_1.load()
                
                liste_representant.extend(liste_repr)
                liste_cles[elt].extend(liste_repr)
        
    else :
        liste_cles.update({num_ARN : []})
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant = mon_depickler_1.load()
                liste_cles[num_ARN].extend(liste_representant)
    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
    
        compter_incomplets = 0
        liste_incomplets = []            
        for elt in liste_representant :
            #print(elt)
            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(elt[0], elt[1]), 'rb') as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                graphe = mon_depickler.load()
                
                i = 1
                ch_incomp = False
                deja_vu = []
                while i < 5 and not ch_incomp :
                    nb_elts = 1
                    compteur = i
                    deja_vu.append(compteur)
                    liaison_B53 = True
                    while liaison_B53 :
                        liaison_B53 = False
                        temp = compteur
                        for voisin in graphe.successors(compteur) :
                            for arc in graphe[compteur][voisin] :
                                if voisin not in [1,2,3,4] and voisin not in deja_vu and graphe[compteur][voisin][arc]["label"] == 'B53' :
                                    liaison_B53 = True
                                    temp = voisin
                                    nb_elts += graphe.nodes[voisin]["poids"]
                                    deja_vu.append(voisin)
                                    
                        for voisin in graphe.predecessors(compteur) :
                            for arc in graphe[voisin][compteur] :
                                if voisin not in [1,2,3,4] and voisin not in deja_vu and graphe[voisin][compteur][arc]["label"] == 'B53' :
                                    liaison_B53 = True
                                    temp = voisin
                                    nb_elts += graphe.nodes[voisin]["poids"]
                                    deja_vu.append(voisin)
                        compteur = temp
                    i += 1
                    if nb_elts < 4 :
                        print(elt)
                        print(i)
                        print(nb_elts)
                        ch_incomp = True
                    if elt == ('5afi', 15) :
                        print(nb_elts)
                
                if ch_incomp :
                    compter_incomplets += 1 
                    print(list(rmsd.keys())[0])
                    nom1 = "fichier_%s_%s_taille_4.pdb"%(elt[0], elt[1])
                    nom2 = "fichier_5dm6_3_taille_4.pdb"
                    if (nom1, nom2) in rmsd.keys() :
                        print(rmsd[(nom1, nom2)])
                    else :
                        print(rmsd[(nom2, nom1)])
                    liste_incomplets.append(elt)
        print(num_ARN)    
        print(compter_incomplets)
        print(liste_incomplets)       
        print(len(liste_representant))
        return liste_incomplets

def renvoi_num_ARN(elt):
    types = ['23S', '18S', '16S', 'Ribozyme', 'Riboswitch', 'SRP', '28S', '25S', 'Intron', 'arnt_16S_arnm', 'arnt_16S']
    type_elt = "" 
    for typ in types : 
        with open(NEW_EXTENSION_PATH_TAILLE+"groupe_%s.pickle"%(typ), 'rb') as fichier_type :
            mon_depickler = pickle.Unpickler(fichier_type)
            groupe = mon_depickler.load()
            compteur = 0
            if elt in groupe :
                compteur += 1
                type_elt = typ
            
            if compteur > 1 : 
                print("bizarre")
    return type_elt

def donne_sim_score_rmsd_plusieurs(paires):
    
    res = []
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
            
        with open("/media/coline/Maxtor/dico_new_110320_sim_par_branche_0.65_avec_liaison_near.pickle", 'rb') as fichier_dico_sim :
            mon_depickler_2 = pickle.Unpickler(fichier_dico_sim)
            dico_sim = mon_depickler_2.load()
                
#         with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_dico_sim_2 :
#             mon_depickler_3 = pickle.Unpickler(fichier_dico_sim_2)
#             dico_sim_2 = mon_depickler_3.load()
            
            for paire in paires :
                if paire[0][0] != '6hrm' and paire[1][0] != '6hrm' :
                    nom1 = "fichier_%s_%d_taille_4.pdb"%(paire[0][0], paire[0][1])
                    nom2 = "fichier_%s_%d_taille_4.pdb"%(paire[1][0], paire[1][1])
                    
                    liste_num_ARN = []
                    for elt in paire :
                        liste_num_ARN.append(renvoi_num_ARN(elt))
                        
                    scores = []
                    if liste_num_ARN[0] == liste_num_ARN[1] or (liste_num_ARN[0] in ['16S', '18S'] and liste_num_ARN[1] in ['16S', '18S']) or (liste_num_ARN[0] in ['23S', '25S', '28S'] and liste_num_ARN[1] in ['23S', '25S', '28S']) :
                        if liste_num_ARN[0] in ['16S', '18S'] :
                            num_ARN = ['16S', '18S']
                        elif liste_num_ARN[0] in ['23S', '25S', '28S'] :
                            num_ARN = ['23S', '25S', '28S']
                        else :
                            num_ARN = liste_num_ARN[0]
                        
                        if "needle_%s_%s_%s_%s_seq%d.txt"%(paire[0][0],paire[0][1],paire[1][0],paire[1][1], 1) in os.listdir("Nouvelles_donnees/alignements_%s/"%num_ARN) :
                            paire_vraie = paire
                        else :
                            paire_vraie = (paire[1], paire[0])
                        for k in range(1,4) : 
                            scores.append(recup_score("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq%d.txt"%(num_ARN,paire_vraie[0][0],paire_vraie[0][1],paire_vraie[1][0],paire_vraie[1][1], k)))
                    elif (liste_num_ARN[0] in ['16S', '18S'] and liste_num_ARN[1] in ['23S', '25S', '28S']) or (liste_num_ARN[1] in ['16S', '18S'] and liste_num_ARN[0] in ['23S', '25S', '28S']) :
                        if "needle_%s_%s_%s_%s_seq%d.txt"%(paire[0][0],paire[0][1],paire[1][0],paire[1][1], 1) in os.listdir("Nouvelles_donnees/alignements_%s_%s/"%(['23S', '25S', '28S'],['16S', '18S'])) :
                            paire_vraie = paire
                        else :
                            paire_vraie = (paire[1], paire[0])
                        
                        for k in range(1,4) : 
                            scores.append(recup_score("Nouvelles_donnees/alignements_%s_%s/needle_%s_%s_%s_%s_seq%d.txt"%(['23S', '25S', '28S'],['16S', '18S'] ,paire_vraie[0][0],paire_vraie[0][1],paire_vraie[1][0],paire_vraie[1][1], k)))
                    else :
                        if "needle_%s_%s_%s_%s_seq%d.txt"%(paire[0][0],paire[0][1],paire[1][0],paire[1][1], 1) in os.listdir("Nouvelles_donnees/alignements_%s/"%([['23S', '25S', '28S'], ['16S', '18S'], 'Ribozyme', 'Riboswitch', 'Intron', 'arnt_16S_arnm', 'arnt_16S'])) :
                            paire_vraie = paire
                        else :
                            paire_vraie = (paire[1], paire[0])
                        for k in range(1,4) : 
                            scores.append(recup_score("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq%d.txt"%([['23S', '25S', '28S'], ['16S', '18S'], 'Ribozyme', 'Riboswitch', 'Intron', 'arnt_16S_arnm', 'arnt_16S'] ,paire_vraie[0][0],paire_vraie[0][1],paire_vraie[1][0],paire_vraie[1][1], k)))
                    print(liste_num_ARN)
                    print("Score d'alignement de sequence :")
                    score_min = min(scores)
                    print(score_min)
                    
                    print("RMSD :")
                    if (nom1, nom2) in rmsd.keys() :
                        print(rmsd[(nom1,nom2)])
                        val_rmsd = rmsd[(nom1,nom2)]
                    else :
                        print(rmsd[(nom2,nom1)])
                        val_rmsd = rmsd[(nom2,nom1)]
                    
                    print("Similarite :")
                    if paire in dico_sim.keys() :
                        print(dico_sim[paire]["sim"])
                        #print(dico_sim_2[paire]["sim"])
                        val_sim = dico_sim[paire]["sim"]
                        #val_sim_new = dico_sim_2[paire]["sim"]
                    else :
                        print(dico_sim[(paire[1], paire[0])]["sim"])
                        #print(dico_sim_2[(paire[1], paire[0])]["sim"])
                        val_sim = dico_sim[(paire[1], paire[0])]["sim"]
                        #val_sim_new = dico_sim_2[(paire[1], paire[0])]["sim"]
                    
                    res.append((paire, liste_num_ARN[0], liste_num_ARN[1], score_min, val_rmsd, val_sim))
    return res

def donne_sim_score_rmsd(paire):
    print(paire)
    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
        nom1 = "fichier_%s_%d_taille_4.pdb"%(paire[0][0], paire[0][1])
        nom2 = "fichier_%s_%d_taille_4.pdb"%(paire[1][0], paire[1][1])
        
#         with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_dico_sim :
#             mon_depickler_2 = pickle.Unpickler(fichier_dico_sim)
#             dico_sim = mon_depickler_2.load()
            
        with open("/media/coline/Maxtor/dico_new_100320_sim_par_branche_0.65.pickle", 'rb') as fichier_dico_sim_2 :
            mon_depickler_3 = pickle.Unpickler(fichier_dico_sim_2)
            dico_sim = mon_depickler_3.load()
            
            liste_num_ARN = []
            for elt in paire :
                liste_num_ARN.append(renvoi_num_ARN(elt))
            
            scores = []
            if liste_num_ARN[0] == liste_num_ARN[1] or (liste_num_ARN[0] in ['16S', '18S'] and liste_num_ARN[1] in ['16S', '18S']) or (liste_num_ARN[0] in ['23S', '25S', '28S'] and liste_num_ARN[1] in ['23S', '25S', '28S']) :
                if liste_num_ARN[0] in ['16S', '18S'] :
                    num_ARN = ['16S', '18S']
                elif liste_num_ARN[0] in ['23S', '25S', '28S'] :
                    num_ARN = ['23S', '25S', '28S']
                else :
                    num_ARN = liste_num_ARN[0]
                
                if "needle_%s_%s_%s_%s_seq%d.txt"%(paire[0][0],paire[0][1],paire[1][0],paire[1][1], 1) in os.listdir("Nouvelles_donnees/alignements_%s/"%num_ARN) :
                    paire_vraie = paire
                else :
                    paire_vraie = (paire[1], paire[0])
                for k in range(1,4) : 
                    scores.append(recup_score("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq%d.txt"%(num_ARN,paire_vraie[0][0],paire_vraie[0][1],paire_vraie[1][0],paire_vraie[1][1], k)))
            elif (liste_num_ARN[0] in ['16S', '18S'] and liste_num_ARN[1] in ['23S', '25S', '28S']) or (liste_num_ARN[1] in ['16S', '18S'] and liste_num_ARN[0] in ['23S', '25S', '28S']) :
                if "needle_%s_%s_%s_%s_seq%d.txt"%(paire[0][0],paire[0][1],paire[1][0],paire[1][1], 1) in os.listdir("Nouvelles_donnees/alignements_%s_%s/"%(['23S', '25S', '28S'],['16S', '18S'])) :
                    paire_vraie = paire
                else :
                    paire_vraie = (paire[1], paire[0])
                
                for k in range(1,4) : 
                    scores.append(recup_score("Nouvelles_donnees/alignements_%s_%s/needle_%s_%s_%s_%s_seq%d.txt"%(['23S', '25S', '28S'],['16S', '18S'] ,paire_vraie[0][0],paire_vraie[0][1],paire_vraie[1][0],paire_vraie[1][1], k)))
            else :
                if "needle_%s_%s_%s_%s_seq%d.txt"%(paire[0][0],paire[0][1],paire[1][0],paire[1][1], 1) in os.listdir("Nouvelles_donnees/alignements_%s/"%([['23S', '25S', '28S'], ['16S', '18S'], 'Ribozyme', 'Riboswitch', 'Intron', 'arnt_16S_arnm', 'arnt_16S'])) :
                    paire_vraie = paire
                else :
                    paire_vraie = (paire[1], paire[0])
                for k in range(1,4) : 
                    scores.append(recup_score("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq%d.txt"%([['23S', '25S', '28S'], ['16S', '18S'], 'Ribozyme', 'Riboswitch', 'Intron', 'arnt_16S_arnm', 'arnt_16S'] ,paire_vraie[0][0],paire_vraie[0][1],paire_vraie[1][0],paire_vraie[1][1], k)))
            print(liste_num_ARN)
            print("Score d'alignement de sequence :")
            score_min = min(scores)
            print(score_min)
            
            print("RMSD :")
            if (nom1, nom2) in rmsd.keys() :
                print(rmsd[(nom1,nom2)])
                val_rmsd = rmsd[(nom1,nom2)]
            else :
                print(rmsd[(nom2,nom1)])
                val_rmsd = rmsd[(nom2,nom1)]
            
            print("Similarite :")
            if paire in dico_sim.keys() :
                print(dico_sim[paire]["sim"])
                #print(dico_sim_2[paire]["sim"])
                val_sim = dico_sim[paire]["sim"]
                #val_sim_new = dico_sim_2[paire]["sim"]
            else :
                print(dico_sim[(paire[1], paire[0])]["sim"])
                #print(dico_sim_2[(paire[1], paire[0])]["sim"])
                val_sim = dico_sim[(paire[1], paire[0])]["sim"]
                #val_sim_new = dico_sim_2[(paire[1], paire[0])]["sim"]
             
            return score_min, val_rmsd, val_sim

''' 14/02/20 '''
def distrib_sim_rmsd_nombre(num_ARN, liste_homologues, seuil_rmsd, seuil_sim, nb_lignes, nb_colonnes):
    for groupe in liste_homologues :
        if ('2r8s', 1) in groupe :
            groupe.append(('6d8o', 2))
    
    with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
        mon_depickler = pickle.Unpickler(fichier_rmsd)
        rmsd = mon_depickler.load()
        
        mini = min([x for x in rmsd.values() if x != None])
        maxi = max([x for x in rmsd.values() if x != None])
        
    #with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_sim :    
    #with open("/media/coline/Maxtor/dico_new_110320_sim_par_branche_0.65_avec_liaison_near.pickle", 'rb') as fichier_sim :
    #with open("/media/coline/Maxtor/dico_new_200320_sim_par_branche_0.65_avec_CWW.pickle", 'rb') as fichier_sim :
    with open("Resultats/dico_heuri_sim_par_branche_0.65_sans_liaison_near.pickle", 'rb') as fichier_sim :
    
        mon_depickler = pickle.Unpickler(fichier_sim)
        dico_sim = mon_depickler.load()
        
    #liste_cles = {}
    if isinstance(num_ARN, list) :
        
        liste_representant = []
        dico_pos = {}
        dico_score = {}
        
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
            mon_depickler = pickle.Unpickler(fichier_score)
            dico_score.update(mon_depickler.load())
        
        if num_ARN in [['23S', '25S', '28S'], ['16S', '18S']] :
            print("rapala")
            with open("dico_min_pos_score_%s_seuil_30.pickle"%(num_ARN), "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos.update(mon_depickler.load())
                
            
        for elt in num_ARN :
            
            #liste_cles.update({elt : []})
            if isinstance(elt, list) :
                
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_score.update(mon_depickler.load()) 
                
                if elt in [['23S', '25S', '28S'], ['16S', '18S']] :
                    print("rapala")
                    with open("dico_min_pos_score_%s_seuil_30.pickle"%(elt), "rb") as fichier_score_pos :
                        mon_depickler = pickle.Unpickler(fichier_score_pos)
                        dico_pos.update(mon_depickler.load())
                        
                    
                for e in elt :
                    with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%e, 'rb') as fichier_num_arn :
                        mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                        liste_repr = mon_depickler_1.load()
                        
                        liste_representant.extend(liste_repr)
            else :
                with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%elt, 'rb') as fichier_num_arn :
                    mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                    liste_representant.extend(mon_depickler_1.load())
                    
                with open("dico_min_tot_score_%s.pickle"%elt, "rb") as fichier_score :
                    mon_depickler = pickle.Unpickler(fichier_score)
                    dico_score.update(mon_depickler.load()) 
                #liste_cles[elt].extend(liste_repr)
        
    else :
        #liste_cles.update({num_ARN : []})
        with open("Nouvelles_donnees/liste_representant_%s_res_3a.pickle"%num_ARN, 'rb') as fichier_num_arn :
                mon_depickler_1 = pickle.Unpickler(fichier_num_arn)
                liste_representant = mon_depickler_1.load()
                #liste_cles[num_ARN].extend(liste_representant)
                
        with open("dico_min_pos_score_%s_seuil_30.pickle"%num_ARN, "rb") as fichier_score_pos :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_pos = mon_depickler.load()
                
        with open("dico_min_tot_score_%s.pickle"%num_ARN, "rb") as fichier_score :
                mon_depickler = pickle.Unpickler(fichier_score_pos)
                dico_score = mon_depickler.load()

    with open("fichier_forte_sim_faible_rmsd_non_homol_rapoulou.pickle", "w", newline="") as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(["Paire","Type ARN 1", "Type ARN 2", "Sim", "RMSD", "Score d'alignement", "NB liaison near 1", "NB liaison near 2"])
        
        with open("fichier_diff_type.pickle", "w", newline="") as fichier_csv_2 :
            csvwriter2 = csv.writer(fichier_csv_2)
            csvwriter2.writerow(["Type ARN 1", "Type ARN 2"])
        
            # 23S
            liste_verte = [(('1vqo', 12), ('6nd6', 7)),(('4y4o', 40), ('4u4r', 8)),(('6eri', 9), ('4u4r', 8)), (('6hma', 15), ('1vq8', 15)), (('5dm6', 2), ('6eri', 1)), (('4u27', 5), ('1vqp', 19)),(('4u27', 5), ('3ccr', 4)),(('4u27', 5), ('3cc7', 4)),(('5dm6', 7), ('1vq8', 11)),(('5dm6', 7), ('3ccr', 12)),(('5dm6', 7), ('3ccl', 11)),(('5dm6', 7), ('3ccs', 12)),(('5dm6', 7), ('3cc7', 12)),(('5dm6', 7), ('3ccq', 12)),(('5dm6', 7), ('3cce', 9)),(('5dm6', 7), ('3ccu', 11)),(('4ybb', 8), ('1vq8', 11)),(('4ybb', 8), ('3ccr', 12)),(('4ybb', 8), ('3ccl', 11)),(('4ybb', 8), ('3ccs', 12)),(('4ybb', 8), ('3cc7', 12)),(('4ybb', 8), ('3ccq', 12)),(('4ybb', 8), ('3cce', 9)),(('4ybb', 8), ('3ccu', 11)),(('3cc2', 2), ('6hma', 1)),(('5dm6', 13), ('1vqp', 19)),(('5dm6', 13), ('3ccr', 4)),(('5dm6', 13), ('3cc7', 4)), (('6hma', 9), ('4u4r', 8)), (('4y4o', 25), ('1vq8', 5)),(('4y4o', 25), ('3cc2', 2)),(('6hma', 15), ('4u4r', 16)), (('1vqp', 19), ('6eri', 9)), (('3ccr', 4), ('6eri', 9)),(('6eri', 9), ('3cc7', 4)),(('5afi', 15), ('4u4r', 28)),(('6hma', 9), ('1vqp', 19)),(('6hma', 9), ('3ccr', 4)),(('6hma', 9), ('3cc7', 4)),(('1vq8', 11), ('6hma', 6)),(('3ccr', 12), ('6hma', 6)),(('6hma', 6), ('3ccl', 11)),(('6hma', 6), ('3ccs', 12)),(('6hma', 6), ('3cc7', 12)),(('6hma', 6), ('3ccq', 12)),(('6hma', 6), ('3cce', 9)),(('6hma', 6), ('3ccu', 11)),(('4y4o', 40), ('1vqp', 19)),(('4y4o', 40), ('3ccr', 4)),(('4y4o', 40), ('3cc7', 4)),(('4y4o', 3), ('3ccl', 11)),(('4y4o', 3), ('3cce', 9)),(('4v51', 22), ('3ccl', 11)),(('4v51', 22), ('3cce', 9)),(('5e81', 12), ('3ccl', 11)),(('5e81', 12), ('3cce', 9)),(('3ccl', 11), ('4v90', 13)),(('4v90', 13), ('3cce', 9)),(('1vqo', 12), ('4u4r', 15)),(('4ybb', 1), ('1vq8', 5)),(('4y4o', 3), ('3ccr', 12)),(('4y4o', 3), ('3ccs', 12)),(('4y4o', 3), ('3cc7', 12)),(('4y4o', 3), ('3ccu', 11)),(('4v51', 22), ('3ccr', 12)),(('4v51', 22), ('3ccs', 12)),(('4v51', 22), ('3cc7', 12)),(('4v51', 22), ('3ccu', 11)),(('5e81', 12), ('3ccr', 12)),(('5e81', 12), ('3ccs', 12)),(('5e81', 12), ('3cc7', 12)),(('5e81', 12), ('3ccu', 11)),(('3ccr', 12), ('4v90', 13)),(('4v90', 13), ('3ccs', 12)),(('4v90', 13), ('3cc7', 12)),(('4v90', 13), ('3ccu', 11)),(('4y4o', 3), ('3ccq', 12)),(('4v51', 22), ('3ccq', 12)),(('5e81', 12), ('3ccq', 12)),(('4v90', 13), ('3ccq', 12)),(('4ybb', 1), ('3cc2', 2)),(('4y4o', 3), ('1vq8', 11)),(('4v51', 22), ('1vq8', 11)),(('5e81', 12), ('1vq8', 11)),(('1vq8', 5), ('6hma', 1)),(('1vq8', 11), ('4v90', 13)), (('5dm6', 9), ('1vq8', 21)), (('6hma', 2), ('4u4r', 19)), (('4ybb', 6), ('1yhq', 24)), (('6hma', 15), ('3cc2', 6))]
            # 16S
            liste_verte.extend([(('5nwy', 17), ('4y4o', 4)), (('6ek0', 2), ('6az1', 1)), (('4u27', 54), ('4y4o', 22)), (('4ybb', 35), ('4u3u', 22))])
            ## homologues cherches a la main avec bonne sim et rmsd
            liste_verte.extend([(('4y4o', 13), ('4u3u', 22)),(('3cc2', 1), ('4u4r', 36)),(('1vqo', 14), ('4u4r', 36)),(('6eri', 3), ('4u3u', 2)),(('4ybb', 6), ('4u3u', 2)),(('4y4o', 9), ('4u4r', 11)),(('4v67', 23), ('4u4r', 11)),(('4v67', 23), ('6ek0', 6)),(('4y4o', 9), ('6ek0', 6)),(('4y4o', 8), ('3cc2', 1)),(('3t1y', 8), ('4u4r', 11)),(('3t1y', 8), ('6ek0', 6)),(('2qex', 19), ('4u3u', 2)),(('1yhq', 24), ('4u3u', 2)),(('4y4o', 13), ('6az1', 6)),(('4y4o', 8), ('1vqo', 14)),(('5ngm', 1), ('4u4r', 11)),(('6eri', 15), ('4u4r', 11)),(('5ngm', 1), ('6ek0', 6)),(('6eri', 15), ('6ek0', 6)),(('4ybb', 9), ('4u4r', 11)),(('5j7l', 10), ('4u4r', 11)),(('4u27', 32), ('4u4r', 11)),(('5j7l', 10), ('6ek0', 6)),(('4ybb', 9), ('6ek0', 6)),(('4u27', 32), ('6ek0', 6)),(('1vq8', 19), ('4u3u', 2)),(('4ybb', 35), ('6az1', 6)),(('2vqe', 5), ('4u4r', 11)),(('2vqe', 5), ('6ek0', 6))])
            
            ## homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim  a 0.75 et pas avec la rmsd a 2
            liste_verte.extend([(('4y4o', 25),('1vq8', 5)),(('4y4o', 25),('3cc2', 2)),(('4ybb', 1),('1vq8', 5)),(('4ybb', 1),('3cc2', 2)), (('6hma', 1),('1vq8', 5)),(('6hma', 1),('3cc2', 2))])
            
            ## homologues du 05/03/20 : cherches dans les bleus totaux faux neg, faux pos et vrais pos
            liste_verte.extend([(('4u27', 54), ('6eri', 14)), (('2r8s', 1), ('1l8v', 1)), (('4y4o', 51), ('4u4r', 21)), (('4ybb', 7), ('4u4r', 21)),(('5ngm', 4), ('4u4r', 21)),(('4y4o', 11), ('4u4r', 24)), (('4y4o', 51), ('6az1', 3)), (('4ybb', 7), ('6az1', 3)), (('5ngm', 4), ('6az1', 3)), (('1vqp', 19), ('4u4r', 8)), (('4u4u', 11), ('6az1', 1)), (('4y4o', 8), ('4u4r', 36)), (('3ccr', 4), ('4u4r', 8)),(('3cc7', 4), ('4u4r', 8)), (('4y4o', 11), ('6eri', 10)), (('1gid', 1), ('1l8v', 1))])
            
            # 23S
            liste_rouge = [(('1vq8', 18), ('5dm6', 5)),(('5dm6', 5), ('3ccl', 14)),(('5dm6', 5), ('3cce', 12)),(('1vq8', 1), ('6eri', 2)), (('5wfs', 19), ('4u4r', 6)),(('4ybb', 12), ('4y4o', 35)),(('4ybb', 12), ('4wsd', 40)),(('4y4o', 34), ('2qex', 12)),(('4y4o', 28), ('5wfs', 6)),(('4v67', 43), ('2qex', 12)), (('4y4o', 50), ('6ek0', 4)), (('5dm6', 4), ('4ybb', 17)), (('4w2g', 52), ('6ek0', 10)), (('4y4o', 47), ('1vq8', 19)),(('4y4o', 47), ('1yhq', 24)),(('4y4o', 47), ('2qex', 19)),(('4y4o', 47), ('4u4r', 13)),(('4ybb', 13), ('1vqo', 8)),(('4w2f', 23), ('4u4r', 36)),(('1vqo', 8), ('5dm6', 1)),(('1vqo', 8), ('1mms', 1)),(('1vq8', 19), ('6eri', 1)),(('1yhq', 24), ('6eri', 1)),(('6eri', 1), ('2qex', 19)), (('5dm6', 14), ('4y4o', 8)),(('4ybb', 17), ('4y4o', 8)), (('4ybb', 30), ('2zjr', 15)),(('2zjr', 15), ('5wfs', 11)), (('2zjr', 15), ('6eri', 8)), (('4y4o', 12), ('2zjr', 15)), (('4w2g', 52), ('4y4o', 28)), (('4y4o', 28), ('4v67', 7))]
            #liste_rouge.extend([(('4ybb', 12), ('4y4o', 38)),(('4ybb', 12), ('4ybb', 21)),(('4ybb', 12), ('6ek0', 7)),(('4v67', 7), ('4y4o', 38)),(('4v67', 7), ('4ybb', 21)),(('4v67', 7), ('6ek0', 7)),(('5wfs', 19), ('4y4o', 38)),(('5wfs', 19), ('4ybb', 21)),(('5wfs', 19), ('6ek0', 7))])
            ## non homologues meme type rmsd 2-2.5, sim 0.7-0.75
            liste_rouge.extend([(('4y4o', 23), ('6h4n', 24)),(('2zjr', 3), ('4u4r', 18)),(('4u27', 54), ('4y4o', 38)),(('4u27', 54), ('6ek0', 7)),(('4u27', 54), ('4u3u', 7)), (('4y4o', 23), ('4u27', 3)),(('4y4o', 23), ('5afi', 24))])
            
            ## les non homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim à 0.75 mais pas avec la rmsd à 2
            liste_rouge.extend([(('4ybb', 12), ('6ek0', 10)),(('4w2g', 52), ('6ek0', 10)),(('2zjr', 3),('6ek0', 10)),(('2zjr', 3),('4u4r', 18)),(('4w2g', 52),('4u4r', 18)),(('4v67', 7), ('4u4r', 18)),(('5wfs', 19), ('4u4r', 18)),(('4v67', 7),('6ek0', 10)),(('5wfs', 19),('6ek0', 10)),(('4y4o', 23),('6h4n', 24)),(('4y4o', 23), ('4u27', 3)),(('4y4o', 23),('5afi', 24))])
            #liste_rouge = []
            
            ## non-homologues du 05/03/20 : cherches dans les bleus totaux faux neg, faux pos et vrais pos (voir fichier excel)
            liste_rouge.extend([(('4u27', 54), ('4ybb', 21)), (('6eri', 14), ('6eri', 15)),(('5ngm', 1), ('6eri', 14)),(('4u27', 54), ('6eri', 15)), (('6eri', 14), ('4u4r', 21)), (('4y4o', 23), ('1vq8', 12)),(('4y4o', 23), ('3cc2', 19)),(('4y4o', 23), ('3cpw', 14)), (('4y4o', 22), ('6ek0', 7)),(('4y4o', 22), ('4u3u', 7)),(('1yhq', 9), ('5afi', 24)),(('5afi', 24), ('2qex', 12)),(('4ybb', 37), ('4ybb', 1)),(('4ybb', 1), ('4y4o', 19)),(('4ybb', 1), ('1vq8', 16)),(('4ybb', 1), ('6h4n', 7)),(('4ybb', 1), ('5wfs', 19)),(('4ybb', 1), ('6eri', 16)),(('4ybb', 1), ('6qul', 4)),(('4ybb', 1), ('3ccm', 16)),(('4ybb', 1), ('4u4r', 6)),(('4ybb', 1), ('6ek0', 3)),(('4ybb', 12), ('1vq8', 16)),(('4w2g', 52), ('1vq8', 16)),(('4w2g', 52), ('3cd6', 14)),(('4w2g', 52), ('2qex', 20)),(('2zjr', 3), ('1vq8', 16)),(('4u27', 3), ('1yhq', 9)),(('4u27', 3), ('2qex', 12)),(('1vq8', 8), ('1vq8', 12)),(('1vq8', 8), ('3cc2', 19)),(('1vq8', 8), ('3cpw', 14)),(('5dm6', 4), ('6hma', 8)),(('4ybb', 12), ('6qul', 4)),(('4w2g', 52), ('6qul', 4)),(('6hma', 8), ('6hma', 14)),(('6hma', 8), ('1vq8', 18)),(('6hma', 8), ('3ccl', 14)),(('6hma', 8), ('6eri', 12)),(('6hma', 8), ('3cce', 12)),(('6hma', 8), ('4u4r', 1)),(('2zjr', 3), ('6qul', 4)),(('1vq8', 8), ('5afi', 24)),(('1vq8', 12), ('1yhq', 9)),(('1vq8', 12), ('2qex', 12)),(('1yhq', 9), ('3cc2', 19)),(('1yhq', 9), ('6h4n', 24)),(('1yhq', 9), ('3cpw', 14)),(('3cc2', 19), ('2qex', 12)),(('6h4n', 24), ('2qex', 12)),(('3cpw', 14), ('2qex', 12)),(('5wfs', 19), ('6qul', 4)), (('4y4o', 22), ('4y4o', 38)),(('4y4o', 22), ('4ybb', 21)), (('4u27', 3), ('1vq8', 8)),(('1vq8', 8), ('6h4n', 24)), (('4ybb', 37), ('4y4o', 25)),(('4ybb', 37), ('6hma', 1)),(('4y4o', 25), ('1vq8', 16)),(('4y4o', 25), ('3ccr', 17)),(('4y4o', 25), ('3cd6', 14)),(('4y4o', 25), ('4u4r', 6)),(('4y4o', 19), ('6hma', 1)),(('1vq8', 16), ('6hma', 1)),(('1yhq', 17), ('6hma', 1)),(('6h4n', 7), ('6hma', 1)),(('3ccr', 17), ('6hma', 1)),(('6hma', 1), ('6hma', 10)),(('6hma', 1), ('3cd6', 14)),(('6hma', 1), ('3ccs', 16)),(('6hma', 1), ('6eri', 16)),(('6hma', 1), ('3cc7', 17)),(('6hma', 1), ('6qul', 4)),(('6hma', 1), ('3ccq', 18)),(('6hma', 1), ('3ccm', 16)),(('6hma', 1), ('3ccu', 16)),(('6hma', 1), ('2qex', 20)),(('6hma', 1), ('4u4r', 6)),(('6hma', 1), ('6ek0', 3)),(('5dm6', 4), ('5dm6', 9)),(('5dm6', 4), ('5afi', 17)),(('5dm6', 9), ('6eri', 12)),(('4y4o', 28), ('6eri', 12)),(('5afi', 17), ('6eri', 12)),(('4ybb', 12), ('6ek0', 3)),(('4w2g', 52), ('1yhq', 17)),(('4w2g', 52), ('6hma', 10)),(('4w2g', 52), ('3cc7', 17)),(('4w2g', 52), ('6ek0', 3)),(('4y4o', 25), ('4y4o', 19)),(('4y4o', 25), ('1yhq', 17)),(('4y4o', 25), ('6hma', 10)),(('4y4o', 25), ('3ccs', 16)),(('4y4o', 25), ('6eri', 16)),(('4y4o', 25), ('3cc7', 17)),(('4y4o', 25), ('6qul', 4)),(('4y4o', 25), ('3ccq', 18)),(('4y4o', 25), ('3ccm', 16)),(('4y4o', 25), ('3ccu', 16)),(('4y4o', 25), ('6ek0', 3)),(('2zjr', 3), ('1yhq', 17)),(('2zjr', 3), ('6hma', 10)),(('2zjr', 3), ('3cc7', 17)),(('2zjr', 3), ('6ek0', 3)),(('5wfs', 19), ('6ek0', 3)),(('5dm6', 9), ('6hma', 14)),(('5dm6', 9), ('1vq8', 18)),(('5dm6', 9), ('3ccl', 14)),(('5dm6', 9), ('3cce', 12)),(('5dm6', 9), ('4u4r', 1)),(('4y4o', 28), ('6hma', 14)),(('4y4o', 28), ('1vq8', 18)),(('4y4o', 28), ('4u4r', 1)), (('4ybb', 12), ('1vq8', 3)),(('4ybb', 12), ('3cc2', 21)),(('4ybb', 12), ('3cpw', 3)), (('5dm6', 4), ('4y4o', 28)),(('4ybb', 12), ('6h4n', 7)),(('4w2g', 52), ('6h4n', 7)),(('4w2f', 37), ('4y4o', 28)),(('4y4o', 28), ('4y4o', 43)),(('4y4o', 28), ('3ccl', 14)),(('4y4o', 28), ('3cce', 12)),(('4y4o', 25), ('6h4n', 7)), (('2zjr', 3), ('6h4n', 7)), (('4u27', 54), ('5ngm', 1)), (('6eri', 14), ('6ek0', 7)),(('6eri', 14), ('4u3u', 7))])
            
            liste_gris = [(('4l81', 2), ('3iqr', 1)),(('4l81', 2), ('5fk2', 1)), (('4l81', 2), ('5fk4', 1)),(('4l81', 2), ('3iqn', 1)),(('5fjc', 1), ('4l81', 2)),(('4l81', 2), ('5fke', 1)),(('5fk1', 1), ('4l81', 2)),(('4l81', 2), ('2ygh', 1)),(('2ydh', 1), ('4l81', 2)),(('4l81', 2), ('5fkd', 1)),(('4l81', 2), ('5fkf', 1))]
            liste_vrais_positifs = []
            
            distrib_sim = []
            distrib_rmsd = []
            distrib_sim_rouge = []
            distrib_sim_verte = []
            distrib_rmsd_rouge = []
            distrib_rmsd_verte = []
            distrib_sim_pos = []
            distrib_rmsd_pos = []
            
            compteur = 0
            compter = 0
            compteer = 0
            dico_type = {}
            for elt in liste_representant :
                dico_type.update({elt : renvoi_num_ARN(elt)})
                
            matrice_nombre_pos = [[ 0 for _ in range(nb_colonnes+1) ] for _ in range(nb_lignes+1) ] 
            matrice_nombre_vert = [[ 0 for _ in range(nb_colonnes+1) ] for _ in range(nb_lignes+1) ] 
            matrice_nombre_rouge = [[ 0 for _ in range(nb_colonnes+1) ] for _ in range(nb_lignes+1) ]
            matrice_nombre_bleu = [[ 0 for _ in range(nb_colonnes+1) ] for _ in range(nb_lignes+1) ] 
            matrice_nombre_homologues = [[ 0 for _ in range(nb_colonnes+1) ] for _ in range(nb_lignes+1) ] 
            matrice_gris = [[ 0 for _ in range(nb_colonnes+1) ] for _ in range(nb_lignes+1) ] 
            compte_rouge = 0
            
            print("combien d'elts ?")
            print(len(liste_representant))
            
            for i in range(len(liste_representant)) :
                for j in range(i+1, len(liste_representant)) :
                    if liste_representant[i][0] != '6hrm' and  liste_representant[j][0] != '6hrm' and liste_representant[i] not in liste_pbs and liste_representant[j] not in liste_pbs :
                        
                        
                        ok = True
                        for groupe in liste_homologues :
                            if liste_representant[i] in groupe and liste_representant[j] in groupe :
                                ok = False
                                
                        if ok :
                            
    #                         type1 = renvoi_num_ARN(liste_representant[i])
    #                         type2 = renvoi_num_ARN(liste_representant[j])
                            type1 = dico_type[liste_representant[i]]
                            type2 = dico_type[liste_representant[j]]
                            diff_type = False
                            if (type1 != type2 and type1 not in ["23S","28S", "25S", "16S", "18S"] and type2 not in ["23S","28S", "25S", "16S", "18S"]) or ((type1 in ["23S","28S", "25S"] and type2 not in ["23S","28S","25S"]) or (type2 in ["23S","28S", "25S"] and type1 not in ["23S","28S","25S"]) ) or ((type1 in ["16S","18S"] and type2 not in  ["16S","18S"]) or (type2 in ["16S", "18S"] and type1 not in ["16S", "18S"])) :
                                diff_type = True
                                print(type1)
                                print(type2)
                                csvwriter2.writerow([type1, type2])
                            nom1 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[i][0], liste_representant[i][1])
                            nom2 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[j][0], liste_representant[j][1])
                            compteur += 1
                            if (nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)] <= seuil_rmsd:
                                
                                if (liste_representant[i], liste_representant[j]) in dico_pos.keys() or (liste_representant[j], liste_representant[i]) in dico_pos.keys() :
                                   
                                    
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                        matrice_nombre_pos[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
                                        distrib_sim_pos.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                        distrib_rmsd_pos.append(rmsd[(nom1,nom2)])
                                        
                                    elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        matrice_nombre_pos[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1                                    
                                        distrib_sim_pos.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                        distrib_rmsd_pos.append(rmsd[(nom1,nom2)])
                                elif (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
                                    
                                    
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
                                        matrice_nombre_vert[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
                                        distrib_sim_verte.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                        distrib_rmsd_verte.append(rmsd[(nom1,nom2)])
                                    elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        matrice_nombre_vert[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
    
                                        distrib_sim_verte.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                        distrib_rmsd_verte.append(rmsd[(nom1,nom2)])
                                elif (liste_representant[i], liste_representant[j]) in liste_gris or (liste_representant[j], liste_representant[i]) in liste_gris :
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
                                        matrice_gris[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1

                                    elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        matrice_gris[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1

  
                                
                                #elif (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge  or  diff_type:
                                else :
                                    if liste_representant[i] in [('6eri', 17), ('5e81', 27)] and liste_representant[j] in [('6eri', 17), ('5e81', 27)] :
                                        print(liste_representant[i], liste_representant[j])
                                        print(int(((rmsd[(nom1,nom2)]-mini)/(maxi-mini))*nb_lignes))
                                        print(((rmsd[(nom1,nom2)]-mini)/(maxi-mini))*nb_lignes)
                                        #exit()
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                        if dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] > 0.7 : 
                                            print("gros tas")
                                            print((liste_representant[i], liste_representant[j]))
                                        print(int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes))
                                        print(int(rmsd[(nom1,nom2)]/16*nb_lignes))
                                        print(len(matrice_nombre_rouge))
                                        print(len(matrice_nombre_rouge[0]))
                                        matrice_nombre_rouge[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
                                        compte_rouge += 1
                                        distrib_sim_rouge.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                        distrib_rmsd_rouge.append(rmsd[(nom1,nom2)])
                                        
                                        if dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.75 and rmsd[(nom1, nom2)] <= 2.5  :
                                        #if (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] <= 0.65 and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= 0.6 and rmsd[(nom1, nom2)] <= 6.5 and rmsd[(nom1, nom2)] >= 6.0) : # or (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] <= 0.6 and rmsd[(nom1, nom2)] <= 2.5):
                                        #if rmsd[(nom1, nom2)] <= 2 :
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
                                                mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                                graphe1 = mon_depickler_1.load()
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
                                                mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                                graphe2 = mon_depickler_2.load()
                                            csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])
                                            liste_vrais_positifs.append((liste_representant[i], liste_representant[j]))
                                        
                                    elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        if dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] > 0.7 : 
                                            print("gros tas")
                                            print((liste_representant[i], liste_representant[j]))
                                        matrice_nombre_rouge[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
                                        compte_rouge += 1
                                        distrib_sim_rouge.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                        distrib_rmsd_rouge.append(rmsd[(nom1,nom2)])
                                        
                                        if dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.75 and rmsd[(nom1, nom2)] <= 2.5  :
                                        #if (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] <= 0.65 and dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.6 and rmsd[(nom1, nom2)] <= 6.5 and rmsd[(nom1, nom2)] >= 6.0) : #or (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] <= 0.6 and rmsd[(nom1, nom2)] <= 2.5):
                                        #if rmsd[(nom1, nom2)] <= 2 :
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
                                                mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                                graphe1 = mon_depickler_1.load()
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
                                                mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                                graphe2 = mon_depickler_2.load()
      
                                            csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])
                                            liste_vrais_positifs.append((liste_representant[i], liste_representant[j]))
    
#                                 else :
#                                      
#                                     if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
#     #                                     if rmsd[(nom1, nom2)] <= 4 and dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
#     #                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)] ])
#                                         distrib_sim.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
#                                         print(int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])*nb_colonnes)
#                                         print(int(rmsd[(nom1,nom2)]*nb_lignes))
#                                         matrice_nombre_bleu[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
#       
#                                         distrib_rmsd.append(rmsd[(nom1,nom2)])
#                                           
#                                           
#                                         if (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.75 and rmsd[(nom1, nom2)] >= 2) or (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.6 and rmsd[(nom1, nom2)] <= 2.5):
#                                         #if dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.7 :
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
#                                                 mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
#                                                 graphe1 = mon_depickler_1.load()
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
#                                                 mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#                                                 graphe2 = mon_depickler_2.load()
#                                             csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])
#   
#                                           
#                                     elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
#                                         #if rmsd[(nom1, nom2)] <= 4 and dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
#                                         #    csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)] ])
#                                         matrice_nombre_bleu[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
#       
#                                         distrib_sim.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
#                                         distrib_rmsd.append(rmsd[(nom1,nom2)])
#                                           
#                                           
#                                         if (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.75 and rmsd[(nom1, nom2)] >= 2) or (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.6 and rmsd[(nom1, nom2)] <= 2.5):
#                                         #if dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.7 :
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
#                                                 mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
#                                                 graphe1 = mon_depickler_1.load()
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
#                                                 mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#                                                 graphe2 = mon_depickler_2.load()
#        
#                                             csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])
                                        
                            elif (nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)] <= seuil_rmsd :

                                if (liste_representant[i], liste_representant[j]) in dico_pos.keys() or (liste_representant[j], liste_representant[i]) in dico_pos.keys() :
                                    
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
                                        matrice_nombre_pos[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
    
                                        distrib_sim_pos.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                        distrib_rmsd_pos.append(rmsd[(nom2,nom1)])
                                    elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim:
                                        matrice_nombre_pos[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
                                        distrib_sim_pos.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                        distrib_rmsd_pos.append(rmsd[(nom2,nom1)])
                                elif (liste_representant[i], liste_representant[j]) in liste_verte or (liste_representant[j], liste_representant[i]) in liste_verte :
#                                     if liste_representant[i] in [('4u4r', 11), ('4v67', 23)] and liste_representant[j] in [('4u4r', 11), ('4v67', 23)] :
#                                         print("raaaa")
#                                         exit()
                                    
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                        matrice_nombre_vert[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
                                        distrib_sim_verte.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                        distrib_rmsd_verte.append(rmsd[(nom2,nom1)])
                                    elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        matrice_nombre_vert[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
                                        distrib_sim_verte.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                        distrib_rmsd_verte.append(rmsd[(nom2,nom1)])
                                elif (liste_representant[i], liste_representant[j]) in liste_gris or (liste_representant[j], liste_representant[i]) in liste_gris :
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                        matrice_gris[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1

                                    elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        matrice_gris[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1

                                
                                #elif (liste_representant[i], liste_representant[j]) in liste_rouge or (liste_representant[j], liste_representant[i]) in liste_rouge or diff_type :
                                else :    
                                    if liste_representant[i] in [('6eri', 17), ('5e81', 27)] and liste_representant[j] in [('6eri', 17), ('5e81', 27)] :
                                        print(liste_representant[i], liste_representant[j])
                                        print(int(((rmsd[(nom2,nom1)]-mini)/(maxi-mini))*nb_lignes))
                                        print(((rmsd[(nom2,nom1)]-mini)/(maxi-mini))*nb_lignes)
                                        #exit()
                                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim:
                                        matrice_nombre_rouge[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
                                        distrib_sim_rouge.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
                                        distrib_rmsd_rouge.append(rmsd[(nom2,nom1)])
                                        compte_rouge += 1
                                        
                                        if dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.75 and rmsd[(nom2, nom1)] <= 2.5 :
                                        #if rmsd[(nom2, nom1)] <= 2 :
                                        #if (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] <= 0.65 and dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.6 and rmsd[(nom2, nom1)] <= 6.5 and rmsd[(nom2, nom1)] >= 6.0)  : #or (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] <= 0.6 and rmsd[(nom2, nom1)] <= 2.5):
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
                                                mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                                graphe1 = mon_depickler_1.load()
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
                                                mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                                graphe2 = mon_depickler_2.load()
      
                                            csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])
                                            liste_vrais_positifs.append((liste_representant[i], liste_representant[j]))
                                    
                                    elif  (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        matrice_nombre_rouge[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
                                        distrib_sim_rouge.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
                                        distrib_rmsd_rouge.append(rmsd[(nom2,nom1)])
                                        compte_rouge += 1
                                        
                                        if dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.75 and rmsd[(nom2, nom1)] <= 2.5 :
                                        #if rmsd[(nom2, nom1)] <= 2  :
                                        #if (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] <= 0.65 and dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.6 and rmsd[(nom2, nom1)] <= 6.5 and rmsd[(nom2, nom1)] >= 6.0) : #or (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] <= 0.6 and rmsd[(nom2, nom1)] <= 2.5):
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
                                                mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
                                                graphe1 = mon_depickler_1.load()
                                            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
                                                mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
                                                graphe2 = mon_depickler_2.load()
      
                                            csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])
                                            liste_vrais_positifs.append((liste_representant[i], liste_representant[j]))
                                        
                                        
                                    
                                
#                                 else :                            
#                                      
#                                     if (liste_representant[i], liste_representant[j]) in dico_sim.keys() and dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
#     #                                     if rmsd[(nom2, nom1)] <= 4 and dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
#     #                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)] ])
#                                         matrice_nombre_bleu[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi, 0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
#                                         distrib_sim.append(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"])
#                                         distrib_rmsd.append(rmsd[(nom2,nom1)])
#                                           
#                                         if (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.75 and rmsd[(nom2, nom1)] >= 2) or (dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.6 and rmsd[(nom2, nom1)] <= 2.5):
#                                         #if dico_sim[ (liste_representant[i], liste_representant[j])]["sim"] >= 0.7 :
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
#                                                 mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
#                                                 graphe1 = mon_depickler_1.load()
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
#                                                 mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#                                                 graphe2 = mon_depickler_2.load()
#        
#                                             csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])
#                                           
#                                     elif (liste_representant[j], liste_representant[i]) in dico_sim.keys() and dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
#     #                                     if rmsd[(nom2, nom1)] <= 4 and dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
#     #                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)] ])
#                                         print(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
#                                         print(rmsd[(nom2,nom1)])
#                                         print()
#                                         matrice_nombre_bleu[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
#                                         distrib_sim.append(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"])
#                                         distrib_rmsd.append(rmsd[(nom2,nom1)])
#                                           
#                                         if (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.75 and rmsd[(nom2, nom1)] >= 2) or (dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.6 and rmsd[(nom2, nom1)] <= 2.5):
#                                         #if dico_sim[ (liste_representant[j], liste_representant[i])]["sim"] >= 0.7 :
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0],liste_representant[i][1]), 'rb') as fichier_graphe1 :
#                                                 mon_depickler_1 = pickle.Unpickler(fichier_graphe1)
#                                                 graphe1 = mon_depickler_1.load()
#                                             with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0],liste_representant[j][1]), 'rb') as fichier_graphe2 :
#                                                 mon_depickler_2 = pickle.Unpickler(fichier_graphe2)
#                                                 graphe2 = mon_depickler_2.load()
#        
#                                             csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])], nb_liaison_near(graphe1), nb_liaison_near(graphe2) ])

            #                 else : 
            #                     print(nom1, nom2)
                        else :
                            nom1 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[i][0], liste_representant[i][1])
                            nom2 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[j][0], liste_representant[j][1])
                            if (nom1, nom2) in rmsd.keys() and rmsd[(nom1, nom2)] != None :
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
                                    matrice_nombre_homologues[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
#                                     if rmsd[(nom1, nom2)] >= 2.5 :
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                else :
                                    matrice_nombre_homologues[nb_lignes - int(((rmsd[(nom1,nom2)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
#                                     if rmsd[(nom1, nom2)] >= 2.5 :
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom1, nom2)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                                            
                            elif (nom2, nom1) in rmsd.keys() and rmsd[(nom2, nom1)] != None :
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
                                    matrice_nombre_homologues[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[i], liste_representant[j])]["sim"]*nb_colonnes)] += 1
#                                     if rmsd[(nom2, nom1)] >= 2.5 :
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[i], liste_representant[j])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])

                                else :
                                    matrice_nombre_homologues[nb_lignes - int(((rmsd[(nom2,nom1)])/(round(maxi,0)))*nb_lignes)][int(dico_sim[ (liste_representant[j], liste_representant[i])]["sim"]*nb_colonnes)] += 1
#                                     if rmsd[(nom2, nom1)] >= 2.5 :
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]),type1, type2, dico_sim[ (liste_representant[j], liste_representant[i])]["sim"], rmsd[(nom2, nom1)], dico_score[(liste_representant[i], liste_representant[j])] ])

                            compter += 1
                        print(compteur)
                    else :
                        compteer += 1
    print(len(distrib_sim_rouge))     
    print(compte_rouge)           
    print("min rmsd")
    print(mini)
    print("max rmsd")
    print(maxi)
    print(max(distrib_rmsd_rouge))
    print([k  for (k, val) in rmsd.items() if val != None and val == 12.198324250127802])
    print(dico_pos)
    print(liste_vrais_positifs)
    fig, ax = plt.subplots(figsize=(20,12))
    im = ax.imshow(matrice_nombre_rouge, interpolation='nearest', cmap=plt.cm.get_cmap('Reds'), vmin=0, vmax=10)
    ax.figure.colorbar(im, ax=ax, pad=0.1)
    # We want to show all ticks...
    ax.set(xticks=np.arange(len(matrice_nombre_rouge[0])),
           yticks=np.arange(len(matrice_nombre_rouge)),
           # ... and label them with the respective list entries
           xticklabels=[0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0], yticklabels=np.arange(int(maxi), -1, -0.5),
           ylabel='RMSD',
           xlabel='Similarité', 
           ) 
    
    plt.title("Distribution de la similarite en fonction de la RMSD", pad=20.0)
    ax.tick_params(axis="x", bottom=True, labelbottom=True)
    # Rotate and align bottom ticklabels
    plt.setp([tick.label1 for tick in ax.xaxis.get_major_ticks()], rotation=45,
             ha="right", va="center", rotation_mode="anchor")
    # Rotate the tick labels and set their alignment.
#     plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
#              rotation_mode="anchor")
    
    ax.tick_params(axis="y", left=True, labelleft=True)
    liste_pos = []
    for tick in ax.get_xticks():
        liste_pos.append(tick-0.5)
    ax.set_xticks(liste_pos)
    
    liste_pos = []
    for tick in ax.get_yticks():
        liste_pos.append(tick+0.5)
    ax.set_yticks(liste_pos)
    
    # Loop over data dimensions and create text annotations.
    fmt = '1.0f'
    thresh = 5
    for i in range(len(matrice_nombre_rouge)):
        for j in range(len(matrice_nombre_rouge[i])):
            if matrice_nombre_rouge[i][j] != -1 :
                ax.text(j, i, format(matrice_nombre_rouge[i][j], fmt),
                        ha="center", va="center",
                        color="white" if matrice_nombre_rouge[i][j] > thresh else "black")
                
    for i in range(len(matrice_gris)):
        for j in range(len(matrice_gris[i])):
            if matrice_gris[i][j] != -1 and matrice_gris[i][j] != 0 :
                ax.text(j, i, format(matrice_gris[i][j], fmt),
                        ha="right", va="bottom",
                        color="grey")
    fig.tight_layout()
    plt.show()
    return liste_vrais_positifs
    
def nts_modifies(types_arn):
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
    print(len(liste_tout))
    compteur = 0
    for elt in liste_tout :
        with open(NEW_EXTENSION_PATH_TAILLE + "fichier_%s_%s.pickle"%(elt[0],elt[1]),'rb') as fichier_graphe :
            mon_depickler = pickle.Unpickler(fichier_graphe)
            graphe = mon_depickler.load()
            
            with open("Graphs/"+elt[0]+".pickle", 'rb') as fichier_graphes :
                mon_depickler_g  = pickle.Unpickler(fichier_graphes)
                graphes = mon_depickler_g.load()
                
                trouve = False
                for noeud, data in graphe.nodes(data=True) :
                    if data["type"] not in [None, -1] :
                        for i in range(data["position"][0], data["position"][1]+1) :
                            if elt == ('6gyv', 1) :
                                print(graphes.nodes[(data["num_ch"], i)])
                            if graphes.nodes[(data["num_ch"], i)]["nt"] != graphes.nodes[(data["num_ch"], i)]["real_nt"] :
                                trouve = True
        if trouve : 
            compteur += 1
            print(elt)
    print(compteur)

def recup_groupes_vrais_positifs(liste, liste_homologues, seuil_sim, seuil_rmsd):
    graphe = nx.Graph()
    res = donne_sim_score_rmsd_plusieurs(liste)
    compteur = 0
    for elt in liste :
        #score_min, val_rmsd, val_sim = donne_sim_score_rmsd(elt)
        graphe.add_edge(elt[0], elt[1], aln=res[compteur][3], rmsd=res[compteur][4], sim=res[compteur][5])
        compteur += 1
    
    print(graphe.number_of_nodes())
    print(graphe.number_of_edges())
    #exit()
    a_ajouter = []
    for noeud in graphe.nodes() :
        for groupe in liste_homologues :
            if noeud in groupe :
                for elt in groupe :
                    if elt not in graphe.nodes() :
                        a_ajouter.append(elt)
                        
#     for elt in a_ajouter :
#         graphe.add_node(elt)
    
    composantes = recherche_composante_connexe(graphe)
    print(len(composantes))
    #exit()
    with open("Nouvelles_donnees/fichier_csv_vrais_positifs_sim_par_branche_0.65_%s_%s.csv"%(seuil_sim, seuil_rmsd),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["source", "target", "weight", "aln", "rmsd", "sim"])
                
                compteur = 0
                for u,v,data in graphe.edges(data=True) :
                    csvwriter.writerow([u,v, 3, data["aln"], round(data["rmsd"], 2), round(data["sim"], 2)])
                print(compteur)
                
#                 for comp in composantes :
#                     for i in range(len(comp)) :
#                         for j in range(i+1, len(comp)) :
#                             if (comp[i], comp[j]) not in graphe.edges() and (comp[j], comp[i]) not in graphe.edges() :
#                                 score_min, rmsd, sim = donne_sim_score_rmsd((comp[i], comp[j]))
#                                 csvwriter.writerow([comp[i], comp[j], 2, score_min, round(rmsd, 2), round(sim, 2)])
                                
#                 for noeud in graphe.nodes() :
#                     if len(graphe[noeud]) == 0 :
#                         for groupe in liste_homologues :
#                             if noeud in groupe :
#                                 for elt in groupe :
#                                     if elt != noeud :
#                                         score_min, rmsd, sim = donne_sim_score_rmsd((noeud, elt)) 
#                                         csvwriter.writerow([noeud, elt, 1, score_min, round(rmsd, 2), round(sim, 2)])
#                                     for voisin in graphe[elt] :
#                                         score_min, rmsd, sim = donne_sim_score_rmsd((noeud, voisin)) 
#                                         csvwriter.writerow([noeud, voisin, 1, score_min, round(rmsd, 2), round(sim, 2)])
                        
                #exit()
               
    with open("Nouvelles_donnees/fichier_csv_vrais_positifs_sim_par_branche_0.65_%s_%s_noeuds.csv"%(seuil_sim, seuil_rmsd),'w') as fichier_csv:
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["id", "libelle", "homologues"])
                
                compteur = 0
                for noeud,data in graphe.nodes(data=True) :
                    nom_struct = " "
                    #print(noeud)
                    with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(noeud[0], noeud[1]), 'rb') as fichier_graphe :
                        mon_depickler = pickle.Unpickler(fichier_graphe)
                        graphe = mon_depickler.load()
                        with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv_mol :
                            csvreader = csv.reader(fichier_csv_mol)
                            
                            for row in csvreader :
                                #print(row)
                                if len(row) == 7 :
                                    if row[1].lower() == noeud[0] :
                                        if graphe.nodes[1]["num_ch"] == row[2] : 
                                            nom_struct += row[3]
                                        elif graphe.nodes[2]["num_ch"] == row[2] : 
                                            nom_struct += row[3]
                                        elif graphe.nodes[5]["num_ch"] == row[2] :
                                            nom_struct += row[3]
                    compte_homol = 0
                    num_homol = -1
                    for groupe in liste_homologues :
                        if noeud in groupe :
                            num_homol = compte_homol
                        compte_homol += 1
                                            
                            #print(nom_struct)
                    compteur += 1
                    
                                            
                    
                    csvwriter.writerow([noeud, nom_struct, num_homol])
                    
                print(compteur)
        
def nb_liaison_near(graphe):
    nb_near = 0
    for u,v, data in graphe.edges(data=True) :   
        if data["near"] == True and (u,v) not in [(1,2), (2,1), (3,4), (4,3), (1,5), (5,1)] :
            print(u,v)
            nb_near += 1
            
    return int(nb_near/2)

def pourcentage_liaison_near(liste_num_ARN):
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 :
                            liste_tout.append((elt, element))
                            
    somme = 0
    somme_4 = 0
    somme_3 = 0
    somme_2 = 0
    somme_1 = 0
    somme_0 = 0
    distrib_nb_liaison_near = []
    for elt in liste_tout :
        with open("Nouvelles_donnees/fichier_%s_%s_4.pickle"%(elt[1][0], elt[1][1]), 'rb') as fichier_graphe :
            mon_depickler_g = pickle.Unpickler(fichier_graphe)
            graphe = mon_depickler_g.load()
            nb_near = nb_liaison_near(graphe)
            
            if elt[1] == ('5ngm', 3) :
                print("rapala")
                print(nb_near)
            
            
            distrib_nb_liaison_near.append(nb_near)
            somme += nb_near
            if nb_near == 4 :
                somme_4 += 1
                print(elt)
            elif nb_near == 3 :
                somme_3 += 1
                print(elt)
            elif nb_near == 2 :
                somme_2 += 1
            elif nb_near == 1 :
                somme_1 += 1
            elif nb_near == 0 :
                somme_0 += 1
    print(somme/len(liste_tout))
    print(min(distrib_nb_liaison_near))
    print(max(distrib_nb_liaison_near))
    print(distrib_nb_liaison_near)
    print(somme_4)
    print(somme_3)
    print(somme_2)
    print(somme_1)
    print(somme_0)
    axs = sns.distplot(distrib_nb_liaison_near, kde=False)
    axs.set_yticks(np.arange(0,300,20))
    #axs.set_xticks(np.arange(0,5,1))
    axs.set_xlabel("Nombre de liaisons near sans compter le motif")
    axs.set_ylabel("Nombre d'occurrences")
    
    plt.show()
    plt.close()
    
    
''' 01/04/20 '''
def distrib_types_liaisons_CWW(liste_num_ARN):
    
    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/liste_representant_%s.pickle"%elt, 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] <= 3 :
                            liste_tout.append(element)
    
    with open("grands_graphes_new_data_taille_4.pickle", 'rb') as fichier_ecriture :
        mon_depickler = pickle.Unpickler(fichier_ecriture)
        dico_graphes = mon_depickler.load()
    
    
    dico_types_CWW = {}
    for elt in liste_tout :
                
                print(dico_graphes[elt])
                for u,v,data in dico_graphes[elt][0].edges(data=True) :
                    if data["label"] == 'CWW' or data["label"] == 'CWWn' :
                            nt1 = dico_graphes[elt][0].nodes[u]["nt"]
                            nt2 = dico_graphes[elt][0].nodes[v]["nt"]
                            if (nt1, nt2) in dico_types_CWW.keys() :
                                dico_types_CWW[(nt1, nt2)] += 1
                            elif (nt2, nt1) in dico_types_CWW.keys() :
                                dico_types_CWW[(nt2, nt1)] += 1
                            else :
                                dico_types_CWW.update({(nt1, nt2) : 1})
    for cle in dico_types_CWW.keys() :
        print(cle, dico_types_CWW[cle])
        
    
    with open("/media/coline/Maxtor/dico_new_200320_sim_par_branche_0.65_avec_CWW_2.pickle", 'rb') as fichier_dico_sim :
        mon_depickler_2 = pickle.Unpickler(fichier_dico_sim)
        dico_sim = mon_depickler_2.load()
   
        dico_types_CWW_super = {}
        for paire in dico_sim.keys() :
            
            with open("Nouvelles_donnees/fichier_%s_%s_4.pickle"%(paire[0][0], paire[0][1]), 'rb') as fichier_graphe :
                mon_depickler_g = pickle.Unpickler(fichier_graphe)
                graphe1 = mon_depickler_g.load()
                
                with open("Nouvelles_donnees/fichier_%s_%s_4.pickle"%(paire[1][0], paire[1][1]), 'rb') as fichier_graphe :
                    mon_depickler_g = pickle.Unpickler(fichier_graphe)
                    graphe2 = mon_depickler_g.load()
                
                    for u,v,data in dico_sim[paire]["graphe"].edges(data=True) :
                        if data["label"] == 'CWW' or data["label"] == 'CWWn' :
                            if graphe1.nodes[u[0]]["poids"] == 1  :
                                nt11 = dico_graphes[paire[0]].nodes[graphe1.nodes[u[0]]["position"][0]]["nt"]
                                nt12 = dico_graphes[paire[0]].nodes[graphe1.nodes[v[0]]["position"][0]]["nt"]
                            
                                nt21 = dico_graphes[paire[1]].nodes[graphe2.nodes[u[1]]["position"][0]]["nt"]
                                nt22 = dico_graphes[paire[1]].nodes[graphe2.nodes[v[1]]["position"][0]]["nt"]
                                
                                if ((nt11,nt12), (nt21, nt22)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt11,nt12), (nt21, nt22))] += 1
                                elif ((nt12,nt11), (nt21, nt22)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt12,nt11), (nt21, nt22))] += 1
                                elif ((nt11,nt12), (nt22, nt21)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt11,nt12), (nt22, nt21))] += 1
                                elif ((nt12,nt11), (nt22, nt21)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt12,nt11), (nt22, nt21))] += 1    
                                
                                elif ((nt21, nt22), (nt11,nt12)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt21, nt22), (nt11,nt12))] += 1
                                elif ((nt21, nt22), (nt12,nt11)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt21, nt22), (nt12,nt11))] += 1
                                elif ((nt22, nt21), (nt11,nt12)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt22, nt21), (nt11,nt12))] += 1
                                elif ((nt22, nt21), (nt12,nt11)) in dico_types_CWW_super.keys() :
                                    dico_types_CWW_super[((nt22, nt21), (nt12,nt11))] += 1     
                                else :
                                    dico_types_CWW_super.update({((nt11,nt12), (nt21, nt22)) : 1})
                                    
                            else :
                                compteur = 0
                                for i in range(graphe1.nodes[u[0]]["position"][0], graphe1.nodes[u[0]]["position"][1]+1) :
                                    nt11 = dico_graphes[paire[0]].nodes[i]["nt"]
                                    nt12 = dico_graphes[paire[0]].nodes[graphe1.nodes[v[0]]["position"][1]-compteur]["nt"]
                                    
                                    nt21 = dico_graphes[paire[1]].nodes[graphe2.nodes[u[1]]["position"][0]]["nt"]
                                    nt22 = dico_graphes[paire[1]].nodes[graphe2.nodes[v[1]]["position"][0]]["nt"]
                                    
                                    compteur += 1
                            
                    
    
if __name__ == '__main__':
    
    
    
    
    
    #donne_sim_score_rmsd((('4v67', 7), ('1vq8', 3)))
    #exit()
#     # 23S
#     liste_verte = [(('1vqo', 12), ('6nd6', 7)),(('4y4o', 40), ('4u4r', 8)),(('6eri', 9), ('4u4r', 8)), (('6hma', 15), ('1vq8', 15)), (('5dm6', 2), ('6eri', 1)), (('4u27', 5), ('1vqp', 19)),(('4u27', 5), ('3ccr', 4)),(('4u27', 5), ('3cc7', 4)),(('5dm6', 7), ('1vq8', 11)),(('5dm6', 7), ('3ccr', 12)),(('5dm6', 7), ('3ccl', 11)),(('5dm6', 7), ('3ccs', 12)),(('5dm6', 7), ('3cc7', 12)),(('5dm6', 7), ('3ccq', 12)),(('5dm6', 7), ('3cce', 9)),(('5dm6', 7), ('3ccu', 11)),(('4ybb', 8), ('1vq8', 11)),(('4ybb', 8), ('3ccr', 12)),(('4ybb', 8), ('3ccl', 11)),(('4ybb', 8), ('3ccs', 12)),(('4ybb', 8), ('3cc7', 12)),(('4ybb', 8), ('3ccq', 12)),(('4ybb', 8), ('3cce', 9)),(('4ybb', 8), ('3ccu', 11)),(('3cc2', 2), ('6hma', 1)),(('5dm6', 13), ('1vqp', 19)),(('5dm6', 13), ('3ccr', 4)),(('5dm6', 13), ('3cc7', 4)), (('6hma', 9), ('4u4r', 8)), (('4y4o', 25), ('1vq8', 5)),(('4y4o', 25), ('3cc2', 2)),(('6hma', 15), ('4u4r', 16)), (('1vqp', 19), ('6eri', 9)), (('3ccr', 4), ('6eri', 9)),(('6eri', 9), ('3cc7', 4)),(('5afi', 15), ('4u4r', 28)),(('6hma', 9), ('1vqp', 19)),(('6hma', 9), ('3ccr', 4)),(('6hma', 9), ('3cc7', 4)),(('1vq8', 11), ('6hma', 6)),(('3ccr', 12), ('6hma', 6)),(('6hma', 6), ('3ccl', 11)),(('6hma', 6), ('3ccs', 12)),(('6hma', 6), ('3cc7', 12)),(('6hma', 6), ('3ccq', 12)),(('6hma', 6), ('3cce', 9)),(('6hma', 6), ('3ccu', 11)),(('4y4o', 40), ('1vqp', 19)),(('4y4o', 40), ('3ccr', 4)),(('4y4o', 40), ('3cc7', 4)),(('4y4o', 3), ('3ccl', 11)),(('4y4o', 3), ('3cce', 9)),(('4v51', 22), ('3ccl', 11)),(('4v51', 22), ('3cce', 9)),(('5e81', 12), ('3ccl', 11)),(('5e81', 12), ('3cce', 9)),(('3ccl', 11), ('4v90', 13)),(('4v90', 13), ('3cce', 9)),(('1vqo', 12), ('4u4r', 15)),(('4ybb', 1), ('1vq8', 5)),(('4y4o', 3), ('3ccr', 12)),(('4y4o', 3), ('3ccs', 12)),(('4y4o', 3), ('3cc7', 12)),(('4y4o', 3), ('3ccu', 11)),(('4v51', 22), ('3ccr', 12)),(('4v51', 22), ('3ccs', 12)),(('4v51', 22), ('3cc7', 12)),(('4v51', 22), ('3ccu', 11)),(('5e81', 12), ('3ccr', 12)),(('5e81', 12), ('3ccs', 12)),(('5e81', 12), ('3cc7', 12)),(('5e81', 12), ('3ccu', 11)),(('3ccr', 12), ('4v90', 13)),(('4v90', 13), ('3ccs', 12)),(('4v90', 13), ('3cc7', 12)),(('4v90', 13), ('3ccu', 11)),(('4y4o', 3), ('3ccq', 12)),(('4v51', 22), ('3ccq', 12)),(('5e81', 12), ('3ccq', 12)),(('4v90', 13), ('3ccq', 12)),(('4ybb', 1), ('3cc2', 2)),(('4y4o', 3), ('1vq8', 11)),(('4v51', 22), ('1vq8', 11)),(('5e81', 12), ('1vq8', 11)),(('1vq8', 5), ('6hma', 1)),(('1vq8', 11), ('4v90', 13)), (('5dm6', 9), ('1vq8', 21)), (('6hma', 2), ('4u4r', 19)), (('4ybb', 6), ('1yhq', 24)), (('6hma', 15), ('3cc2', 6))]
#     # 16S
#     liste_verte.extend([(('5nwy', 17), ('4y4o', 4)), (('6ek0', 2), ('6az1', 1)), (('4u27', 54), ('4y4o', 22)), (('4ybb', 35), ('4u3u', 22))])
#     ## homologues cherches a la main avec bonne sim et rmsd
#     liste_verte.extend([(('4y4o', 13), ('4u3u', 22)),(('3cc2', 1), ('4u4r', 36)),(('1vqo', 14), ('4u4r', 36)),(('6eri', 3), ('4u3u', 2)),(('4ybb', 6), ('4u3u', 2)),(('4y4o', 9), ('4u4r', 11)),(('4v67', 23), ('4u4r', 11)),(('4v67', 23), ('6ek0', 6)),(('4y4o', 9), ('6ek0', 6)),(('4y4o', 8), ('3cc2', 1)),(('3t1y', 8), ('4u4r', 11)),(('3t1y', 8), ('6ek0', 6)),(('2qex', 19), ('4u3u', 2)),(('1yhq', 24), ('4u3u', 2)),(('4y4o', 13), ('6az1', 6)),(('4y4o', 8), ('1vqo', 14)),(('5ngm', 1), ('4u4r', 11)),(('6eri', 15), ('4u4r', 11)),(('5ngm', 1), ('6ek0', 6)),(('6eri', 15), ('6ek0', 6)),(('4ybb', 9), ('4u4r', 11)),(('5j7l', 10), ('4u4r', 11)),(('4u27', 32), ('4u4r', 11)),(('5j7l', 10), ('6ek0', 6)),(('4ybb', 9), ('6ek0', 6)),(('4u27', 32), ('6ek0', 6)),(('1vq8', 19), ('4u3u', 2)),(('4ybb', 35), ('6az1', 6)),(('2vqe', 5), ('4u4r', 11)),(('2vqe', 5), ('6ek0', 6))])
#        
#     ## homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim  a 0.75 et pas avec la rmsd a 2
#     liste_verte.extend([(('4y4o', 25),('1vq8', 5)),(('4y4o', 25),('3cc2', 2)),(('4ybb', 1),('1vq8', 5)),(('4ybb', 1),('3cc2', 2)), (('6hma', 1),('1vq8', 5)),(('6hma', 1),('3cc2', 2))])
#        
#     ## homologues du 05/03/20 : cherches dans les bleus totaux faux neg, faux pos et vrais pos
#     liste_verte.extend([(('4u27', 54), ('6eri', 14)), (('2r8s', 1), ('1l8v', 1)), (('4y4o', 51), ('4u4r', 21)), (('4ybb', 7), ('4u4r', 21)),(('5ngm', 4), ('4u4r', 21)),(('4y4o', 11), ('4u4r', 24)), (('4y4o', 51), ('6az1', 3)), (('4ybb', 7), ('6az1', 3)), (('5ngm', 4), ('6az1', 3)), (('1vqp', 19), ('4u4r', 8)), (('4u4u', 11), ('6az1', 1)), (('4y4o', 8), ('4u4r', 36)), (('3ccr', 4), ('4u4r', 8)),(('3cc7', 4), ('4u4r', 8)), (('4y4o', 11), ('6eri', 10)), (('1gid', 1), ('1l8v', 1))])
#        
#     # 23S
#     liste_rouge = [(('1vq8', 18), ('5dm6', 5)),(('5dm6', 5), ('3ccl', 14)),(('5dm6', 5), ('3cce', 12)),(('1vq8', 1), ('6eri', 2)), (('5wfs', 19), ('4u4r', 6)),(('4ybb', 12), ('4y4o', 35)),(('4ybb', 12), ('4wsd', 40)),(('4y4o', 34), ('2qex', 12)),(('4y4o', 28), ('5wfs', 6)),(('4v67', 43), ('2qex', 12)), (('4y4o', 50), ('6ek0', 4)), (('5dm6', 4), ('4ybb', 17)), (('4w2g', 52), ('6ek0', 10)), (('4y4o', 47), ('1vq8', 19)),(('4y4o', 47), ('1yhq', 24)),(('4y4o', 47), ('2qex', 19)),(('4y4o', 47), ('4u4r', 13)),(('4ybb', 13), ('1vqo', 8)),(('4w2f', 23), ('4u4r', 36)),(('1vqo', 8), ('5dm6', 1)),(('1vqo', 8), ('1mms', 1)),(('1vq8', 19), ('6eri', 1)),(('1yhq', 24), ('6eri', 1)),(('6eri', 1), ('2qex', 19)), (('5dm6', 14), ('4y4o', 8)),(('4ybb', 17), ('4y4o', 8)), (('4ybb', 30), ('2zjr', 15)),(('2zjr', 15), ('5wfs', 11)), (('2zjr', 15), ('6eri', 8)), (('4y4o', 12), ('2zjr', 15)), (('4w2g', 52), ('4y4o', 28)), (('4y4o', 28), ('4v67', 7))]
#     #liste_rouge.extend([(('4ybb', 12), ('4y4o', 38)),(('4ybb', 12), ('4ybb', 21)),(('4ybb', 12), ('6ek0', 7)),(('4v67', 7), ('4y4o', 38)),(('4v67', 7), ('4ybb', 21)),(('4v67', 7), ('6ek0', 7)),(('5wfs', 19), ('4y4o', 38)),(('5wfs', 19), ('4ybb', 21)),(('5wfs', 19), ('6ek0', 7))])
#     ## non homologues meme type rmsd 2-2.5, sim 0.7-0.75
#     liste_rouge.extend([(('4y4o', 23), ('6h4n', 24)),(('2zjr', 3), ('4u4r', 18)),(('4u27', 54), ('4y4o', 38)),(('4u27', 54), ('6ek0', 7)),(('4u27', 54), ('4u3u', 7)), (('4y4o', 23), ('4u27', 3)),(('4y4o', 23), ('5afi', 24))])
#        
#     ## les non homologues trouves parmi les 23S et cie qu'on met ensemble avec notre sim à 0.75 mais pas avec la rmsd à 2
#     liste_rouge.extend([(('4ybb', 12), ('6ek0', 10)),(('4w2g', 52), ('6ek0', 10)),(('2zjr', 3),('6ek0', 10)),(('2zjr', 3),('4u4r', 18)),(('4w2g', 52),('4u4r', 18)),(('4v67', 7), ('4u4r', 18)),(('5wfs', 19), ('4u4r', 18)),(('4v67', 7),('6ek0', 10)),(('5wfs', 19),('6ek0', 10)),(('4y4o', 23),('6h4n', 24)),(('4y4o', 23), ('4u27', 3)),(('4y4o', 23),('5afi', 24))])
#     #liste_rouge = []
#        
#     ## non-homologues du 05/03/20 : cherches dans les bleus totaux faux neg, faux pos et vrais pos (voir fichier excel)
#     liste_rouge.extend([(('4u27', 54), ('4ybb', 21)), (('6eri', 14), ('6eri', 15)),(('5ngm', 1), ('6eri', 14)),(('4u27', 54), ('6eri', 15)), (('6eri', 14), ('4u4r', 21)), (('4y4o', 23), ('1vq8', 12)),(('4y4o', 23), ('3cc2', 19)),(('4y4o', 23), ('3cpw', 14)), (('4y4o', 22), ('6ek0', 7)),(('4y4o', 22), ('4u3u', 7)),(('1yhq', 9), ('5afi', 24)),(('5afi', 24), ('2qex', 12)),(('4ybb', 37), ('4ybb', 1)),(('4ybb', 1), ('4y4o', 19)),(('4ybb', 1), ('1vq8', 16)),(('4ybb', 1), ('6h4n', 7)),(('4ybb', 1), ('5wfs', 19)),(('4ybb', 1), ('6eri', 16)),(('4ybb', 1), ('6qul', 4)),(('4ybb', 1), ('3ccm', 16)),(('4ybb', 1), ('4u4r', 6)),(('4ybb', 1), ('6ek0', 3)),(('4ybb', 12), ('1vq8', 16)),(('4w2g', 52), ('1vq8', 16)),(('4w2g', 52), ('3cd6', 14)),(('4w2g', 52), ('2qex', 20)),(('2zjr', 3), ('1vq8', 16)),(('4u27', 3), ('1yhq', 9)),(('4u27', 3), ('2qex', 12)),(('1vq8', 8), ('1vq8', 12)),(('1vq8', 8), ('3cc2', 19)),(('1vq8', 8), ('3cpw', 14)),(('5dm6', 4), ('6hma', 8)),(('4ybb', 12), ('6qul', 4)),(('4w2g', 52), ('6qul', 4)),(('6hma', 8), ('6hma', 14)),(('6hma', 8), ('1vq8', 18)),(('6hma', 8), ('3ccl', 14)),(('6hma', 8), ('6eri', 12)),(('6hma', 8), ('3cce', 12)),(('6hma', 8), ('4u4r', 1)),(('2zjr', 3), ('6qul', 4)),(('1vq8', 8), ('5afi', 24)),(('1vq8', 12), ('1yhq', 9)),(('1vq8', 12), ('2qex', 12)),(('1yhq', 9), ('3cc2', 19)),(('1yhq', 9), ('6h4n', 24)),(('1yhq', 9), ('3cpw', 14)),(('3cc2', 19), ('2qex', 12)),(('6h4n', 24), ('2qex', 12)),(('3cpw', 14), ('2qex', 12)),(('5wfs', 19), ('6qul', 4)), (('4y4o', 22), ('4y4o', 38)),(('4y4o', 22), ('4ybb', 21)), (('4u27', 3), ('1vq8', 8)),(('1vq8', 8), ('6h4n', 24)), (('4ybb', 37), ('4y4o', 25)),(('4ybb', 37), ('6hma', 1)),(('4y4o', 25), ('1vq8', 16)),(('4y4o', 25), ('3ccr', 17)),(('4y4o', 25), ('3cd6', 14)),(('4y4o', 25), ('4u4r', 6)),(('4y4o', 19), ('6hma', 1)),(('1vq8', 16), ('6hma', 1)),(('1yhq', 17), ('6hma', 1)),(('6h4n', 7), ('6hma', 1)),(('3ccr', 17), ('6hma', 1)),(('6hma', 1), ('6hma', 10)),(('6hma', 1), ('3cd6', 14)),(('6hma', 1), ('3ccs', 16)),(('6hma', 1), ('6eri', 16)),(('6hma', 1), ('3cc7', 17)),(('6hma', 1), ('6qul', 4)),(('6hma', 1), ('3ccq', 18)),(('6hma', 1), ('3ccm', 16)),(('6hma', 1), ('3ccu', 16)),(('6hma', 1), ('2qex', 20)),(('6hma', 1), ('4u4r', 6)),(('6hma', 1), ('6ek0', 3)),(('5dm6', 4), ('5dm6', 9)),(('5dm6', 4), ('5afi', 17)),(('5dm6', 9), ('6eri', 12)),(('4y4o', 28), ('6eri', 12)),(('5afi', 17), ('6eri', 12)),(('4ybb', 12), ('6ek0', 3)),(('4w2g', 52), ('1yhq', 17)),(('4w2g', 52), ('6hma', 10)),(('4w2g', 52), ('3cc7', 17)),(('4w2g', 52), ('6ek0', 3)),(('4y4o', 25), ('4y4o', 19)),(('4y4o', 25), ('1yhq', 17)),(('4y4o', 25), ('6hma', 10)),(('4y4o', 25), ('3ccs', 16)),(('4y4o', 25), ('6eri', 16)),(('4y4o', 25), ('3cc7', 17)),(('4y4o', 25), ('6qul', 4)),(('4y4o', 25), ('3ccq', 18)),(('4y4o', 25), ('3ccm', 16)),(('4y4o', 25), ('3ccu', 16)),(('4y4o', 25), ('6ek0', 3)),(('2zjr', 3), ('1yhq', 17)),(('2zjr', 3), ('6hma', 10)),(('2zjr', 3), ('3cc7', 17)),(('2zjr', 3), ('6ek0', 3)),(('5wfs', 19), ('6ek0', 3)),(('5dm6', 9), ('6hma', 14)),(('5dm6', 9), ('1vq8', 18)),(('5dm6', 9), ('3ccl', 14)),(('5dm6', 9), ('3cce', 12)),(('5dm6', 9), ('4u4r', 1)),(('4y4o', 28), ('6hma', 14)),(('4y4o', 28), ('1vq8', 18)),(('4y4o', 28), ('4u4r', 1)), (('4ybb', 12), ('1vq8', 3)),(('4ybb', 12), ('3cc2', 21)),(('4ybb', 12), ('3cpw', 3)), (('5dm6', 4), ('4y4o', 28)),(('4ybb', 12), ('6h4n', 7)),(('4w2g', 52), ('6h4n', 7)),(('4w2f', 37), ('4y4o', 28)),(('4y4o', 28), ('4y4o', 43)),(('4y4o', 28), ('3ccl', 14)),(('4y4o', 28), ('3cce', 12)),(('4y4o', 25), ('6h4n', 7)), (('2zjr', 3), ('6h4n', 7)), (('4u27', 54), ('5ngm', 1)), (('6eri', 14), ('6ek0', 7)),(('6eri', 14), ('4u3u', 7))])
   
       
    with open("groupes_homologues_60_2.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        groupes_homologues = mon_depickler.load()
        
        liste_elts = [('3cc2', 19), ('6eri', 17), ('1vq8', 12), ('3cpw', 14)]
        for elt in liste_elts : 
            with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(elt[0], elt[1]), 'rb') as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                graphe  = mon_depickler.load()   
                 
                seq = recup_sequences_autour_motif(elt[0], graphe, 4)
                compteur = 0
                for groupe in groupes_homologues :
                    if elt in groupe :
                        print(compteur)
                    compteur += 1
                #print(elt)
                print(seq[2])
                #print(seq[1])
                #print(seq[2])
        #exit()
        #pourcentage_liaison_near(['23S', '18S', '16S', 'Ribozyme', 'Riboswitch', 'SRP', '28S', '25S', 'Intron', 'arnt_16S_arnm', 'arnt_16S'])
        #exit()
#         liste = [(('5dm6', 3), ('4ybb', 5)), (('5dm6', 3), ('6hrm', 12)), (('5dm6', 4), ('4u4r', 1)), (('5dm6', 10), ('4v51', 49)), (('5dm6', 10), ('5ngm', 15)), (('5dm6', 10), ('6eri', 17)), (('5dm6', 9), ('1vq8', 21)), (('5dm6', 9), ('3cc2', 17)), (('4ybb', 6), ('6eri', 3)), (('4ybb', 14), ('4y4o', 15)), (('4ybb', 14), ('1vq8', 13)), (('4ybb', 14), ('3cc2', 14)), (('4ybb', 14), ('4v51', 5)), (('4ybb', 14), ('4v9f', 8)), (('4ybb', 14), ('6eri', 23)), (('4ybb', 14), ('2qex', 8)), (('4ybb', 14), ('5tbw', 7)), (('4ybb', 17), ('1vqo', 8)), (('4ybb', 19), ('4w2g', 53)), (('4ybb', 19), ('4y4o', 37)), (('4ybb', 19), ('6hma', 17)), (('4ybb', 19), ('3ccl', 22)), (('4ybb', 19), ('3cd6', 19)), (('4ybb', 19), ('3cce', 20)), (('4ybb', 20), ('5dm6', 5)), (('4u26', 22), ('5j5b', 15)), (('4ybb', 30), ('6ek0', 10)), (('4ybb', 30), ('4u4r', 18)), (('4ybb', 54), ('1vq8', 21)), (('4ybb', 54), ('3cc2', 17)), (('4ybb', 12), ('1vq8', 3)), (('4ybb', 12), ('1vq8', 10)), (('4ybb', 12), ('1vq8', 15)), (('4ybb', 12), ('3cc2', 4)), (('4ybb', 12), ('3cc2', 6)), (('4ybb', 12), ('3cc2', 21)), (('4ybb', 12), ('3cpw', 3)), (('4ybb', 12), ('2qex', 24)), (('4ybb', 12), ('6ek0', 7)), (('4ybb', 12), ('4u3u', 7)), (('4ybb', 12), ('6az1', 5)), (('4ybb', 12), ('4u4r', 18)), (('4ybb', 12), ('4faw', 1)), (('4ybb', 27), ('5dm6', 5)), (('4y4o', 23), ('4u27', 3)), (('4y4o', 23), ('5afi', 24)), (('4y4o', 23), ('6h4n', 24)), (('4y4o', 40), ('6hma', 9)), (('4y4o', 34), ('4v67', 43)), (('4y4o', 12), ('2zjr', 3)), (('4y4o', 12), ('4y4o', 18)), (('4y4o', 12), ('4u3u', 7)), (('4y4o', 12), ('6ek0', 10)), (('4y4o', 12), ('4u4r', 9)), (('4y4o', 12), ('4u4r', 18)), (('4w2g', 52), ('4u4r', 18)), (('4y4o', 15), ('5mdv', 21)), (('4y4o', 15), ('6qul', 18)), (('4y4o', 28), ('1vq8', 21)), (('4y4o', 28), ('3cc2', 17)), (('4y4o', 28), ('4y1n', 1)), (('4y4o', 25), ('1vq8', 5)), (('4y4o', 25), ('3cc2', 2)), (('4w2g', 53), ('6hrm', 13)), (('4y4o', 53), ('4v51', 49)), (('4y4o', 53), ('6ek0', 12)), (('5mdv', 21), ('1vq8', 13)), (('5mdv', 21), ('3cc2', 14)), (('5mdv', 21), ('4v51', 5)), (('5mdv', 21), ('4v9f', 8)), (('5mdv', 21), ('6eri', 23)), (('5mdv', 21), ('2qex', 8)), (('5mdv', 21), ('5tbw', 7)), (('4ybb', 13), ('4w2f', 23)), (('4ybb', 13), ('5dm6', 1)), (('4ybb', 13), ('6eri', 4)), (('5dm6', 6), ('1vq8', 14)), (('5dm6', 14), ('1vqo', 8)), (('4y4o', 37), ('6hrm', 13)), (('4w2f', 23), ('5d8h', 1)), (('4w2f', 23), ('1y39', 2)), (('6hma', 15), ('4y4o', 18)), (('6hma', 16), ('1vq8', 14)), (('6hma', 2), ('4v51', 49)), (('6hma', 2), ('5ngm', 15)), (('6hma', 2), ('6eri', 17)), (('6hma', 8), ('1vq8', 21)), (('6hma', 8), ('3cc2', 17)), (('6hma', 11), ('1vqo', 8)), (('4y4o', 8), ('4u4r', 36)), (('4w2f', 4), ('1vqo', 6)), (('4v51', 29), ('1vq8', 14)), (('4v51', 49), ('1vq8', 4)), (('4v51', 49), ('5wdt', 12)), (('2zjr', 3), ('1vq8', 3)), (('2zjr', 3), ('1vq8', 15)), (('2zjr', 3), ('3cc2', 6)), (('2zjr', 3), ('3cc2', 21)), (('2zjr', 3), ('4y4o', 18)), (('2zjr', 3), ('3ccq', 3)), (('2zjr', 3), ('3cpw', 3)), (('2zjr', 3), ('6eri', 8)), (('2zjr', 3), ('5wfs', 11)), (('2zjr', 3), ('2qex', 16)), (('2zjr', 3), ('6ek0', 7)), (('2zjr', 3), ('4u3u', 7)), (('2zjr', 3), ('6hrm', 14)), (('2zjr', 3), ('4u4r', 16)), (('2zjr', 3), ('4u4r', 18)), (('5e81', 12), ('6eri', 13)), (('4wsd', 4), ('1vq8', 14)), (('1vq8', 3), ('4v67', 7)), (('1vq8', 3), ('5wfs', 19)), (('1vq8', 3), ('6hrm', 8)), (('1vq8', 3), ('4u3u', 7)), (('1vq8', 3), ('1u9s', 1)), (('1vq8', 4), ('5ngm', 15)), (('1vq8', 4), ('6eri', 17)), (('1vq8', 5), ('6hma', 1)), (('1vq8', 5), ('6hrm', 2)), (('1vq8', 5), ('6ek0', 10)), (('1vqo', 8), ('6hrm', 20)), (('1vqo', 8), ('6ek0', 1)), (('1vqo', 8), ('4u4r', 12)), (('1vq8', 10), ('4v67', 7)), (('1vq8', 10), ('5wfs', 19)), (('1vq8', 10), ('1u9s', 1)), (('1vqo', 12), ('4u4r', 15)), (('1vqo', 14), ('4u4r', 36)), (('1vq8', 13), ('4v51', 5)), (('1vq8', 13), ('5tbw', 7)), (('1vq8', 14), ('4v67', 5)), (('1vq8', 15), ('4v67', 7)), (('1vq8', 15), ('5wfs', 19)), (('1vq8', 15), ('6ek0', 7)), (('1vq8', 15), ('4u3u', 7)), (('1vq8', 15), ('4u27', 54)), (('1vq8', 15), ('1u9s', 1)), (('1vq8', 15), ('6ek0', 10)), (('1vq8', 15), ('4u4r', 9)), (('1vq8', 15), ('4u4r', 18)), (('1vqo', 21), ('4u4r', 28)), (('1vq8', 19), ('6eri', 3)), (('1vq8', 21), ('5afi', 17)), (('1vq8', 21), ('6eri', 11)), (('1vq8', 21), ('6hrm', 5)), (('1vq8', 21), ('6ek0', 4)), (('1vq8', 21), ('4u4r', 4)), (('5dm6', 5), ('5afi', 15)), (('5dm6', 5), ('4v67', 43)), (('5dm6', 5), ('6hma', 13)), (('1vy7', 20), ('6i7v', 27)), (('1vy7', 20), ('4u4r', 13)), (('3cc2', 2), ('6hma', 1)), (('3cc2', 2), ('6hrm', 2)), (('3cc2', 2), ('6ek0', 10)), (('3cc2', 4), ('4v67', 7)), (('3cc2', 4), ('5wfs', 19)), (('3cc2', 4), ('1u9s', 1)), (('3cc2', 6), ('4v67', 7)), (('3cc2', 6), ('5wfs', 19)), (('3cc2', 6), ('6ek0', 7)), (('3cc2', 6), ('4u3u', 7)), (('3cc2', 6), ('4u27', 54)), (('3cc2', 6), ('1u9s', 1)), (('3cc2', 6), ('6ek0', 10)), (('3cc2', 6), ('4u4r', 9)), (('3cc2', 6), ('4u4r', 18)), (('3cc2', 14), ('4v51', 5)), (('3cc2', 14), ('5tbw', 7)), (('3cc2', 17), ('5afi', 17)), (('3cc2', 17), ('6eri', 11)), (('3cc2', 17), ('6hrm', 5)), (('3cc2', 17), ('6ek0', 4)), (('3cc2', 17), ('4u4r', 4)), (('3cc2', 21), ('4v67', 7)), (('3cc2', 21), ('5wfs', 19)), (('3cc2', 21), ('6hrm', 8)), (('3cc2', 21), ('4u3u', 7)), (('3cc2', 21), ('1u9s', 1)), (('1yhq', 24), ('6eri', 3)), (('4y4o', 18), ('3ccq', 3)), (('4y4o', 18), ('6eri', 8)), (('4y4o', 18), ('5wfs', 11)), (('4y4o', 18), ('6hrm', 14)), (('4y4o', 18), ('4u4r', 16)), (('6h4n', 7), ('2qex', 16)), (('6h4n', 7), ('6hrm', 8)), (('6h4n', 7), ('4u27', 54)), (('6h4n', 7), ('6eri', 14)), (('6h4n', 7), ('3mut', 1)), (('3cc2', 1), ('4u4r', 36)), (('1vqo', 2), ('6ek0', 5)), (('5dm6', 15), ('6i7v', 27)), (('4v67', 7), ('3cpw', 3)), (('4v67', 7), ('2qex', 24)), (('4v67', 7), ('6ek0', 7)), (('4v67', 7), ('4u3u', 7)), (('4v67', 7), ('6az1', 5)), (('4v67', 7), ('4u4r', 18)), (('4v67', 7), ('4faw', 1)), (('1vqo', 6), ('6hma', 4)), (('1vqo', 6), ('6eri', 20)), (('4v51', 5), ('4v9f', 8)), (('4v51', 5), ('2qex', 8)), (('4v9f', 8), ('5tbw', 7)), (('6nd6', 7), ('4u4r', 15)), (('3ccr', 4), ('4u4r', 8)), (('6hma', 1), ('4faw', 1)), (('3ccq', 3), ('5ngm', 4)), (('3ccq', 3), ('2nz4', 3)), (('3ccq', 3), ('1u9s', 1)), (('6hma', 17), ('6hrm', 13)), (('6hma', 7), ('5ngm', 3)), (('6hma', 4), ('4u4q', 28)), (('3cpw', 3), ('5wfs', 19)), (('3cpw', 3), ('6hrm', 8)), (('3cpw', 3), ('4u3u', 7)), (('3cpw', 3), ('1u9s', 1)), (('3ccl', 22), ('6ek0', 5)), (('6eri', 1), ('6hrm', 25)), (('6eri', 1), ('4ybb', 4)), (('6eri', 1), ('4y4o', 54)), (('6eri', 1), ('6eri', 7)), (('6eri', 8), ('4u3u', 7)), (('6eri', 8), ('6ek0', 10)), (('6eri', 8), ('4u4r', 9)), (('6eri', 8), ('4u4r', 18)), (('3cd6', 19), ('6ek0', 5)), (('5wfs', 11), ('4u3u', 7)), (('5wfs', 11), ('6ek0', 10)), (('5wfs', 11), ('4u4r', 9)), (('5wfs', 11), ('4u4r', 18)), (('5wfs', 19), ('2qex', 24)), (('5wfs', 19), ('6ek0', 7)), (('5wfs', 19), ('4u3u', 7)), (('5wfs', 19), ('6az1', 5)), (('5wfs', 19), ('4u4r', 18)), (('5wfs', 19), ('4faw', 1)), (('6eri', 23), ('6qul', 18)), (('6eri', 3), ('6hrm', 15)), (('6eri', 3), ('4u3u', 2)), (('3cc7', 4), ('4u4r', 8)), (('6eri', 4), ('1y39', 2)), (('3cce', 20), ('6ek0', 5)), (('3ccm', 21), ('6ek0', 5)), (('1vqo', 16), ('2qex', 16)), (('2qex', 8), ('5tbw', 7)), (('2qex', 16), ('6hrm', 8)), (('2qex', 16), ('4u27', 54)), (('2qex', 16), ('6eri', 14)), (('2qex', 16), ('6ek0', 10)), (('2qex', 24), ('1u9s', 1)), (('6hrm', 1), ('6hrm', 8)), (('6hrm', 1), ('4v67', 23)), (('6hrm', 1), ('6eri', 14)), (('6hrm', 8), ('6ek0', 7)), (('6hrm', 8), ('4u3u', 7)), (('6hrm', 8), ('5ngm', 1)), (('6hrm', 8), ('6eri', 15)), (('6hrm', 8), ('3mut', 1)), (('6hrm', 8), ('4u4r', 18)), (('6hrm', 11), ('6eri', 19)), (('6hrm', 24), ('4woi', 62)), (('6hrm', 25), ('4u4r', 7)), (('6hrm', 26), ('4y4o', 11)), (('6hrm', 27), ('4u4r', 10)), (('6hrm', 27), ('6az1', 2)), (('6hrm', 27), ('4v7l', 15)), (('6ek0', 7), ('6eri', 17)), (('6ek0', 7), ('6eri', 14)), (('6ek0', 7), ('1u9s', 1)), (('6ek0', 7), ('6ek0', 10)), (('6ek0', 7), ('5nwy', 12)), (('4u4r', 7), ('6eri', 7)), (('4u3u', 7), ('6eri', 14)), (('4u3u', 7), ('6hrm', 14)), (('4u3u', 7), ('1u9s', 1)), (('4u3u', 7), ('6ek0', 10)), (('4u3u', 7), ('4u4r', 16)), (('4u4r', 21), ('6az1', 3)), (('4u4r', 17), ('5t2a', 8)), (('4u4r', 17), ('4wsd', 1)), (('4u4r', 17), ('3t1y', 6)), (('4u4r', 17), ('4woi', 43)), (('4u4r', 17), ('4v67', 25)), (('4u4r', 17), ('4v90', 2)), (('4u4r', 17), ('4v8b', 31)), (('4u4r', 17), ('4v9h', 14)), (('4u4r', 17), ('5e81', 1)), (('4u4r', 10), ('4v7l', 15)), (('4u4r', 10), ('6h4n', 17)), (('4u4r', 11), ('4v67', 23)), (('4u3u', 22), ('4ybb', 35)), (('6ek0', 6), ('4v67', 23)), (('6az1', 2), ('4v7l', 15)), (('6az1', 2), ('3t1y', 10)), (('6az1', 2), ('6h4n', 17)), (('6az1', 5), ('1u9s', 1)), (('6az1', 5), ('6ek0', 10)), (('6az1', 5), ('4faw', 1)), (('4ybb', 46), ('6eri', 17)), (('4ybb', 2), ('4v7l', 15)), (('4ybb', 51), ('6eri', 19)), (('4u27', 42), ('4y4o', 26)), (('4ybb', 7), ('4woi', 62)), (('4u27', 32), ('4v67', 23)), (('4u27', 54), ('3mut', 1)), (('4y4o', 58), ('5ngm', 15)), (('4y4o', 58), ('4w2f', 40)), (('4y4o', 58), ('5mdv', 30)), (('4y4o', 58), ('5e81', 27)), (('4y4o', 58), ('4wqf', 14)), (('4y4o', 58), ('5f8k', 15)), (('4y4o', 58), ('5wdt', 15)), (('4y4o', 58), ('5nwy', 12)), (('4y4o', 9), ('4v67', 23)), (('4y4o', 31), ('4v7l', 15)), (('4y4o', 11), ('6eri', 10)), (('4ybb', 9), ('4v67', 23)), (('5j7l', 10), ('4v67', 23)), (('2vqe', 5), ('4v67', 23)), (('5ngm', 4), ('2nz4', 3)), (('5ngm', 4), ('4u4r', 28)), (('4wsd', 1), ('4u26', 1)), (('5ndk', 13), ('4v7l', 15)), (('4v8d', 26), ('4v7l', 15)), (('3t1y', 8), ('4v67', 23)), (('4v7l', 15), ('4u27', 2)), (('4v7l', 15), ('5f8k', 24)), (('4v7l', 15), ('3t1y', 10)), (('4v7l', 15), ('4v67', 6)), (('4v7l', 15), ('6h4n', 17)), (('4v7l', 15), ('5ngm', 18)), (('4v7l', 15), ('5e81', 10)), (('4v7l', 15), ('6eri', 18)), (('5f8k', 24), ('6ek0', 4)), (('5f8k', 24), ('4u4r', 4)), (('3t1y', 6), ('4u26', 1)), (('3t1y', 10), ('4v90', 11)), (('4woi', 43), ('4u26', 1)), (('4v67', 23), ('5ngm', 1)), (('4v67', 23), ('6eri', 15)), (('4v67', 25), ('4u26', 1)), (('4v90', 2), ('4u26', 1)), (('4v8b', 31), ('4u26', 1)), (('4v9h', 14), ('4u26', 1)), (('4u26', 1), ('5e81', 1)), (('5ngm', 1), ('6eri', 14)), (('5ngm', 9), ('6eri', 19)), (('5ngm', 15), ('4w2f', 40)), (('5ngm', 15), ('5mdv', 30)), (('5ngm', 15), ('5e81', 27)), (('5ngm', 15), ('5f8k', 15)), (('5ngm', 15), ('5wdt', 15)), (('5ngm', 15), ('5nwy', 12)), (('5ngm', 18), ('6ek0', 4)), (('5ngm', 18), ('4u4r', 4)), (('6eri', 17), ('4w2f', 40)), (('6eri', 17), ('5mdv', 30)), (('6eri', 17), ('5e81', 27)), (('6eri', 17), ('4wqf', 14)), (('6eri', 17), ('5f8k', 15)), (('6eri', 17), ('4woi', 62)), (('6eri', 17), ('5wdt', 15)), (('6eri', 17), ('5nwy', 12)), (('6eri', 14), ('6eri', 15)), (('6eri', 14), ('3mut', 1)), (('6eri', 14), ('4u4r', 18)), (('5e81', 10), ('6ek0', 4)), (('5e81', 10), ('4u4r', 4)), (('6hrm', 2), ('4faw', 1)), (('6hrm', 14), ('6ek0', 10)), (('6hrm', 14), ('4u4r', 9)), (('6hrm', 14), ('4u4r', 18)), (('2r8s', 1), ('1gid', 1)), (('2nz4', 3), ('1u9s', 1)), (('1u9s', 1), ('3mut', 1)), (('1u9s', 1), ('4u4r', 18)), (('1u9s', 1), ('4faw', 1)), (('4l81', 2), ('4aob', 1)), (('4aob', 1), ('5fkf', 1)), (('4aob', 1), ('5fkd', 1)), (('6ek0', 10), ('4u4r', 16)), (('6ek0', 10), ('4u4r', 18)), (('6ek0', 10), ('4faw', 1)), (('4u4r', 16), ('4u4r', 9)), (('4u4r', 16), ('4u4r', 18)), (('4u4r', 18), ('4faw', 1)), (('1vy5', 21), ('5ib7', 31)), (('4w2f', 40), ('4wqf', 14)), (('4wqf', 14), ('4woi', 62)), (('4wqf', 14), ('5nwy', 12))]
#             
#         with open("csv_sim_modifie_0.65.csv", 'w', newline='') as fichier_csv :
#             csvwriter = csv.writer(fichier_csv)
#             csvwriter.writerow(["Paire", "Type1", "Type2", "Score", "RMSD", "Sim1", "Sim2", "Homologues", "Liste"])
#                 
#             res = donne_sim_score_rmsd_plusieurs(liste)
#                 
#             for elt in res :
#                 homol = False
#                 for groupe in groupes_homologues :
#                     if elt[0][0] in groupe and elt[0][1] in groupe :
#                         homol = True
#                 if homol :
#                     csvwriter.writerow([elt[0], elt[1], elt[2], elt[3], elt[4], elt[5], elt[6], homol])
#                 elif elt[0] in liste_verte or (elt[0][1], elt[0][0]) in liste_verte : 
#                     csvwriter.writerow([elt[0], elt[1], elt[2], elt[3], elt[4], elt[5], elt[6], homol,"vert"])
#                 elif elt[0] in liste_rouge or (elt[0][1], elt[0][0]) in liste_rouge : 
#                     csvwriter.writerow([elt[0], elt[1], elt[2], elt[3], elt[4], elt[5], elt[6], homol,"rouge"])
#                 else : 
#                     csvwriter.writerow([elt[0], elt[1], elt[2], elt[3], elt[4], elt[5], elt[6],  homol,"bleu"])
#                      
#              
#              
#     exit()
    types = ['23S', '18S', '16S', 'Ribozyme', 'Riboswitch', 'SRP', '28S', '25S', 'Intron', 'arnt_16S_arnm', 'arnt_16S']
    typ = [['23S', '25S', '28S'], ['16S', '18S'], 'Ribozyme', 'Riboswitch', 'Intron', 'arnt_16S_arnm', 'arnt_16S']
    
    #distrib_types_liaisons_CWW(types)
    #exit()
    #nts_modifies(types)
    #exit()
    chaines_incomp = []
    for t in typ : 
        chaines_incomp.extend(chaines_incompletes(t))
    print(chaines_incomp)
#     
#     with open(PATH_MMCIF+"fichiers_rmsd_taille_%d_que_carbone_3_new_data_new.pickle"%4, 'rb') as fichier_rmsd :
#         mon_depickler = pickle.Unpickler(fichier_rmsd)
#         rmsd = mon_depickler.load()
#          
#         compter_normal = 0
#         compter_pas_normal = 0
#         for cle in rmsd.keys() :
#             if rmsd[cle] == None :
#                 if (cle[0].split("_")[1], int(cle[0].split("_")[2])) not in chaines_incomp and (cle[1].split("_")[1], int(cle[1].split("_")[2])) not in chaines_incomp :
#                     if (cle[0].split("_")[1], int(cle[0].split("_")[2])) not in [('6az1', 5), ('6ek0', 9), ('4wqf', 14), ('6ek0', 1), ('6ek0', 8), ('6ek0', 9), ('6az1', 4), ('5wdt', 11), ('6gyv', 1), ('4y4o', 18), ('5dm6', 6)] and (cle[1].split("_")[1], int(cle[1].split("_")[2])) not in [('6az1', 5), ('6ek0', 9), ('4wqf', 14), ('6ek0', 1), ('6ek0', 8), ('6ek0', 9), ('6az1', 4), ('5wdt', 11), ('6gyv', 1), ('4y4o', 18), ('5dm6', 6)] :
#                         print(cle)
#                         compter_pas_normal += 1
#                 else :
#                     compter_normal += 1
#         print(compter_normal)
#         print(compter_pas_normal)
#         
#         
#     liste = [('6az1', 5), ('6ek0', 9), ('4wqf', 14), ('6ek0', 1), ('6ek0', 8), ('6ek0', 9), ('6az1', 4), ('5wdt', 11), ('6gyv', 1), ('4y4o', 18), ('5dm6', 6)]
#     for elt in liste :
#         print(renvoi_num_ARN(elt))
#     exit()
    
    #traitement_aln_global_inter_types(typ, "min", 2)
    #exit()
#     with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_030220.pickle", 'rb') as fichier_dico_sim :
#         mon_depickler_2 = pickle.Unpickler(fichier_dico_sim)
#         dico_sim = mon_depickler_2.load()
        
    #     distrib_sim_rmsd_inter_types(typ[0], typ[1]) ## distrib sim/rmsd 16S/23S
#         traitement_aln_global_inter_types(typ[0], typ[1], "min") ## creations des fichiers de stockage des scores et des RMSD des paires 23S/16S et distrib RMSD/score
        #recup_certaines_aretes_inter_types(types, typ[0], typ[1], 0, 4) ## creations des fichiers gephi du sous-graphe induit par les aretes inf à 4 en RMSD pour les 23S/16S
        #recup_certaines_aretes_plusieurs_types(types, typ, 0, 4) ## creations des fichiers gephi du sous-graphe induit par les aretes inf à 4 en RMSD pour tous les inter-types
        #exit()
        #liste_aretes_qui_ny_sont_pas = recup_groupes_non_homologues(types, typ[0], 70, 2) ## creations des fichiers gephi du sous-graphe induit par les aretes sup à 4 en rmsd et inf à 70 en score (grosse patate) des 23S
        #exit()
#         for i in np.arange(10, 101, 10) :
#             traitement_aln_global(['16S', '18S'], dico_sim,i, 0.6, "min", [])
        #exit()
#         with open("csv_diff_clustering_sim_rmsd_new_120320.csv", 'a', newline="") as fichier_csv :  
#             csvwriter = csv.writer(fichier_csv)
#             csvwriter.writerow(["Seuil sim", "Seuil eps", "Indice Jaccard", "Infos"])
#             jaccard = []
#             for j in np.arange(0.6,0.8, 0.05) :
#                 for i in np.arange(0.1, 1.1, 0.1) :
#                     jaccard = recup_groupes_homologues(types, ['23S', '25S', '28S'], 60, 2, j, i)
#                     csvwriter.writerow([j, i, jaccard[0], jaccard[1]])
#             print(jaccard)
# #         for elt in jaccard :
# #             print(elt[0], elt[1][0], elt[1][1], elt[1][4], elt[1][5])
#         print(jaccard)
        #exit()
        
        ## distrib sim/RMSD sans les homologues presumes avec un seuil de 80 en score et 2.5 en rmsd 
#         groupes_homologues = []
#         for groupe in typ :
#             homologues,num_pbs = recup_groupes_homologues(types,groupe , 55.5, 80, 2.5, 4.62, 55.5, 0.034, 0.65)
#             #distrib_sim_rmsd(groupe, homologues)
#             groupes_homologues.extend(homologues)
#         print(len(groupes_homologues))
# #         
#  
#         with open("groupes_homologues_0_2.pickle", 'wb') as fichier_pickle :
#             mon_pickler = pickle.Pickler(fichier_pickle)
#             mon_pickler.dump(groupes_homologues)


        #distrib_sim_rmsd(typ, groupes_homologues, 50, 0)
    liste_vrais_positifs = distrib_sim_rmsd_nombre(typ, groupes_homologues, 50, 0, 32, 20)
        #recup_groupes_vrais_positifs(liste_vrais_positifs, groupes_homologues, 0.75, 2.5)
    exit()
        #distrib_3D_sim_rmsd_score(typ, groupes_homologues)
        
    #     for t in typ :
    #         homologues,num_pbs = recup_groupes_homologues(types,t , 55.5, 80, 2.5, 4.62, 55.5, 0.034, 0.65)
        
        
        #homologues,num_pbs = recup_groupes_homologues(types,["16S", "18S"] , 55.5, 60, 2, 4.62, 55.5, 0.034, 0.65)
        #traitement_aln_global(["16S", "18S"], dico_sim, 30, 0.6, "min")
        
        ## alignement global des sequences, affichage des distrib rmsd/score et des distrib sim/RMSD de tous les types 
        #alignement_global_plusieurs_types(typ)
        #traitement_aln_global_inter_types(typ, "min")
        #distrib_sim_rmsd_plusieurs_types(typ)
        
        
        
    #exit()
        
        #creation_fichier_input_cdhit("23S")
        #traitement_res_cdhit("23S")
        
        #alignement_global_deux_types(['23S', '25S', '28S'], ['16S', '18S'])
        
        
        #distrib_3D_sim_rmsd_score(['23S', '25S', '28S'])
        #distrib_sim_rmsd(['16S', '18S'])
        #exit()
        
        
#     types_arn = ["18S", "16S","Ribozyme", "Riboswitch", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"]
#   
# 
#     
#     homologues,num_pbs = recup_groupes_homologues(types,["23S", "25S", "28S"] , 55.5, 60, 2.5, 4.62, 55.5, 0.034, 0.65)
#     #traitement_aln_global(["16S", "18S"], dico_sim, 30, 0.6, "min")
#     exit()
# #         
# #     for typ in types_arn :
# #         alignement_global(["23S","25S","28S"])
#     
#     #recup_groupes_homologues(typ, 60, 80, 2.5, 4.62)
#     groupes_homologues = []
#     compter = []
#     seuil1_score = 55.5
#     seuil2_score = 80
#     seuil1_rmsd = 2.5
#     seuil2_rmsd = 4.62
# 
#     seuil3_score = 55.5
#     a_droite = 0.034 
#     b_droite = 0.65
#     taille_pas_bon = 0
#     for t in typ :
#         homologues,num_pbs = recup_groupes_homologues(types, t, seuil1_score, seuil2_score, seuil1_rmsd, seuil2_rmsd, seuil3_score, a_droite, b_droite)
#         groupes_homologues.extend(homologues)
#         print(num_pbs)
#         print(taille_pas_bon)
#         for elt in num_pbs :
#             if taille_pas_bon > 0 :
#                 compter.append(elt+taille_pas_bon)
#             else :
#                 compter.append(elt)
#         taille_pas_bon = len(groupes_homologues)
#     print(groupes_homologues)
#     print(compter)
#     print(len(groupes_homologues))
#     nb_groupes = 0
#     for groupe in groupes_homologues :
#         if len(groupe) > 1 :
#             nb_groupes += 1
#     print(nb_groupes)
#     #with open("groupes_homologues_version_5fev_seuil_%s.pickle"%(str(seuil3_score)+"_"+str(a_droite)+"_"+str(b_droite)), 'wb') as fichier_hom :
#     with open("groupes_homologues_version_5fev_seuil_score_%s_rmsd_%s.pickle"%(str(seuil1_score)+"_"+str(seuil2_score), str(seuil1_rmsd)+"_"+str(seuil2_rmsd)), 'wb') as fichier_hom :
#         mon_pickler = pickle.Pickler(fichier_hom)
#         mon_pickler.dump(groupes_homologues)
#      
#     #with open("groupes_homologues_pas_bon_version_5fev_seuil_%s.pickle"%(str(seuil3_score)+"_"+str(a_droite)+"_"+str(b_droite)), 'wb') as fichier_hom_2 :
#     with open("groupes_homologues_pas_bon_version_5fev_seuil_score_%s_rmsd_%s.pickle"%(str(seuil1_score)+"_"+str(seuil2_score), str(seuil1_rmsd)+"_"+str(seuil2_rmsd)), 'wb') as fichier_hom_2 :
#         mon_pickler = pickle.Pickler(fichier_hom_2)
#         mon_pickler.dump(compter)  
#     
#     exit()
#     chaines_incompletes(types)
#     
#     with open("groupes_homologues_version_5fev_seuil_%s.pickle"%(str(seuil3_score)+"_"+str(a_droite)+"_"+str(b_droite)), 'rb') as fichier_hom :
#         mon_depickler = pickle.Unpickler(fichier_hom)
#         groupes_homologues_1 = mon_depickler.load()
#         
#     with open("groupes_homologues_version_5fev_seuil_score_%s_rmsd_%s.pickle"%(str(seuil1_score)+"_"+str(seuil2_score), str(seuil1_rmsd)+"_"+ str(seuil2_rmsd)), 'rb') as fichier_hom_2 :
#         mon_depickler_2 = pickle.Unpickler(fichier_hom_2)
#         groupes_homologues_2 = mon_depickler_2.load()
#     
#     
#     print(calcul_sim_jaccard(groupes_homologues_1, groupes_homologues_2))
#     #exit()
#     
#     #distrib_sim_rmsd(typ)
#     #exit()
#     typ = ['23S', '25S', '28S']
#     #for typ in types :
#     with open("Nouvelles_donnees/alignements_%s/fichier_nb_pos_sim_en_dessous_de_100.csv"%typ, 'w', newline='') as fichier_csv_ecriture :
#             csvwriter = csv.writer(fichier_csv_ecriture)
#             csvwriter.writerow(["seuil", "nb en dessous de 100"])
#             for seuil_pos in np.arange(30,40,10) :
                #for seuil_sim in np.arange(0,1.1,0.1) :
                #liste_en_dessous_de_100, liste_en_dessous_de_100_sim, liste_en_dessous_de_100_pos, liste_en_dessous_de_100_rmsd = traitement_aln_global(typ, dico_sim, seuil_pos, 0.7, "min")
    
        #traitement_aln_global_inter_types(['23S', '25S', '28S'], ['16S', '18S'], "min")
#                 liste_en_dessous_de_100.insert(0,len(liste_en_dessous_de_100))
#                 liste_en_dessous_de_100.insert(0,seuil_pos)
#                 liste_en_dessous_de_100_sim.insert(0," ")
#                 liste_en_dessous_de_100_sim.insert(0, " ")
#                 liste_en_dessous_de_100_pos.insert(0, " ")
#                 liste_en_dessous_de_100_pos.insert(0, " ")
#                 liste_en_dessous_de_100_rmsd.insert(0, " ")
#                 liste_en_dessous_de_100_rmsd.insert(0, " ")
#                 csvwriter.writerow(liste_en_dessous_de_100)
#                 csvwriter.writerow(liste_en_dessous_de_100_sim)
#                 csvwriter.writerow(liste_en_dessous_de_100_pos)
#                 csvwriter.writerow(liste_en_dessous_de_100_rmsd)
                
    
    #rechercher_homologues_par_blast("16S")
    #definit_groupes_homologues_par_blast("16S")
    
    #stocker_resolution()
    #stocker_resolution_qui_manque()
    
    #recup_type_molecule("Graphs_carnaval/7left.pickle", "fichier_molecules_7left.csv")

#     types_arn = ["23S", "16S", "18S", "25S", "28S", '5S', "LSU-alpha", "LSU-beta", "Riboswitch", "Ribozyme", "Intron", "ARNm"]
#     for typ in types_arn :
#         print(typ)
# #         #groupes_ARN(typ, "fichier_molecules_7left.csv", "Graphs_carnaval/7left.pickle", "motif_7left")
#         with open(NEW_EXTENSION_PATH_TAILLE+"%s/liste_representant_%s_res_3a.pickle"%("motif_7left", typ), 'rb') as fichier_sortie :
#             mon_depickler = pickle.Unpickler(fichier_sortie)
#             groupes = mon_depickler.load()
#             print(len(groupes))

#         #filtrage_sequences_identiques(typ, "motif_7left")
#     entre_plusieurs_chaines("Graphs_carnaval/7left.pickle")
    
#     for elt in types_arn :
#         recherche_representant(elt, "motif_7left")
#         distrib_resolution(elt, "motif_7left")
        
    
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_5dm6_8.pickle", "rb") as fichier_extension :
#                         mon_depickler_ext = pickle.Unpickler(fichier_extension)
#                         extension3 = mon_depickler_ext.load()
#     
#     seq11,seq12, seq13 = recup_sequences_autour_motif('5dm6', extension3, 3)
#     print(seq12)
    
    #types_arn = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"]
    #for typ in types_arn :
        #distrib_resolution(typ)
    #rechercher_homologues_par_blast("18S")
    #definit_groupes_homologues_par_blast("18S")
#     
#     with open("groupes_16S_homologues.pickle", "rb") as fichier_homologues :
#         mon_depickler = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler.load()
#         
#         with open("Nouvelles_donnees/liste_representant_16S_res_3a.pickle", 'rb') as fichier_repr :
#             mon_depickler = pickle.Unpickler(fichier_repr)
#             liste_representant = mon_depickler.load() 
#         
#         for groupe in groupes_homologues :
#             new_groupe = []
#             for elt in groupe :'Kluyveromyces lactis'
#                 if elt in liste_representant :
#                     new_groupe.append(elt)
#             if len(new_groupe) < 5 and len(new_groupe) > 2:
#                 for elt in new_groupe :
#                     print(elt)
#                     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_%s_%s.pickle"%(elt[0], elt[1]), "rb") as fichier_extension :
#                                 mon_depickler_ext = pickle.Unpickler(fichier_extension)
#                                 extension1 = mon_depickler_ext.load()
#                                 
#                                 seq11,seq12, seq13 = recup_sequences_autour_motif(elt[0], extension1, 30)
#                                 print(seq11)
#                                 print(seq12)
#                                 print(seq13)
#                 break
                         
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_4y4p_20.pickle", "rb") as fichier_extension :
#                         mon_depickler_ext = pickle.Unpickler(fichier_extension)
#                         extension2 = mon_depickler_ext.load()
#                          
#     with open(NEW_EXTENSION_PATH_TAILLE+"fichier_4y4o_39.pickle", "rb") as fichier_extension :
#                         mon_depickler_ext = pickle.Unpickler(fichier_extension)
#                         extension3 = mon_depickler_ext.load()
     
#     seq11,seq12, seq13 = recup_sequences_autour_motif('4y4o', extension1, 30)
#     seq21,seq22, seq23 = recup_sequences_autour_motif('4y4p', extension2, 30)
#     seq31,seq32, seq33 = recup_sequences_autour_motif('4y4o', extension3, 30)
#     
#     print(seq11)
#     print(seq21)
#     print(seq31)
#     
#     print(seq12)
#     print(seq22)
#     print(seq32)
#     
#     print(seq13)
#     print(seq23)
#     print(seq33)
    
    
#     print(seq11)
#     print(seq21)
#     print(seq31)
    
    #recherche_representant("arnt_16s_arnm")
    #stocker_resolution_qui_manque()
    #recherche_representant("arnt_16s")
    #filtrage_sequences_identiques("arnt_16s_arnm")
    #verif_filtrage("18S")
    #groupes_ARN("23S")
    #entre_plusieurs_chaines()
    #chercher_liaison_entre_chaines_par_sequence()
    
    
#     with open("new_groupes_homologues_28S.pickle", 'rb') as fichier_homologues :
#         mon_depickler = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler.load()
#         
#         print(len(groupes_homologues))
#         
#     with open("groupes_28S_homologues_sequences.pickle", 'rb') as fichier_homologues :
#         mon_depickler = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler.load()
#         
#         print(len(groupes_homologues))
# #         
# #         print(groupes_identiques)
# #     groupes_ARN_creation_homologues_et_identiques("18S")
#     #verif_homologues_version_2("28S", "Nouvelles_donnees/Groupes_homologues/28S")
#     #regrouper_par_homologues_critere_sequence("18S")
# #     types_oublies = ['4.5S', '4.8S']
#     types = ['4.5S', '4.8S', '25S', '18S', 'ARNt', 'Ribozyme', 'Riboswitch', 'SRP', "Intron", "ARNm", '23S', '16S']
#     for typ in types :
#         #os.makedirs("Nouvelles_donnees/Groupes_homologues/%s"%type, exist_ok=True)
#         verif_homologues_version_2(typ,"Nouvelles_donnees/Groupes_homologues/%s"%typ)
        #regrouper_par_homologues_critere_sequence(type)
       # groupes_ARN_creation_homologues_et_identiques(type)
        
#     liste_deja_fait = []
#     compteur = 0
#     for fic in os.listdir("Nouvelles_donnees/Resultats") :
#         if "test" not in fic and "pickle" in fic and ("dico_graphe" in fic or "dico_sim" in fic) :
#             with open("Nouvelles_donnees/Resultats/"+fic, 'rb') as fichier :
#                 mon_depickler = pickle.Unpickler(fichier)
#                 dico = mon_depickler.load()
#                 print(compteur)
#                 print(fic)
#                 for cle in dico.keys() :
#                     if cle not in liste_deja_fait :
#                         liste_deja_fait.append(cle)
#                 compteur += 1
#                         
#     print(len(liste_deja_fait))
#     with open("Nouvelles_donnes/liste_deja_fait.pickle", 'wb') as fichier_ecriture :
#         mon_pickler = pickle.Pickler(fichier_ecriture)
#         mon_pickler.dump(liste_deja_fait)from Bio.Blast import NCBIWWW


    #entre_deux_chaines()
    #nombre_de_groupes_homologues("4.5S")
    #nombre_de_groupes_identiques("4.5S")
    #recup_type_molecule()
    #identiques_extensions()
        
#     liste_types = ['23S', '28S', '16S', '18S', '25S', '4.5S', '4.8S', "ARNt", "Riboswitch", "Ribozyme", "Intron", "SRP", "ARNm"]
#     for elt in liste_types :
        #groupes_ARN(elt)
    
#     with open("groupes_25S_homologues.pickle", 'rb') as fichier_homologues :
#         mon_depickler_1 = pickle.Unpickler(fichier_homologues)
#         groupes_homologues = mon_depickler_1.load()
#         
#         for groupe in groupes_homologues :
#             print(len(groupe))
        
# #         with open("groupes_23S_identiques.pickle", 'rb') as fichier :
# #             mon_depickler_1 = pickle.Unpickler(fichier)
# #             groupes_identiques = mon_depickler_1.load()
#          
#         compteur = 0
#         compter = 0
#         compter_idem = 0
#         entier = []
#         plusieurs_fois = []
#  
#         liste_plus_grande_diff = []
#         compteur = 0
#         for groupes in groupes_homologues :
#             print(len(groupes))
#             if compteur == 6 :
#                 print(groupes)
#             compteur += 1
#             for groupe in groupes :
#         #for groupe in groupes_homologues :
#             #print(groupe)
#                 plus_grande_difference = 0
#                      
#                 for i in range(len(groupe)) : 
#                     with open(NEW_EXTENSION_PATH_TAILLE + "fichier_"+str(groupe[i][0])+"_"+str(groupe[i][1])+".pickle", 'rb') as fichier_1 :
#                         mon_depickler_1 = pickle.Unpickler(fichier_1)
#                         extension1 = mon_depickler_1.load()
#                         for j in range(i+1, len(groupe)) :
#                             with open(NEW_EXTENSION_PATH_TAILLE + "fichier_"+str(groupe[j][0])+"_"+str(groupe[j][1])+".pickle", 'rb') as fichier_2 :
#                                 mon_depickler_2 = pickle.Unpickler(fichier_2)
#                                 extension2 = mon_depickler_2.load()
#                                 if not homologue(extension1, extension2) :
#     #                                 print(extension1.nodes[1]["position"])
#     #                                 print(extension2.nodes[1]["position"])
#     #                                 print("ramous")    
#                                         if abs(extension1.nodes[1]["position"][0]-extension2.nodes[1]["position"][0]) > plus_grande_difference :
#                                             plus_grande_difference = abs(extension1.nodes[1]["position"][0]-extension2.nodes[1]["position"][0])
#                                          
#                 print("plus grande diff")
#                 print(plus_grande_difference)
#                 liste_plus_grande_diff.append(plus_grande_difference)
#              
# #                 for elt1 in groupe :
# #                     if elt in entier :
# #     #                     print("plusieurs_fois")
# #                         if elt not in plusieurs_fois :
# #                             plusieurs_fois.append(elt)
# #     #                         print(elt)
# #     #                         print(plusieurs_fois)
# #                     else :
# #                         entier.append(elt)
# #                         
# #                     if elt == ('4ujc', 2) or elt == ('4d67', 3):
# #                         print(groupe)
#             compter += len(groupe)
#  
#             compteur += 1
#             print(len(groupe))
#         print(len(groupes_homologues))        
#         print(compter_idem)    
#         print(compter)
#         print(compteur)
#         print(len(entier))
#         print(len(plusieurs_fois))
#         print(liste_plus_grande_diff)
#         
#         compter = 0
#         for groupe in groupes_identiques :
#             #for elt in groupe :
#             compter += len(groupe)
#                 
#         print(compter)
        #print(plusieurs_fois)
        
        #print(groupes_homologues)
    
    
    
#     with open("fichiers_idem_taille_4.pickle", 'rb') as fichier :
#         mon_depickler = pickle.Unpickler(fichier)
#         groupes = mon_depickler.load() 
#         
#         print(groupes)
        