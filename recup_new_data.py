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
    
    
    m = p.match(graphe1.nodes[(ext1.nodes[1]["num_ch"],ext1.nodes[1]["position"][0])]["fr3d"])
    if m :
        num = m.group()
    else :
        print("probleme5")
    pos11 = int(num)
    
    m = p.match(graphe2.nodes[(ext2.nodes[1]["num_ch"],ext2.nodes[1]["position"][0])]["fr3d"])
    if m :
        num = m.group()
    else :
        print("probleme5")
    pos21 = int(num)
    
    m = p.match(graphe1.nodes[(ext1.nodes[2]["num_ch"],ext1.nodes[2]["position"][0])]["fr3d"])
    if m :
        num = m.group()
    else :
        print("probleme5")
    pos12 = int(num)
    
    m = p.match(graphe2.nodes[(ext2.nodes[2]["num_ch"],ext2.nodes[2]["position"][0])]["fr3d"])
    if m :
        num = m.group()
    else :
        print("probleme5")
    pos22 = int(num)

    m = p.match(graphe1.nodes[(ext1.nodes[5]["num_ch"],ext1.nodes[5]["position"][0])]["fr3d"])
    if m :
        num = m.group()
    else :
        print("probleme5")
    pos15 = int(num)
    
    m = p.match(graphe2.nodes[(ext2.nodes[5]["num_ch"],ext2.nodes[5]["position"][0])]["fr3d"])
    if m :
        num = m.group()
    else :
        print("probleme5")
    pos25 = int(num)

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
    print("composantes connexes")
    print(composantes_connexes)
    print(len(composantes_connexes))
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
- distribution des valeurs de RMSD en fonction des valeurs de score
- creation de fichiers csv pour Gephi stockant le graphe complet des occurrences de ce type d'ARN avec les valeurs de score d'alignement et de RMSD en ponderation des aretes
- creation de fichiers csv pour Gephi stockant le graphe des occurrences de ce type d'ARN avec uniquement les aretes au-dessus du seuil de 70 pour le score et de 2.5 pour la RMSD
(- affichage du graphe avec les aretes au-dessus d'un seuil qu'on definit
- clustering DBSCAN sur le graphe complet)
num_ARN : type d'ARN
seuil_hom : seuil de nombres de nts pour les positions similaires
dico_sim, seuil_sim : pas utilise'''
def traitement_aln_global(num_ARN, dico_sim, seuil_hom, seuil_sim):
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
        
    with open("fichier_groupe_10_2_des_bizarres.csv", 'w', newline='') as fichier_csv :
        csvwriter = csv.writer(fichier_csv) 
        csvwriter.writerow(["Num paires", "Score min", "RMSD", "Diff des positions"])
        
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
        dico_min_pos_scores = {}
        
        distrib_min_cluster = []
        distrib_min_cluster_rmsd = []

        #clusters = traitement_res_cdhit(num_ARN)
        for i in range(len(liste_representant)) :
                for j in range(i+1, len(liste_representant)) :
                    scores = []
                    for k in range(1,4) :
                        scores.append(recup_score("Nouvelles_donnees/alignements_%s/needle_%s_%s_%s_%s_seq%d.txt"%(num_ARN,liste_representant[i][0],liste_representant[i][1], liste_representant[j][0],liste_representant[j][1], k)))
                    #print(scores)
                    distrib_min_tot.append(min(scores))
                    dico_min_tot.update({(liste_representant[i], liste_representant[j]) : min(scores)})
                    
                    
                    
                    nom1 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[i][0], liste_representant[i][1])
                    nom2 = "fichier_%s_%d_taille_4.pdb"%(liste_representant[j][0], liste_representant[j][1])
                    
                    if (nom1, nom2) in rmsd.keys() :
#                         if (nom1 == "fichier_4u27_10_taille_4.pdb" and nom2 == "fichier_4y4p_20_taille_4.pdb") or (nom2 == "fichier_4y4p_20_taille_4.pdb" and nom1 == "fichier_4u27_10_taille_4.pdb") :
#                             print("celuila")
#                             print(rmsd[(nom1, nom2)])
#                             print(scores)
#                             exit(0)
#                         if min(scores) > 150 and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)]  > 12   :
#                                             print("petit rat mou")
#                                             print(min(scores))
#                                             print(rmsd[(nom1, nom2)])
#                                             print((liste_representant[i], liste_representant[j]))
#                                             csvwriter.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)]])
#                             
#                         if min(scores) < 40 :
#                             print("petit rat")
#                             print(min(scores))
#                             print(rmsd[(nom1, nom2)])
#                             print((liste_representant[i], liste_representant[j]))
                        distrib_min_rmsd.append(rmsd[(nom1, nom2)])
                        dico_min_tot_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom1, nom2)]})
                        
#                         
#                         for cluster in clusters :
#                             if liste_representant[i] in cluster and liste_representant[j] in cluster :
#                                 distrib_min_cluster.append(min(scores))
#                                 distrib_min_cluster_rmsd.append(rmsd[(nom1, nom2)])
                        
                    else :
#                         if min(scores) > 150 and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)]  > 12   :
#                                             print("petit rat mou")
#                                             print(min(scores))
#                                             print(rmsd[(nom2, nom1)])
#                                             print((liste_representant[i], liste_representant[j]))
#                                             csvwriter.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom2, nom1)]])

#                         if (nom1 == "fichier_1vq8_19_taille_4.pdb" and nom2 == "fichier_1yhq_24_taille_4.pdb") or (nom2 == "fichier_1vq8_19_taille_4.pdb" and nom1 == "fichier_1yhq_24_taille_4.pdb") :
#                             print("celuila")
#                             print(rmsd[(nom2, nom1)])
#                             print(scores)
                            #exit(0)
#                         if min(scores) < 80 and min(scores) > 70 and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)] < 3 :
#                             print("petit rat")
#                             print(min(scores))
#                             print(rmsd[(nom2, nom1)])
#                             print((liste_representant[i], liste_representant[j]))
                        distrib_min_rmsd.append(rmsd[(nom2, nom1)])
                        dico_min_tot_rmsd.update({(liste_representant[i], liste_representant[j]) : rmsd[(nom2, nom1)]})
                        
#                         for cluster in clusters :
#                             if liste_representant[i] in cluster and liste_representant[j] in cluster :
#                                 distrib_min_cluster.append(min(scores))
#                                 distrib_min_cluster_rmsd.append(rmsd[(nom2, nom1)])
#                         print(liste_representant[i], liste_representant[j])
                    
        
#                         toutes_elevees = True
#                         for elt in scores :
#                             if elt < 100 :
#                                 toutes_elevees = False
                     
#                         if toutes_elevees :
#                             print(liste_representant[i])
#                             print(liste_representant[j])
#                             if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
#                                 print(dico_sim[(liste_representant[i], liste_representant[j])]["sim"])
#                             else :
#                                 print(dico_sim[(liste_representant[j], liste_representant[i])]["sim"])
                    
                    with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[i][0], liste_representant[i][1]), "rb") as fichier_1 :
                        mon_depickler_1 = pickle.Unpickler(fichier_1)
                        graphe1 = mon_depickler_1.load()
                         
                    with open("Nouvelles_donnees/fichier_%s_%s_2.pickle"%(liste_representant[j][0], liste_representant[j][1]), "rb") as fichier_2 :
                        mon_depickler_2 = pickle.Unpickler(fichier_2)
                        graphe2 = mon_depickler_2.load()
                        
                        #pos_similaire = pos_similaire_fr3d(graphe1, graphe2, liste_representant[i], liste_representant[j], seuil_hom)
                        pos_similaire = positions_similaires(graphe1, graphe2, seuil_hom)
                        if pos_similaire[0] :
                            distrib_diff_pos.append(abs(abs(graphe1.nodes[1]["position"][0]-graphe1.nodes[2]["position"][0])-abs(graphe2.nodes[1]["position"][0]-graphe2.nodes[2]["position"][0])))
                            if liste_representant[i][0] != liste_representant[j][0] :
#                                 if min(scores) == 39.5 :
#                                     print(graphe1.nodes[1]["position"])
#                                     print(graphe2.nodes[1]["position"])
#                                     print(graphe1.nodes[2]["position"])
#                                     print(graphe2.nodes[2]["position"])
#                                     print(liste_representant[i], liste_representant[j])
                                if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
                                    #if dico_sim[(liste_representant[i], liste_representant[j])]["sim"] >= seuil_sim :
                                        distrib_min_pos_sim.append(min(scores))
                                        dico_min_pos_scores.update({(liste_representant[i], liste_representant[j]) : (min(scores))})
                                        dico_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (min(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                else :
                                    #if dico_sim[(liste_representant[j], liste_representant[i])]["sim"] >= seuil_sim :
                                        distrib_min_pos_sim.append(min(scores))
                                        dico_min_pos_scores.update({(liste_representant[i], liste_representant[j]) : (min(scores))})
                                        dico_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (min(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
#                                 print("\n")
#                                 print(len(dico_min_pos_rmsd))
#                                 print(len(dico_min_pos_scores))
                                if (nom1, nom2) in rmsd.keys() :
#                                     if min(scores) > 85 and min(scores) < 120 and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)]  > 3  and rmsd[(nom1, nom2)]  < 5  :
#                                         print("petit rat mou")
#                                         print(min(scores))
#                                         print(rmsd[(nom1, nom2)])
#                                         print(pos_similaire)
#                                         print((liste_representant[i], liste_representant[j]))
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)], pos_similaire])
                                    entree = False
                                    ajoute = False
                                    not_ajoute = False
                                    liste_a_enlever = []
                                    for cle in dico_min_pos_rmsd.keys() :
                                        if not (liste_representant[i] == cle[0] and liste_representant[j] == cle[1]) and not (liste_representant[i] == cle[1] and  liste_representant[j] == cle[0]) and ((liste_representant[i] in cle and (liste_representant[j][0] == cle[0][0] or liste_representant[j][0] == cle[1][0])) or (liste_representant[j] in cle and (liste_representant[i][0] == cle[0][0] or liste_representant[i][0] == cle[1][0])) ):
                                            entree = True
#                                             print("rapala")
#                                             print(cle)
#                                             print(dico_min_pos_rmsd[cle])
#                                             print(dico_min_pos_sim[cle])
#                                             print(nom1, nom2)
#                                             print(rmsd[(nom1, nom2)])
#                                             print(dico_min_pos_sim[(liste_representant[i], liste_representant[j])])
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
#                                     if min(scores) > 85 and min(scores) < 120 and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)]  > 3  and rmsd[(nom2, nom1)]  < 5  :
#                                         print("petit rat mou")
#                                         print(min(scores))
#                                         print(rmsd[(nom1, nom2)])
#                                         print((liste_representant[i], liste_representant[j]))
#                                         csvwriter.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)], pos_similaire])
                                    entree = False
                                    ajoute = False
                                    not_ajoute = False
                                    liste_a_enlever = [] #                                             print("rapala")
#                                             print(cle)
#                                             print(dico_min_pos_rmsd[cle])
#                                             print(dico_min_pos_sim[cle])
#                                             print(nom1, nom2)
#                                             print(rmsd[(nom1, nom2)])
#                                             print(dico_min_pos_sim[(liste_representant[i], liste_representant[j])])
                                    for cle in dico_min_pos_rmsd.keys() :
                                        if not (liste_representant[i] == cle[0] and liste_representant[j] == cle[1]) and not (liste_representant[i] == cle[1] and  liste_representant[j] == cle[0]) and ((liste_representant[i], liste_representant[j]) != cle and (liste_representant[j], liste_representant[i]) != cle and (liste_representant[i] in cle and (liste_representant[j][0] == cle[0][0] or liste_representant[j][0] == cle[1][0])) or (liste_representant[j] in cle and (liste_representant[i][0] == cle[0][0] or liste_representant[i][0] == cle[1][0]))) :
                                            entree = True
#                                             print("rapala")
#                                             print(cle)
#                                             print(dico_min_pos_rmsd[cle])
#                                             print(dico_min_pos_sim[cle])
#                                             print(nom1, nom2)
#                                             print(rmsd[(nom2, nom1)])
#                                             print(dico_min_pos_sim[(liste_representant[i], liste_representant[j])])
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
                                
#                             if min(scores) <= 100 :
#                                 print("rat")
#                                 print((liste_representant[i], liste_representant[j]))
#                                 print(min(scores))
#                                 print(scores)
#                                 liste_seq_min[scores.index(min(scores))] += 1
#                                 liste_en_dessous_de_100.append((liste_representant[i], liste_representant[j], scores))
#                                 if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
#                                     liste_en_dessous_de_100_sim.append(dico_sim[(liste_representant[i], liste_representant[j])]["sim"])
#                                 else :
#                                     liste_en_dessous_de_100_sim.append(dico_sim[(liste_representant[j], liste_representant[i])]["sim"])
#                                 liste_en_dessous_de_100_pos.append(pos_similaire)
#                                 
#                                 if (nom1, nom2) in rmsd.keys() :
#                                     liste_en_dessous_de_100_rmsd.append(rmsd[(nom1, nom2)])
#                                 else :
#                                     liste_en_dessous_de_100_rmsd.append(rmsd[(nom2, nom1)])
                        else :
                                dico_not_min_pos_sim.update({(liste_representant[i], liste_representant[j]) : (min(scores), pos_similaire[1], pos_similaire[2], pos_similaire[3]) })
                                
                                if (nom1, nom2) in rmsd.keys() :
                                    if min(scores) > 150 and min(scores) < 200 and rmsd[(nom1, nom2)] != None and rmsd[(nom1, nom2)]  < 2.5   :
                                                        print("petit rat mou")
                                                        print(min(scores))
                                                        print(rmsd[(nom1, nom2)])
                                                        print((liste_representant[i], liste_representant[j]))
                                                        csvwriter.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom1, nom2)], pos_similaire])
                         
                                else :
                                    if min(scores) > 150 and min(scores) < 200 and rmsd[(nom2, nom1)] != None and rmsd[(nom2, nom1)]  < 2.5   :
                                            print("petit rat mou")
                                            print(min(scores))
                                            print(rmsd[(nom2, nom1)])
                                            print((liste_representant[i], liste_representant[j]))
                                            csvwriter.writerow([(liste_representant[i], liste_representant[j]), min(scores), rmsd[(nom2, nom1)], pos_similaire])
                     
                    
                           
                    if (liste_representant[i], liste_representant[j]) in dico_sim.keys() :
                        if dico_sim[(liste_representant[i], liste_representant[j])]["sim"] == 1.0 :
                            if min(scores) == 72 :
                                print("gros rat")
                                print(min(scores))
                                print(liste_representant[i], liste_representant[j])
#                             print(liste_representant[i])
#                             print(liste_representant[j])
#                             print(dico_sim[(liste_representant[i], liste_representant[j])]["sim"])
#                             print(scores)
                            distrib_min.append(min(scores))
                            dico_min.update({(liste_representant[i], liste_representant[j]) : min(scores) })
                    else :
                        if dico_sim[(liste_representant[j], liste_representant[i])]["sim"] == 1.0 :
                            if min(scores) == 72 :
                                print("gros rat")
                                print(min(scores))
                                print(liste_representant[i], liste_representant[j])
#                             print(liste_representant[i])
#                             print(liste_representant[j])
#                             print(dico_sim[(liste_representant[j], liste_representant[i])]["sim"])
#                             print(scores)
                            distrib_min.append(min(scores))
                            dico_min.update({(liste_representant[j], liste_representant[i]) : min(scores) })
        
        
#         print(distrib_min_tot)
#         print(distrib_min_pos_sim)
#         print(min(distrib_min_pos_sim))
#         print(distrib_diff_pos)

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
        plt.xlim(min_tot, 351)

        axs.set_xlabel("Score minimum d'alignement global")
        axs.set_ylabel("Nombre de paires")
        axs.set_title("Distribution des valeurs de scores d'alignement global \n (minimum des 3 valeurs) pour les %s \n (avec seuil positions similaires : %d min scores : %d)"%(num_ARN, seuil_hom, min_pos_sim))
        #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_seuil_pos_sim_clustering_%d.png"%(num_ARN, seuil_hom))
        
        calcul_kmeans_sur_distrib(dico_min_tot, dico_min_tot_rmsd)
        
        plt.show()
        plt.close()
        
        distrib_min_tot.remove(-20.0)
        distrib_min_pos_sim.remove(-20.0)
        distrib_min_cluster.remove(-20.0)
        print(len(dico_min_pos_scores))
        print(dico_min_pos_scores.values())
        print(len(dico_min_pos_rmsd))
        print(dico_min_pos_rmsd.values())
        plt.figure(figsize=(12,6))
        axs = plt.gca()
        plt.scatter(dico_min_tot.values(), dico_min_tot_rmsd.values())
        plt.scatter(dico_min_pos_scores.values(), dico_min_pos_rmsd.values())
        plt.scatter(distrib_min_cluster, distrib_min_cluster_rmsd)
        axs.set_xlabel("Score minimum d'alignement global")
        axs.set_ylabel("RMSD (en Angstrom)")
        axs.set_xticks(np.arange(0, 310, 10))
        axs.set_title("Distribution des valeurs de RMSD en fonction des valeurs de score d'alignement global \n (seuil positions similaires : %d)"%seuil_hom)
        #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_seuil_pos_rmsd_clustering_%d.png"%(num_ARN, seuil_hom))
        plt.show()
        plt.close()
        
        #return liste_en_dessous_de_100, liste_en_dessous_de_100_sim, liste_en_dessous_de_100_pos, liste_en_dessous_de_100_rmsd
        
        
#         if len(distrib_diff_pos) > 0 :
#             plt.figure(figsize=(9,6))
#             axs = sns.distplot(distrib_diff_pos, kde=False)
#             axs.set_xlabel("Distance entre les brins du motif")
#             axs.set_ylabel("Nombre de paires")
#             axs.set_title("Distribution des valeurs de distance entre brins du motif \n pour les paires dont les positions sont proches (seuil : %d)"%(seuil_hom))
#             #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_pos_diff_seuil_pos_sim_%d.png"%(num_ARN, seuil_hom))
    
            #plt.show()
            #plt.close()
        
        
        
        
        if len(distrib_min) > 0 :
            zoom = min(distrib_min)
        else :
            zoom = 0
        for elt in distrib_min_tot :
            if elt >= zoom : 
                distrib_min_tot_zoom.append(elt)
        
        plt.figure(figsize=(8,8))
        axs = sns.distplot(distrib_min_tot_zoom, kde=False, bins = len(liste_representant))
        sns.distplot(distrib_min_pos_sim, kde=False, bins = len(liste_representant))
        axs.set_xticks(np.arange(1,352,50))
        axs.set_xlabel("Score minimum d'alignement global")
        axs.set_ylabel("Nombre de paires")
        axs.set_title("Distribution des valeurs de scores d'alignement global \n (minimum des 3 valeurs) pour les %s \n (zoom sur valeurs élevées)"%num_ARN)
        #plt.savefig("Nouvelles_donnees/alignements_%s/distrib_score_zoom_%d.png"%(num_ARN, zoom))
        
        plt.show()
        plt.close()
        
        graph = nx.Graph()
        for cle in dico_min_tot.keys() :
            if cle[0] not in graph.nodes() :
                graph.add_node(cle[0])
            if cle[1] not in graph.nodes() :
                graph.add_node(cle[1])
            #if dico_min_tot[cle] >= 71 :
            nom1 = "fichier_%s_%d_taille_4.pdb"%(cle[0][0], cle[0][1])
            nom2 = "fichier_%s_%d_taille_4.pdb"%(cle[1][0], cle[1][1])
            if (nom1, nom2) in rmsd.keys() :
                graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=rmsd[(nom1, nom2)])
            else :
                graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], rmsd=rmsd[(nom2, nom1)])
#         for cle in dico_not_min_pos_sim.keys() :
#             if (cle[0], cle[1]) not in graph.edges() :
#                 graph.add_edge(cle[0], cle[1], aln=dico_not_min_pos_sim[cle][0], val_pos1=dico_not_min_pos_sim[cle][1], val_pos2 = dico_not_min_pos_sim[cle][2])

        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_no_seuil.csv"%(num_ARN,num_ARN),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["source", "target","weight"])
        
            for u,v,data in graph.edges(data=True) :
                csvwriter.writerow([u,v, data["aln"]])
        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_no_seuil_noeuds.csv"%(num_ARN, num_ARN),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["id", "clustering_cdhit"])
            
            for noeud,data in graph.nodes(data=True) :
                num_cluster = -1
                compteur = 1
#                 for cluster in clusters :
#                     if noeud in cluster :
#                         if num_cluster != -1 : print("bizarre")
#                         num_cluster = compteur
#                     compteur += 1  
                csvwriter.writerow([noeud, num_cluster])
        
        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_seuil_70_2.5.csv"%(num_ARN,num_ARN),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["source", "target", "aln", "rmsd"])
            
            a_enlever = []
            for u,v,data in graph.edges(data=True) :
                if data["aln"] >= 70 and data["rmsd"] != None and data["rmsd"] <= 2.5 :
                    csvwriter.writerow([u,v, data["aln"], round(data["rmsd"],2)])
                else :
                    a_enlever.append((u,v))
        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_pos_sim_%s_test_seuil_70_2.5_noeuds.csv"%(num_ARN, num_ARN),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["id", "clustering_cdhit", "num_ARN"])
            
            for noeud,data in graph.nodes(data=True) :
                num_type = ""
                for cle in liste_cles :
                    if noeud in liste_cles[cle] :
                        num_type = cle
                num_cluster = -1
                compteur = 1
#                 for cluster in clusters :
#                     if noeud in cluster :
#                         if num_cluster != -1 : print("bizarre")
#                         num_cluster = compteur
#                     compteur += 1  
                csvwriter.writerow([noeud, num_cluster, num_type])
        
        ### faire les groupes d'homologues par composantes connexes au-dessus des seuils : score 70, rmsd 2
        graph_seuils = graph.copy()
        for elt in a_enlever :
            graph_seuils.remove_edge(elt[0], elt[1])
            
        composantes = recherche_composante_connexe(graph_seuils)
        print(composantes)
#         if "groupes_homologues_version_4janv.pickle" in os.listdir("/home/coline/eclipse-workspace/a_minor_extension_V2/recup_data/") :
#             with open("groupes_homologues_version_4janv.pickle", 'rb') as fichier_hom :
#                 mon_depickler = pickle.Unpickler(fichier_hom)
#                 homologues = mon_depickler.load()
#              
#                 homologues.extend(composantes)
#             with open("groupes_homologues_version_4janv.pickle", 'wb') as fichier_hom_new :
#                 mon_pickler = pickle.Pickler(fichier_hom_new)
#                 mon_pickler.dump(homologues)
#         else :
#             with open("groupes_homologues_version_4janv.pickle", 'wb') as fichier_hom_new :
#                 mon_pickler = pickle.Pickler(fichier_hom_new)
#                 mon_pickler.dump(composantes)
            
            
        return liste_en_dessous_de_100, liste_en_dessous_de_100_sim, liste_en_dessous_de_100_pos, liste_en_dessous_de_100_rmsd
    
        dico_trie = sorted(dico_min_tot.items(), key=lambda t: t[1], reverse=True)
        graph = nx.Graph()
        
        for cle in dico_min_tot.keys() :
            if cle[0] not in graph.nodes() :
                graph.add_node(cle[0])
            if cle[1] not in graph.nodes() :
                graph.add_node(cle[1])
            if cle in dico_sim :
                #if dico_min_tot[cle] >= 72 :   
                graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], sim=dico_sim[cle]["sim"])
            else :
                graph.add_edge(cle[0], cle[1], aln=dico_min_tot[cle], sim=dico_sim[(cle[1], cle[0])]["sim"])
        
        graphe_renomme = nx.convert_node_labels_to_integers(graph, first_label=0, label_attribute="nom")
        
        matrice_distance = matrice_distance_selon_score_aln(graphe_renomme)
        clustering = DBSCAN(min_samples=2, eps=225.0, metric='precomputed').fit(matrice_distance)
    #         tab_clustering = [[]]
    #         print(len(tab_clustering))
    #         print(len(clustering.labels_))
    #         print(clustering.labels_)
    #         compteur = 0
    #         for elt in clustering.labels_ :
    #             if elt != -1 :
    #                 for _ in range(len(tab_clustering), elt + 1) :
    #                         tab_clustering.append([])
    #                 tab_clustering[elt].append(compteur)
    #             compteur += 1
        
    #         graph_epure = nx.Graph()
    #         for cle in dico_min_tot.keys() :
    #             if cle[0] not in graph_epure.nodes() :
    #                 graph_epure.add_node(cle[0])
    #             if cle[1] not in graph_epure.nodes() :
    #                 graph_epure.add_node(cle[1])
    #             
    #             if dico_min_tot[cle] >= 72 :   
    #                 graph_epure.add_edge(cle[0], cle[1], aln=dico_min_tot[cle])
        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_clustering_homologie_%s.csv"%(num_ARN, num_ARN),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["source", "target", "label"])
        
            for u,v,data in graphe_renomme.edges(data=True) :
                if data["aln"] >= zoom :
                    csvwriter.writerow([u,v,(data["aln"], round(data["sim"],2))])
        
        with open("Nouvelles_donnees/alignements_%s/fichier_csv_clustering_homologie_%s_noeuds.csv"%(num_ARN, num_ARN),'w') as fichier_csv:
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["id", "label", "clustering_dbscan"])
            
            for noeud,data in graphe_renomme.nodes(data=True) :
#                 with open("fichier_infos_molecules_new.csv", 'r', newline='') as fichier_csv :
#                     csvreader = csv.reader(fichier_csv) 
#                     type_arn_precis = ""
#                     for ligne in csvreader :
#                         if len(ligne) == 7 :
#                             #print(ligne[1])
#                             print(data["nom"][0])
#                             if ligne[1] == data["nom"][0].upper() :
#                                 type_arn_precis = ligne[3]
#                                 print("petit rat")
                        
                csvwriter.writerow([noeud, data["nom"], clustering.labels_[noeud]])
    
        print(len(liste_representant))
            
            
    #         print(graph.number_of_nodes())
    #         print(graph.number_of_edges())
    #         
    #         print(len(liste_representant))
    #         print(len(dico_min_tot))
    #         
    #         pos = spring_layout(graph)
    #         nx.draw_networkx_nodes(graph, pos, node_size = 5)
    #         nx.draw_networkx_edges(graph, pos)
    #         plt.show()
    #         print(min(distrib_min))
    #         
    #         composantes = recherche_composante_connexe(graph)
    #         for c in composantes :
    #             cliques = nx.find_cliques(graph.subgraph(c))
    #             for clique in cliques : 
    #                 if len(clique) == len(c) :
    #                     if len(c) < 5 :
    #                         print(c)
    #                         for i in range(len(c)) :
    #                             for j in range(i+1, len(c)) :
    #                                 with open("Nouvelles_donnees/fichier_%s_%s.pickle"%(c[i][0], c[i][1]), "rb") as fichier_1 :
    #                                     mon_depickler_1 = pickle.Unpickler(fichier_1)
    #                                     graphe1 = mon_depickler_1.load()
    #                                     
    #                                 with open("Nouvelles_donnees/fichier_%s_%s.pickle"%(c[j][0], c[j][1]), "rb") as fichier_2 :
    #                                     mon_depickler_2 = pickle.Unpickler(fichier_2)
    #                                     graphe2 = mon_depickler_2.load()
    #                                 
    #                                 if homologue(graphe1, graphe2) :
    #                                     print(c[i])
    #                                     print(c[j])
    #                                     print("ok")
    #                 subgraph = nx.subgraph(graph, c)
    #                 if subgraph.number_of_edges() != (len(c)*(len(c)-1))/2 :
    #                     for i in range(len(c)) :
    #                         for j in range(i+1, len(c)) :
    #                             if (c[i], c[j]) not in graph.edges() :
    #                                 if (c[i], c[j]) in dico_min_tot.keys() :
    #                                     print(dico_min_tot[(c[i], c[j])])
    #                                 else :
    #                                     print(dico_min_tot[(c[j], c[i])])
                                    
            
        for elt in dico_trie :
                if elt[0][0] not in graph.nodes() :
                    graph.add_node(elt[0][0])
                if elt[0][1] not in graph.nodes() :
                    graph.add_node(elt[0][1])
                     
                graph.add_edge(elt[0][0], elt[0][1], aln = elt[1])
                cliques = nx.find_cliques(graph)
                 
                for clique in cliques :
                    c = clique
                    break
                if len(c) != graph.number_of_nodes() :
                    break
        print(graph.nodes.data())
        print(graph.edges.data())


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
              
if __name__ == '__main__':
    #creation_fichier_input_cdhit("23S")
    #traitement_res_cdhit("23S")
    
    types_arn = ["18S", "16S","Ribozyme", "Riboswitch", "28S", "25S", "Intron", "arnt_16s_arnm", "arnt_16s"]
      
    with open("/media/coline/Maxtor/dico_new_avec_derniere_modif_encore_modif_612.pickle", 'rb') as fichier_dico_sim :
        mon_depickler_2 = pickle.Unpickler(fichier_dico_sim)
        dico_sim = mon_depickler_2.load()
#     
#     for typ in types_arn :
#         alignement_global(["23S","25S","28S"])
        #for typ in types_arn :
        with open("Nouvelles_donnees/alignements_['23S', '25S', '28S']/fichier_nb_pos_sim_en_dessous_de_100.csv", 'w', newline='') as fichier_csv_ecriture :
            csvwriter = csv.writer(fichier_csv_ecriture)
            csvwriter.writerow(["seuil", "nb en dessous de 100"])
            for seuil_pos in np.arange(30,40,10) :
                #for seuil_sim in np.arange(0,1.1,0.1) :
                liste_en_dessous_de_100, liste_en_dessous_de_100_sim, liste_en_dessous_de_100_pos, liste_en_dessous_de_100_rmsd = traitement_aln_global(['23S', '25S', '28S'], dico_sim, seuil_pos, 0.7)
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
        