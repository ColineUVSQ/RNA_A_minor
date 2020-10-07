'''
Created on 25 juil. 2020

@author: coline
'''
import pickle
import csv
import os
import urllib
from gemmi import cif
import networkx as nx
from posix import setreuid
from recup_data.extension_new import stockage_motif
from recup_data.recup_new_data import recup_sequences_autour_motif,\
    get_resolution, recup_score
from recup_data.constantes import PATH_MMCIF
from Bio.Emboss.Applications import NeedleCommandline
import random
from random import randint
from recup_data.traitement_clusters_graphe_commun import recherche_graphe_commun,\
    etablissement_expr_reg
import re
from networkx.algorithms import isomorphism
from recup_data.new_algo_comparaison import recup_chaines
import numpy as np
from urllib import request


def stocker_type_arn():
    '''
        stocke les types d'arn de toutes les structures de la PDB (stockées dans le fichier graphs_all_2020_07_with_SSEs.pickle)
        dans un dictionnaire
    '''
    with open("Resultats_sequences/graphs_all_2020_07_with_SSEs.pickle", 'rb') as fichier_pickle_1 :
        mon_depickler = pickle.Unpickler(fichier_pickle_1)
        graphes = mon_depickler.load()
        
        with open("type_arn_tot.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            type_arn = mon_depickler.load()
        

        compteur = 0
        for cle in graphes :
            if type_arn[cle] == -1 : 
            #type_ARN = -1
            
                if cle[0].upper()+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
                            print("petit rat")
                            url = 'http://www.rcsb.org/pdb/files/%s.cif' % cle[0].upper()
                            urllib.request.urlretrieve(url, PATH_MMCIF+cle[0].upper()+".cif")
                try:
                    doc = cif.read_file(PATH_MMCIF+cle[0].upper()+".cif")
                    block = doc.sole_block()
                    cat = block.find_mmcif_category("_entity_poly.") 
                    print(list(cat.tags))
                    id = -1
                    for row in cat :
                                        #print(tuple(row))
                        num_ch = row[6].split(",")
                                        
            #                             print("petit rat")
            #                             print(num_ch)
                        if cle[1] in num_ch :
                            id = row[0]
                    
                    cat = block.find_mmcif_category("_entity.") 
                                      
                    for row in cat :
                                        #print(tuple(row))
                        if row[0] == id :
                            type_ARN = row[3]
                except :
                    print("probleme")
                
                type_arn.update({cle : type_ARN})
            compteur += 1
            print(compteur)
            
        print(type_arn)
        
        with open("type_arn_tot.pickle", 'wb') as fichier_pickle :
            mon_pickler = pickle.Pickler(fichier_pickle)
            mon_pickler.dump(type_arn)



def stocker_resolution_tot():
    ''' pour chaque structure PDB 
        recupere la resolution dans le fichier PDB
        stocke toutes les resolutions dans un dictionnaire 
    '''
    with open("graphs_all_2020_07_with_SSEs.pickle", 'rb') as fichier_pickle_1 :
        mon_depickler = pickle.Unpickler(fichier_pickle_1)
        graphes = mon_depickler.load()
    
        with open("resolutions_tot.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
            
            compteur = 0
            for cle in graphes :
                if cle[0] not in resolutions.keys() :
                    print(cle[0])
                    print(compteur)
                    res = get_resolution(cle[0])
                    if res == None :
                        if cle[0].upper()+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
                            print("petit rat")
                            url = 'http://www.rcsb.org/pdb/files/%s.cif' % cle[0].upper()
                            urllib.request.urlretrieve(url, PATH_MMCIF+cle[0].upper()+".cif")
                        
                        try:
                            doc = cif.read_file(PATH_MMCIF+cle[0].upper()+".cif")
                            
                            block = doc.sole_block()
                            print(block.name)
                
                             
                            cat = block.find("_em_3d_reconstruction.",["resolution"])
                            print(cat) 
                            for row in cat :
                                print(row[0])
                                res = float(row[0])
                        except:
                            print("Error in RNA %s" % PATH_MMCIF+cle[0].upper()+".cif")
                            #exit()
                    
                    resolutions.update({cle[0] : res})
                    compteur += 1
            
    with open("resolutions_tot.pickle", 'wb') as fichier_pickle :
            mon_pickler = pickle.Pickler(fichier_pickle)
            mon_pickler.dump(resolutions)
            

def non_redondantes_seq_RNA3Dmotifatlas_new(liste_tout_aminor):
    '''
    Recherche la redondance dans les structures d'ARN de la PDB, à l'aide des classes de représentants de RNA 3D motif Atlas
    '''
    with open("resolutions_tot.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        resolutions = mon_depickler.load()
        
    with open("Resultats_sequences/dico_organismes.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        orgs = mon_depickler.load()
    
    with open("Resultats_sequences/nrlist_3.136_all.csv", 'r') as fichier_csv_motifs_atlas:
        csvreader = csv.reader(fichier_csv_motifs_atlas)
        
        
        
        liste_plus = []
        dico_redondance = {}
        for row in csvreader :
            if '+' not in row[1] : ## cas de base
                
                repr = (row[1].split("|")[0], row[1].split("|")[2])
                
                representes = []
                autres = row[2].split(",")
                for elt in autres :
                    representes.append((elt.split("|")[0], elt.split("|")[2]))
                dico_redondance.update({repr:representes})
            else : ## cas où le représentant de la classe est composé de plusieurs structures (je sais pas pourquoi ?)
                compteur = 0
                repr = []
                for morceau_row in row[1].split('+') :
                    repr.append((morceau_row.split("|")[0], morceau_row.split("|")[2]))
                
                bizarre = False
                representes = [[]]*len(repr)
                print(representes)
                print(len(repr))
                autres = row[2].split(",")
                for elt in autres :
                    c = 0
                    
                    for elt_2 in elt.split('+') :
                        print(repr)
                        if c > len(representes)-1 :
                            representes.append([])
                        representes[c].append((elt_2.split("|")[0], elt_2.split("|")[2]))
                        c += 1
                c = 0
                for elt in repr :
                    dico_redondance.update({elt : representes[c] })
                    
                    c+= 1
                

        
        print(dico_redondance)
        print(len(dico_redondance))
        for elt in dico_redondance :
            print(elt)
            print(dico_redondance[elt])
            
        
            
    
    with open("Resultats_sequences/graphs_all_2020_07_with_SSEs.pickle", 'rb') as fichier_pickle_1 :
        mon_depickler = pickle.Unpickler(fichier_pickle_1)
        graphes = mon_depickler.load()
        print("Nombre d'ARN en tout :")
        print(len(graphes))
        
        
        liste_non_trouvees = []
        liste_pas_repr = []
        liste_nom_graphes = []
        liste_plusieurs = []
        liste_nom_graphes_aminor = []
        
        
        for cle in liste_tout_aminor :
            with open("Nouvelles_donnees/fichier_%s_%s_10.pickle"%(cle[1][0], cle[1][1]), "rb") as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                extension = mon_depickler.load()
                
                if cle[1][0] == '4v67' :
                    print(extension.nodes[1]["num_ch"])
                    print(extension.nodes[3]["num_ch"])
                    print(extension.nodes[5]["num_ch"])
                
                if (cle[1][0].upper(), extension.nodes[1]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[1]["num_ch"]))
                
                if (cle[1][0].upper(), extension.nodes[3]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[3]["num_ch"]))
                
                if (cle[1][0].upper(), extension.nodes[5]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[5]["num_ch"]))
        
        print(len((liste_nom_graphes_aminor)))          
        #exit() 
        #a_faire = [liste_nom_graphes_aminor, graphes]
        liste_a_faire = graphes
        groupes_redondances_aminor = []
        
        
        c = 0
        for elt in liste_nom_graphes_aminor :
            existe = False
            for cle in dico_redondance :
                if cle == elt or elt in dico_redondance[cle] :
                    existe = True
            if not existe :
                c += 1
                print(elt)
        print(c)
        print(len(liste_nom_graphes_aminor))
        #exit()
        
        
#         for cle in dico_redondance :
#             cle_dans = False
#             if cle in liste_nom_graphes_aminor :
#                 groupes_redondances_aminor.append([cle])
#                 cle_dans = True
#             else :
#                 for elt in dico_redondance[cle] :
#                     if elt in liste_nom_graphes_aminor :
#                         if cle_dans :
#                             groupes_redondances_aminor[len(groupes_redondances_aminor)-1].append(elt)
#                            
#                         else :
#                             groupes_redondances_aminor.append([elt])
#                             cle_dans = True
#         
#         print(groupes_redondances_aminor)
#         for elt in groupes_redondances_aminor :
#             print(len(elt))
#             print(elt)
#         print(len(groupes_redondances_aminor))
        #exit()
        liste_nom_graphes = []
        for cle in dico_redondance :
            if cle in liste_a_faire :
                liste_nom_graphes.append([])
                for elt in dico_redondance[cle] :
                    if elt in graphes :
                        liste_nom_graphes[len(liste_nom_graphes)-1].append(elt)
            else :
                a_ajoute = False
                for elt in dico_redondance[cle] :
                    if elt in liste_nom_graphes_aminor:
                        a_ajoute = True
                        break
                if a_ajoute :
                    liste_nom_graphes.append([])
                    for elt in dico_redondance[cle] :
                        if elt in graphes :
                            liste_nom_graphes[len(liste_nom_graphes)-1].append(elt)
                    
        
        print(liste_nom_graphes)
        print(len(liste_nom_graphes))
        
        return liste_nom_graphes




def non_redondantes_seq_RNA3Dmotifatlas(liste_tout_aminor):
    '''
    Recherche la redondance dans les structures d'ARN de la PDB, à l'aide des classes de représentants de RNA 3D motif Atlas
    '''
    with open("resolutions_tot.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        resolutions = mon_depickler.load()
        
    with open("Resultats_sequences/dico_organismes.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        orgs = mon_depickler.load()
    
    with open("Resultats_sequences/nrlist_3.136_all.csv", 'r') as fichier_csv_motifs_atlas:
        csvreader = csv.reader(fichier_csv_motifs_atlas)
        
        
        
        liste_plus = []
        dico_redondance = {}
        for row in csvreader :
            if '+' not in row[1] : ## cas de base
                
                repr = (row[1].split("|")[0], row[1].split("|")[2])
                
                representes = []
                autres = row[2].split(",")
                for elt in autres :
                    representes.append((elt.split("|")[0], elt.split("|")[2]))
                dico_redondance.update({repr:representes})
            else : ## cas où le représentant de la classe est composé de plusieurs structures (je sais pas pourquoi ?)
                compteur = 0
                repr = row[1]
                #for morceau_row in row[1].split('+') :
                    #repr.append((morceau_row.split("|")[0], morceau_row.split("|")[2]))
                
                bizarre = False
                representes = []
                autres = row[2].split(",")
                for elt in autres :
                    for elt_2 in elt.split('+') :
                        representes.append((elt_2.split("|")[0], elt_2.split("|")[2]))
#                         print(elt)
#                         if compteur < len(elt.split("+")) :
#                             representes.append((elt.split('+')[compteur].split("|")[0], elt.split('+')[compteur].split("|")[2]))
#                         else :
#                             bizarre = True
#                             representes.append((elt.split('+')[compteur%len(elt.split('+'))].split("|")[0], elt.split('+')[compteur%len(elt.split('+'))].split("|")[2]))
#                             if repr not in liste_plus :
#                                 liste_plus.append(repr)
#                             for elt in autres : 
#                                 for elt_2 in elt.split('+') :
#                                     if (elt_2.split("|")[0], elt_2.split("|")[2]) not in liste_plus :
#                                         liste_plus.append((elt_2.split("|")[0], elt_2.split("|")[2]))
#                             break
                    #print(repr)
                    #print(representes)
                    #if not bizarre :
                    #    dico_redondance.update({repr:representes})
                    #compteur += 1
                dico_redondance.update({repr:representes})
        
        print(dico_redondance)
        print(len(dico_redondance))
        for elt in dico_redondance :
            print(elt)
            print(dico_redondance[elt])
            
        
            
    
    with open("Resultats_sequences/graphs_all_2020_07_with_SSEs.pickle", 'rb') as fichier_pickle_1 :
        mon_depickler = pickle.Unpickler(fichier_pickle_1)
        graphes = mon_depickler.load()
        print("Nombre d'ARN en tout :")
        print(len(graphes))
        
        
        liste_non_trouvees = []
        liste_pas_repr = []
        liste_nom_graphes = []
        liste_plusieurs = []
        liste_nom_graphes_aminor = []
        
        
        for cle in liste_tout_aminor :
            with open("Nouvelles_donnees/fichier_%s_%s_10.pickle"%(cle[1][0], cle[1][1]), "rb") as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                extension = mon_depickler.load()
                
                if cle[1][0] == '4v67' :
                    print(extension.nodes[1]["num_ch"])
                    print(extension.nodes[3]["num_ch"])
                    print(extension.nodes[5]["num_ch"])
                
                if (cle[1][0].upper(), extension.nodes[1]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[1]["num_ch"]))
                
                if (cle[1][0].upper(), extension.nodes[3]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[3]["num_ch"]))
                
                if (cle[1][0].upper(), extension.nodes[5]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[5]["num_ch"]))
        
        print(len((liste_nom_graphes_aminor)))          
        #exit() 
        #a_faire = [liste_nom_graphes_aminor, graphes]
        liste_a_faire = graphes
        groupes_redondances_aminor = []
        
        
        c = 0
        for elt in liste_nom_graphes_aminor :
            existe = False
            for cle in dico_redondance :
                if cle == elt or elt in dico_redondance[cle] :
                    existe = True
            if not existe :
                c += 1
                print(elt)
        print(c)
        print(len(liste_nom_graphes_aminor))
        #exit()
        
        
#         for cle in dico_redondance :
#             cle_dans = False
#             if cle in liste_nom_graphes_aminor :
#                 groupes_redondances_aminor.append([cle])
#                 cle_dans = True
#             else :
#                 for elt in dico_redondance[cle] :
#                     if elt in liste_nom_graphes_aminor :
#                         if cle_dans :
#                             groupes_redondances_aminor[len(groupes_redondances_aminor)-1].append(elt)
#                            
#                         else :
#                             groupes_redondances_aminor.append([elt])
#                             cle_dans = True
#         
#         print(groupes_redondances_aminor)
#         for elt in groupes_redondances_aminor :
#             print(len(elt))
#             print(elt)
#         print(len(groupes_redondances_aminor))
        #exit()
        compteur = 0
        compte_false = 0
        compte_true = 0
        compteur_ramousnif = 0
        groupe_a_minor= []
        #for liste in a_faire :
        for cle in liste_a_faire:
            if cle not in liste_nom_graphes_aminor :
                ## si on veut ajouter une contrainte de résolution
                res = resolutions[cle[0]]
                #res = 2.0
                
                
                if res != None and res <= 3.0 and orgs[cle] != None :
                #if cle == ('6QZP', 'L5') :
                    print(cle)
                    trouve = False
                    trouve_repr = []
                    for elt in dico_redondance :
                        if '+' not in elt :
                            if cle == elt :
                                repr_dans_liste = False
                                for redondant in dico_redondance[elt] :
                                        if redondant in liste_nom_graphes_aminor :
                                            repr_dans_liste = True
                                
                                if not repr_dans_liste :
                                    trouve = True
                                else : 
                                    groupe_a_minor.append(elt)
                            elif cle in dico_redondance[elt] :
                                trouve_repr.append(elt)
                        else : ## cas des classes au représentant multiple
                            elt_split = []
                            for elt_2 in elt.split('+') :
                                elt_split.append((elt_2.split("|")[0], elt_2.split('|')[2]))
                            if cle in elt_split :
                                repr_dans_liste = False
                                for redondant in dico_redondance[elt] :
                                    for redondant2 in redondant :
                                        if redondant2 in liste_nom_graphes_aminor :
                                            repr_dans_liste = True
                                if not repr_dans_liste :
                                    trouve = True
                                    print("riiiiiit")
                                else : 
                                    groupe_a_minor.append(elt)
                                        
                            
                            elif cle in dico_redondance[elt] :
                                trouve_repr.append(elt_split)
                            else :
                                for redondant in dico_redondance[elt] :
                                    if cle in redondant :
                                        trouve_repr.append(elt_split)
                                
                
                    print(trouve)
                    print(trouve_repr)
                    
                    if trouve :
                        if cle not in liste_nom_graphes :
                            liste_nom_graphes.append(cle)
                    elif len(trouve_repr) == 0 and cle not in groupe_a_minor :
                        print("aaaaaah")
                        #if cle in liste_plus :
                        liste_non_trouvees.append(cle)
#                         if compteur == 0 :
#                             liste_nom_graphes.append(cle)
                        #exit()
                    elif len(trouve_repr) == 1 :
                        
                        if isinstance(trouve_repr[0], tuple) :
                            
                            if trouve_repr[0] not in graphes.keys() :
                                liste_pas_repr.append(cle)
                        else :
                            ok = True
                            for elt in trouve_repr[0] :
                                if elt not in graphes.keys() :
                                    ok = False
                            if not ok :
                                liste_pas_repr.append(cle)
                    else :
                        print("dans plusieurs")
                        liste_plusieurs.append(cle)
            #compteur += 1
        
#         for elt in liste_nom_graphes_aminor :
#             if elt not in liste_nom_graphes :
#                 liste_nom_graphes.append(elt)
#             else :
#                 print(elt)
                    
        ## les représentants des classes de notre jeu de données
        print(liste_nom_graphes)
        print(len(liste_nom_graphes))
        
        c = 0
        for elt in liste_nom_graphes_aminor :
            if elt in liste_nom_graphes :
                c += 1
                print(elt)
        print(c)
        
        print(groupe_a_minor)
        print(len(groupe_a_minor))
        
        c = 0
        for elt in liste_nom_graphes_aminor :
            if elt in groupe_a_minor :
                c += 1
                print(elt)
        print(c)
        
        ## les structures de notre jeu de données non trouvées dans les classes de la base de données
        print(liste_non_trouvees)
        print(len(liste_non_trouvees))
        print(('5MDV', 2) in liste_non_trouvees)
        print(('5MDV', 2) in liste_nom_graphes)
        print(('5MDV', 2) in liste_pas_repr)
#         c = 0
#         m = 0
#         for elt in liste_non_trouvees :
#             if elt in liste_nom_graphes_aminor :
#                 
#                 trouve1 = False
#                 trouve2 = False
#                 k = 0
#                 while k < len(dico_redondance) and not trouve :
#                         cle = list(dico_redondance.keys())[k]
#                         if cle == elt :
#                                 trouve1 = True
#                         if elt in dico_redondance[cle] :
#                                 trouve2 = True
#                         k += 1
#                 if not trouve1 :
#                     print("rat")    
#                     print(elt)
#                     c += 1
#                 else :
#                     if elt in liste_nom_graphes :
#                         m += 1
#                 
#                 print(elt)
#         print(c)
#         print(m)
        
        ## les structures de notre jeu de données pour lesquelles le représentant de la classe dans la base de données n'existe pas dans le jeu de données
        print(liste_pas_repr)
        print(len(liste_pas_repr))
        c = 0
        for elt in liste_pas_repr:
            if elt in liste_nom_graphes_aminor :
                c += 1
                
        print(c)
        
        ## les structures de notre jeu de données apparaissant dans plusieurs classes (il n'y en a pas)
        print(liste_plusieurs)
        print(len(liste_plusieurs))
        
#         liste_tot = []
#         liste_tot.extend(liste_nom_graphes)
#         liste_tot.append(('2VQE', 'A'))
        
        print(compte_false)
        print(compte_true)
        print(compteur_ramousnif)
        #liste_tot.extend(liste_non_trouvees)
        #liste_tot.extend(liste_pas_repr)
        
        return liste_nom_graphes

def recup_liste_pos_aminor():
    ''' 
        Renvoie les positions des motifs A-minor dans les structures de la PDB
    '''
    with open("all_aminor.pickle", 'rb') as fichier_all_aminor :
        mon_depickler = pickle.Unpickler(fichier_all_aminor)
        all_aminor = mon_depickler.load()
        
        liste_pos_aminor = []
        for cle in all_aminor :
            #print(cle)
            for graphe in all_aminor[cle] :
                
                positions_ajoutees = []
                graphe_motif = nx.MultiDiGraph()
                graphe_motif, positions_ajoutees = stockage_motif(graphe_motif, graphe, positions_ajoutees)
                
                liste_pos_graphe = (cle, [])
                
                liste_pos_graphe[1].append((graphe_motif.nodes[3]["num_ch"], graphe_motif.nodes[3]["position"][0], 1))
                liste_pos_graphe[1].append((graphe_motif.nodes[2]["num_ch"], graphe_motif.nodes[2]["position"][0], 2))
                liste_pos_graphe[1].append((graphe_motif.nodes[5]["num_ch"], graphe_motif.nodes[5]["position"][0], 3))

                liste_pos_aminor.append(liste_pos_graphe)
                print(liste_pos_graphe)
            
    return liste_pos_aminor

def alignement_global(cluster, liste_sequence_idem, num_seq, graphes):
    '''
        calcule les alignements globaux entre les séquences autour de chacun des motifs contenus dans le cluster 
        et les séquences autour des occurrences trouvées dans les structures de la PDB 
        et renvoie le nombre d'occurrences qui ont un score inférieur à la valeur qu'on veut (ici 80)
        
        :param cluster: une liste de noms de graphes d'extension contenus dans un même cluster
        :param liste_sequence_idem: une liste d'occurrences de séquences trouvées dans les structures de la PDB (num_pdb, num_chaine et position de l'occurrence)
        :param num_seq: numéro du brin considéré (1, 2 ou 3)
        :param graphes: un dictionnaire de toutes les structures de la PDB sous forme de graphes
        
        :type cluster: liste de tuples
        :type liste_sequence_idem : liste de tuples
        :type num_seq: int
        :type graphes: dictionnaire de nx.DiGraph
        
    '''
    nb_aln_global_inf_60  = []
    
    
    for elt in cluster : 
        nb_aln_global_inf_60.append(0)
        
        with open("Nouvelles_donnees/fichier_%s_%s_10.pickle"%(elt[0], elt[1]), 'rb') as fichier_graphe1 :
            mon_depickler = pickle.Unpickler(fichier_graphe1)
            extension = mon_depickler.load()
            if "fichier_%s_%s_%s.fa"%(elt[0],extension.nodes[num_seq+2]["num_ch"], extension.nodes[num_seq+2]["position"][0]) not in os.listdir("Resultats_sequences/aln_global/") :

            
    
                seqs = recup_sequences_autour_motif(elt[0], extension, 30, [[3,1], [2,4], [5,5]])
                print(seqs)
                
                with open("Resultats_sequences/aln_global/fichier_%s_%s_%s.fa"%(elt[0],extension.nodes[num_seq+2]["num_ch"], extension.nodes[num_seq+2]["position"][0]), 'w') as fichier_seq1 :
                        fichier_seq1.write(">%s_%s_%s\n"%(elt[0],extension.nodes[num_seq+2]["num_ch"], extension.nodes[num_seq+2]["position"][0]))
                        fichier_seq1.write(seqs[num_seq-1])
                
        for elt2 in liste_sequence_idem :
           

                    print(elt2)
                    for occ in elt2[1] :
                        if "fichier_%s_%s_%s.fa"%(elt2[0][0],elt2[0][1],occ) not in os.listdir("Resultats_sequences/aln_global/") :
                            
                            graphe = graphes[elt2[0]]
                            print(graphe.nodes())
                            print(occ)
                            #print(graphe.nodes[list(graphe.nodes())[0]])
                            #print(graphe.nodes.data())
                            seq = ""
                            for i in range(27, -1, -1) :
                                #print(i)
                                if occ - i > 0 and occ-i in graphe.nodes()  : 
                                    seq += graphe.nodes[occ-i]["nt"]
                                        
                            for i in range(1, 33) :
                                if occ + i < graphe.number_of_nodes() and occ+i in graphe.nodes()  : 
                                    seq += graphe.nodes[occ+i]["nt"]
                                
                            print(seq)
                            
                            with open("Resultats_sequences/aln_global/fichier_%s_%s_%s.fa"%(elt2[0][0],elt2[0][1],occ), 'w') as fichier_seq2 :
                                fichier_seq2.write(">%s_%s_%s\n"%(elt2[0][0],elt2[0][1],occ))
                                fichier_seq2.write(seq)
                                
                                
                        if "needle_%s_%s_%s_%s_%s_%s.txt"%(elt[0],extension.nodes[num_seq+2]["num_ch"], extension.nodes[num_seq+2]["position"][0], elt2[0][0],elt2[0][1],occ) not in os.listdir("Resultats_sequences/aln_global/") :               
                            needle_cline_1 = NeedleCommandline(r"/usr/local/emboss/bin/needle", asequence="Resultats_sequences/aln_global/fichier_%s_%s_%s.fa"%(elt[0],extension.nodes[num_seq+2]["num_ch"], extension.nodes[num_seq+2]["position"][0]), bsequence="Resultats_sequences/aln_global/fichier_%s_%s_%s.fa"%(elt2[0][0],elt2[0][1],occ), gapopen=10, gapextend=0.5, outfile="Resultats_sequences/aln_global/needle_%s_%s_%s_%s_%s_%s.txt"%(elt[0],extension.nodes[num_seq+2]["num_ch"], extension.nodes[num_seq+2]["position"][0], elt2[0][0],elt2[0][1],occ))                
                            stdout_1, stderr_1 = needle_cline_1()
                                
                        score = recup_score("Resultats_sequences/aln_global/needle_%s_%s_%s_%s_%s_%s.txt"%(elt[0],extension.nodes[num_seq+2]["num_ch"], extension.nodes[num_seq+2]["position"][0], elt2[0][0],elt2[0][1],occ))              
                        if score <= 80 :
                            nb_aln_global_inf_60[len(nb_aln_global_inf_60)-1] += 1
    
    return min(nb_aln_global_inf_60)
            

def sequence_aleatoire(graphes, liste_non_redondant):

        seq_alea = ""
        pas_trouve = True
        while pas_trouve :
            nb_cle = randint(0, len(liste_non_redondant)-1)
            cle = liste_non_redondant[nb_cle]
            
            #if graphes[cle].number_of_nodes() > 1000 :
            print(cle)
            if graphes[cle].number_of_nodes() >= 8 : 
                nb_al = randint(1, graphes[cle].number_of_nodes() - 7)
                print(graphes[cle].nodes())
                for i in range(nb_al, nb_al+8) :
                    print(i)
                    print(graphes[cle].nodes[i])
                    seq_alea += graphes[cle].nodes[i]["nt"]
                
                pas_trouve = False
                #break
                    
        print(seq_alea)
        return seq_alea, cle, nb_al
        
def donne_seq(graphes, liste_non_redondant):
    dico_seq = {}
    for cle in graphes :
            #print(cle)
            if cle in liste_non_redondant  :
                print(cle)
                print(graphes[cle].number_of_nodes())
                seq = ""
                for i in range(1, graphes[cle].number_of_nodes()+1) :
                    seq += graphes[cle].nodes[i]["nt"]
                
                print(seq)
                dico_seq.update({cle : seq})
    return dico_seq

def recherche_seq_alea(seq_alea, dico_seq_tot, cle_alea, nb_al):
        
        dico_seq = {}
        
        for cle in dico_seq_tot :
            #print(cle)
                seq = dico_seq_tot[cle]
#                 print(seq[nb_al:])
#                 print(nb_al)
                
                occ = []
                start = -1
                while start != 0 :
                    start = seq.find(seq_alea, max(start,0))
                    if start != -1 :
                        occ.append(start)
                    start += 1
                if len(occ) > 0 :
                    dico_seq.update({cle : occ})
                
        
        compte = 0
        for elt in dico_seq :
            if len(dico_seq[elt]) > 0 :
                compte += len(dico_seq[elt])
            
        return compte, dico_seq


def recherche_tot(graphes, liste_non_redondant, clusters, liste_aminor):
    with open("Resultats_sequences/nrlist_3.136_all.csv", 'r') as fichier_csv_motifs_atlas:
        csvreader = csv.reader(fichier_csv_motifs_atlas)
        dico_redondance = {}
        for row in csvreader :
            if '+' not in row[1] : ## cas de base
                
                repr = (row[1].split("|")[0], row[1].split("|")[2])
                
                representes = []
                autres = row[2].split(",")
                for elt in autres :
                    representes.append((elt.split("|")[0], elt.split("|")[2]))
                dico_redondance.update({repr:representes})
            else : ## cas où le représentant de la classe est composé de plusieurs structures (je sais pas pourquoi ?)
                compteur = 0
                repr = tuple(row[1])
                bizarre = False
                representes = []
                autres = row[2].split(",")
                for elt in autres :
                    for elt_2 in elt.split('+') :
                        representes.append((elt_2.split("|")[0], elt_2.split("|")[2]))

                dico_redondance.update({repr:representes})
    
    with open("sequences_par_cluster.txt", 'r') as fichier_seq :

        communs_aminor = []
        communs_non_aminor = []
        with open("Resultats_sequences/fichier_csv_nb_seq_non_aminor_par_sequence_temp.csv", 'w', newline="") as fichier_csv :
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["Numéro cluster", "Sequences", "Elts cluster et position aminor", "Nb seq 1", "Nb seq 2", "Nb seq 3", "Non A-minor", "Nb 3 seq"])
       
            with open("Resultats_sequences/fichier_csv_nb_seq_aminor_par_sequence_temp.csv", 'w', newline="") as fichier_csv2 :
                csvwriter2 = csv.writer(fichier_csv2)
                csvwriter2.writerow(["Numéro cluster", "Sequences", "Elts cluster et position aminor", "Nb seq 1", "Nb seq 2", "Nb seq 3", "A-minor", "Nb 3 seq"])
    
              
                dico_seq_tot = donne_seq(graphes, liste_non_redondant)
                
                compte_num = 0
                
                ligne = fichier_seq.readline()
                nom_cluster = ""
                while ligne != "" :
                    print(ligne)
                    if ligne[0] == '>' :
                        print(compte_num)
                        nom_cluster = ligne[1:len(ligne)-1]
                    
                        num_cluster = int(nom_cluster.split("_")[1])
                        #if num_cluster == 4 : exit()
                        liste_cluster = []
                        for aminor in clusters[num_cluster-1] :
                            with open("Nouvelles_donnees/fichier_%s_%s_10.pickle"%(aminor[0], aminor[1]), "rb") as fichier_ext :
                                mon_depickler = pickle.Unpickler(fichier_ext)
                                extension = mon_depickler.load()
                            
                                liste_cluster.append((aminor, (extension.nodes[1]["num_ch"], extension.nodes[1]["position"][0]), (extension.nodes[2]["num_ch"], extension.nodes[2]["position"][0]), (extension.nodes[5]["num_ch"], extension.nodes[5]["position"][0])))
    
                        liste_seqs = []
                        fichier_seq.readline()
                        ligne = fichier_seq.readline()
                        while ligne != "" and ligne[0] != '>' :
                            
                            seqs = ligne.split(", ")
                            seqs[2] = seqs[2][:len(seqs[2])-1]
                            print(seqs)
                            
                            deja_fait = False
                            for sequence in liste_seqs :
                                if sequence[0] == seqs[0] and sequence[1] == seqs[1] and sequence[2] == seqs[2] :
                                    deja_fait = True
                            
                            if not deja_fait :
                                
                                liste_seqs.append(seqs)
                                
                                aminors = [[],[],[]]
                                non_aminors = [[],[],[]]
                                nb_trouve = [0,0,0]
                                nb_trouve_non_aminor = [0,0,0]
                                dico_seqs_3_aminor = [[],[],[]]
                                dico_seqs_3_non_aminor = [[],[],[]]
                                
                                
                                for i in range(3) :
                                    compte_occ, dico_seqs = recherche_seq_alea(seqs[i], dico_seq_tot, nom_cluster, -1)
                                
                                        
                                    
                                    #dico_seqs_3[i].extend(list(dico_seqs.keys()))
        #                             if i == 0 :
        #                                 print(compte_occ)
        #                                 for cle in dico_seqs :
        #                                     if len(dico_seqs[cle]) > 0 :
        #                                         print(cle)
        #                                         print(dico_seqs[cle])
                                        #print(dico_seqs)
                                    
                                    for occ in dico_seqs :
                                        print(nb_trouve)
                                        print(nb_trouve_non_aminor)
                                        print(len(dico_seqs[occ]))
                                        print(occ)
                                        print(dico_seqs[occ])
                                        #break
                                        ancien_nb_trouve_i = nb_trouve[i]
                                        liste_a_chercher = []
                                        trouve = False
                                        for elt in dico_redondance :
                                            if occ == elt or occ in elt :
                                                liste_a_chercher.append(dico_redondance[elt])
                                                trouve = True
                                            elif occ in dico_redondance[elt] :
                                                liste_a_chercher.append(dico_redondance[elt])
                                                trouve = True
                                        if not trouve :
                                            liste_a_chercher.append([occ])
                                        print("liste a chercher")
                                        print(len(liste_a_chercher))
                                        #exit()
                                        if i == 0 and num_cluster == 2 and occ == ('4V4B', 'B3') :
                                                print("estla")
                                                print(occ)
                                                print(liste_a_chercher)
                                                print(dico_seqs[occ])
                                                #print(('5DMV', 2) in liste_a_chercher)
                                        
                                        for liste in liste_a_chercher : 
                                            nb_trouve_temp = 0
                                            aminors_temp = []
                                            existe = False
                                            for elt in liste :
                                                if i == 0 and num_cluster == 2 and occ == ('4V4B', 'B3') :
                                                            print("estici")
                                                            print(occ)
                                                            print(liste_a_chercher)
                                                            print(dico_seqs[occ])
                                                            print(elt)
                                                            print(occ == elt)
                                                if occ != elt : 
                                                    
                                                    if elt in graphes.keys() :
                                                        seq = ""
                                                        for j in range(1, graphes[elt].number_of_nodes()+1) :
                                                            seq += graphes[elt].nodes[j]["nt"]
                                                        
                                                        #print(seq)
                                                        
                                                        occ_temp = []
                                                        start = -1
                                                        while start != 0 :
                                                            start = seq.find(seqs[i], max(start,0))
                                                            if start != -1 :
                                                                occ_temp.append(start)
                                                            start += 1
                                                            
                                                        if i == 1 and num_cluster == 3 and len(occ_temp) > 1 :
                                                            print("estici")
                                                            print(occ_temp)
                                                            
                                                            print(occ)
                                                            print(liste_a_chercher)
                                                            print(dico_seqs[occ])
                                                            print(elt)
                                                            print(occ == elt)
                                                        if len(occ_temp) > 0 :
                                                            
                                                            
                                                            for aminor in liste_aminor :
                                #                                 print(aminor)
                                #                                 print(occ)
                                                                if elt[0] == aminor[0].upper() :
                                                                         
                                                                        if elt[1] == aminor[1][i][0] :
                                                                            if i == 0 and num_cluster == 2 and occ == ('4V4B', 'B3') :
                                                                                print("estici")
                                                                                print(occ)
                                                                                print(liste_a_chercher)
                                                                                print(dico_seqs[occ])
                                                                            for occ2 in occ_temp :
                                                                                if occ2 >= aminor[1][i][1] - 5 and occ2 <= aminor[1][i][1] :
                                                                                    if i == 0 and num_cluster == 2 and occ == ('4V4B', 'B3') :
                                                                                        print("estencorela")
                                                                                        print(occ)
                                                                                        print(liste_a_chercher)
                                                                                        print(dico_seqs[occ])
                                                                                    
                                                                                    if (elt, occ2, seqs[i]) not in aminors_temp :
                                                                                        nb_trouve_temp += 1
                                                                                        aminors_temp.append((elt, occ2, seqs[i]))
                                                                                        #if elt not in dico_seqs_3_aminor[i] :
                                                                                        dico_seqs_3_aminor[i].append(elt)

                                                                                            
                                                                                    existe = True
                                                else :
                                                        if i == 0 and num_cluster == 2 and occ == ('4V4B', 'B3') :
                                                            print("estici")
                                                            print(occ)
                                                            print(liste_a_chercher)
                                                            print(dico_seqs[occ])
                                                            print(elt)
                                                        for aminor in liste_aminor :
                            #                                 print(aminor)
                            #                                 print(occ)
                                                            if elt[0] == aminor[0].upper() :
                                                                     
                                                                    if elt[1] == aminor[1][i][0] :
                                                                        
                                                                        for occ2 in dico_seqs[occ] :
                                                                            if occ2 >= aminor[1][i][1] - 5 and occ2 <= aminor[1][i][1] :
                                                                                if i == 0 and num_cluster == 2 and occ == ('4V4B', 'B3') :
                                                                                    print("estencorela")
                                                                                    print(occ)
                                                                                    print(liste_a_chercher)
                                                                                    print(dico_seqs[occ])
                                                                                
                                                                                if (occ, occ2, seqs[i]) not in aminors_temp :
                                                                                    nb_trouve_temp += 1
                                                                                    aminors_temp.append((occ, occ2, seqs[i]))
                                                                                    #if occ not in dico_seqs_3_aminor[i] :
                                                                                    dico_seqs_3_aminor[i].append(occ)
        
                                                                                existe = True
                                                if existe :
                                                    nb_trouve[i] += nb_trouve_temp
                                                    aminors[i].extend(aminors_temp)
                                                    break
                                        if nb_trouve[i] == ancien_nb_trouve_i :
                                            for occ2 in dico_seqs[occ] :
                                                if (occ, occ2, seqs[i]) not in non_aminors[i] :
                                                    non_aminors[i].append((occ, occ2, seqs[i]))
                                                    #if occ not in dico_seqs_3_non_aminor :
                                                    dico_seqs_3_non_aminor[i].append(occ)
                                    nb_trouve_non_aminor[i] = compte_occ - nb_trouve[i] 
                                    
                                    if i == 1 and num_cluster == 3 and nb_trouve_non_aminor[i] == -1 :
                                        print("lalaa")
                                        print(dico_seqs)
                                        print(aminors[i])
                                        #exit()
                                                
                                aminors_temp = list(aminors)
                                
                                commun_aminor = []
                                for cle1 in aminors[0] :
                                    for cle2 in aminors[1] :
                                        for cle3 in aminors[2] :
                                            if cle1[0] == cle2[0] and cle2[0] == cle3[0] :
                                                commun_aminor.append((cle1[0], cle1[1], cle2[1], cle3[1]))
                                                aminors_temp[0].remove(cle1)
                                                aminors_temp[1].remove(cle2)
                                                aminors_temp[2].remove(cle3)
                                                
                                
                                non_aminors_temp = list(non_aminors)
                                non_aminors_temp[0].extend(aminors_temp[0])
                                non_aminors_temp[1].extend(aminors_temp[1])
                                non_aminors_temp[2].extend(aminors_temp[2])
                                
                                commun_non_aminor = []
                                for cle1 in non_aminors[0] :
                                    for cle2 in non_aminors[1] :
                                        for cle3 in non_aminors[2] :
                                            if cle1[0] == cle2[0] and cle2[0] == cle3[0] :
                                                commun_non_aminor.append((cle1[0], cle1[1], cle2[1], cle3[1]))   
                                                
                                                  
                                        
                                #commun_aminor = [x for x in dico_seqs_3_aminor[0] if x in dico_seqs_3_aminor[1] and x in dico_seqs_3_aminor[2]]
                                #commun_non_aminor = [x for x in dico_seqs_3_non_aminor[0] if x in dico_seqs_3_non_aminor[1] and x in dico_seqs_3_non_aminor[2]]

                                communs_aminor.append(commun_aminor)
                                communs_non_aminor.append(commun_non_aminor)
                                
                                commun_aminor_temp = list(commun_aminor)
                                commun_non_aminor_temp = list(commun_non_aminor)
                                
                                ## recherche la partie stem dans les occurrences identifiées comme A-minor
                                
                                c = 0
                                for cle in commun_aminor :
                                    for k in range(3,6) :
                                        for l in range(3,6) :
                                            if 's' not in commun_aminor[c] and (("Stem" in graphes[cle[0]].nodes[cle[2]+k]["part"]  and "Stem" in graphes[cle[0]].nodes[cle[3]+l]["part"] and \
                                                graphes[cle[0]].nodes[cle[2]+k]["part_id"][graphes[cle[0]].nodes[cle[2]+k]["part"].index("Stem")] == graphes[cle[0]].nodes[cle[3]+l]["part_id"][graphes[cle[0]].nodes[cle[3]+l]["part"].index("Stem")]) or \
                                                ((cle[2]+k, cle[3]+l) in graphes[cle[0]].edges() and  graphes[cle[0]].edges[cle[2]+k, cle[3]+l]["label"] == 'CWW') ) :
                                                commun_aminor[c] += tuple("s")
                                    c += 1
                                                            
                                
                                ## recherche la partie stem dans les occurrences identifiées comme non A-minor (la position doit être exacte ici)
                                                   
                                c = 0
                                for cle in commun_non_aminor :
                                    for k in range(4,5):
                                        for l in range(4,5) :
                                            if 's' not in commun_non_aminor[c] and (("Stem" in graphes[cle[0]].nodes[cle[2]+k]["part"]  and "Stem" in graphes[cle[0]].nodes[cle[3]+l]["part"] and \
                                                graphes[cle[0]].nodes[cle[2]+k]["part_id"][graphes[cle[0]].nodes[cle[2]+k]["part"].index("Stem")] == graphes[cle[0]].nodes[cle[3]+l]["part_id"][graphes[cle[0]].nodes[cle[3]+l]["part"].index("Stem")]) or \
                                                ((cle[2]+k, cle[3]+l) in graphes[cle[0]].edges() and  graphes[cle[0]].edges[cle[2]+k, cle[3]+l]["label"] == 'CWW') ) :
                                                commun_non_aminor[c] += tuple("s")
                                    c += 1
                                
                                ## recherche la partie boucle dans les occurrences identifiées comme A-minor
                                c = 0
                                for cle in commun_aminor :            
                                        for k in range(3,6) :
                                            voisin_cww_1 = False
                                            for noeud_voisin in graphes[cle[0]][cle[1]+k] :
                                                if graphes[cle[0]].edges[cle[1]+k, noeud_voisin]["label"] == "CWW" :
                                                    voisin_cww_1 = True
                                            voisin_cww_2 = False
                                            for noeud_voisin in graphes[cle[0]][cle[1]+k+1] :
                                                if graphes[cle[0]].edges[cle[1]+k+1, noeud_voisin]["label"] == "CWW" :
                                                    voisin_cww_2 = True
                                            if (not voisin_cww_1 and not voisin_cww_2) and 'l' not in commun_aminor[c]   :
                                                commun_aminor[c] += tuple("l")
                                                fait = True
                                                print("bouhh")
                                        c += 1
                                                #exit()
                                                
                                ## recherche la partie stem dans les occurrences identifiées comme non A-minor (doit être exacte ici)
                                
                                c = 0
                                for cle in commun_non_aminor :             
                                        for k in range(4,5) :
                                            voisin_cww_1 = False
                                            for noeud_voisin in graphes[cle[0]][cle[1]+k] :
                                                if graphes[cle[0]].edges[cle[1]+k, noeud_voisin]["label"] == "CWW" :
                                                    voisin_cww_1 = True
                                            voisin_cww_2 = False
                                            for noeud_voisin in graphes[cle[0]][cle[1]+k+1] :
                                                if graphes[cle[0]].edges[cle[1]+k+1, noeud_voisin]["label"] == "CWW" :
                                                    voisin_cww_2 = True
                                            if (not voisin_cww_1 and not voisin_cww_2) and 'l' not in commun_non_aminor[c]   :
                                                commun_non_aminor[c] += tuple("l")
                                                fait = True
                                                print("bouhh")
                                        c += 1
                                                #exit()
                                                    #exit()
                                
                                #print(communs)       
                                #exit()
                                csvwriter2.writerow([nom_cluster, seqs, liste_cluster, nb_trouve[0], nb_trouve[1], nb_trouve[2], aminors, len(commun_aminor)])
                                csvwriter.writerow([nom_cluster, seqs, liste_cluster, nb_trouve_non_aminor[0], nb_trouve_non_aminor[1], nb_trouve_non_aminor[2], non_aminors, len(commun_non_aminor)])
                        
                            ligne = fichier_seq.readline()
            
                        compte_num += 1
                    #break
                #break
#         print(communs)
#         for elt in communs :
#             print(len(elt))

        with open("Resultats_sequences/communs_aminor.pickle", 'wb') as fichier_aminor :
            mon_pickler = pickle.Pickler(fichier_aminor)
            mon_pickler.dump(communs_aminor)
        with open("Resultats_sequences/communs_non_aminor.pickle", 'wb') as fichier_non_aminor :
            mon_pickler = pickle.Pickler(fichier_non_aminor)
            mon_pickler.dump(communs_non_aminor)    

def lecture_fichier_csv(nom_fichier_csv):
    with open(nom_fichier_csv, 'r') as fichier_csv :
        csvreader = csv.reader(fichier_csv, delimiter=',')
        dico_val = {}
        compteur = 0
        for row in csvreader :
            if row[0][0] == "c" :
                seqs = [row[1].split("', '")[0][2:], row[1].split("', '")[1], row[1].split("', '")[2][:len(row[1].split("', '")[2])-2]] 
                dico_val.update({(row[0], compteur) : [seqs, row[3], row[4], row[5]]})
                compteur += 1
        return dico_val
       
def diff_occ_aminor(nom_fichier_csv1, nom_fichier_csv2):
        dico_val_1 = lecture_fichier_csv(nom_fichier_csv1)
        dico_val_2 = lecture_fichier_csv(nom_fichier_csv2)
        
        compter = 0
        for cle1 in dico_val_1 :
            for cle2 in dico_val_2 :
                if cle1[0] == cle2[0] :
                    if dico_val_1[cle1][0][0] == dico_val_2[cle2][0][0] and dico_val_1[cle1][0][1] == dico_val_2[cle2][0][1] and dico_val_1[cle1][0][2] == dico_val_2[cle2][0][2] :
                        if dico_val_1[cle1][1] > dico_val_2[cle2][1] or dico_val_1[cle1][2] > dico_val_2[cle2][2] or dico_val_1[cle1][3] > dico_val_2[cle2][3] :
                            print(cle1)
                            print(dico_val_1[cle1])
                            print(cle2)
                            print(dico_val_2[cle2])
                            compter += 1
        print(compter)

def edge_match(d_g1, d_g2):
    if d_g1['label'] == d_g2['label'] :
        return True
    return False


def recup_org_nat(num_pdb, num_ch):
    org = None  
    file = PATH_MMCIF+ num_pdb+".cif"
    if num_pdb+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
        print("petit rat")
        url = 'http://www.rcsb.org/pdb/files/%s.cif' % num_pdb
        request.urlretrieve(url, file)
        print("houhou")
    try :     
        doc = cif.read_file(PATH_MMCIF+num_pdb+".cif")
        block = doc.sole_block()
        
        cat = block.find_mmcif_category("_entity_poly.") 
        #print(list(cat.tags))
        id = -1
        for r in cat :
            #print(tuple(row))
            num_chaine = r[6].split(",")
            
#                             print("petit rat")
#                             print(num_ch)
            if num_ch in num_chaine :
                id = r[0]
        
        
        cat = block.find_mmcif_category("_entity_src_nat.")         
        
             
        for r in cat :
            if r[0] == id :
                org = r[6]
        
        return org
    except RuntimeError :
        print("probleme de fichier cif : %s"%num_pdb)
        return -1
        
def recup_org_gen(num_pdb, num_ch):
    org = None  
    file = PATH_MMCIF+ num_pdb+".cif"
    if num_pdb+".cif" not in os.listdir(PATH_MMCIF) : ## si le fichier n est pas encore stocke en interne, on va le chercher sur la PDB et on le stocke ou on veut
        print("petit rat")
        url = 'http://www.rcsb.org/pdb/files/%s.cif' % num_pdb
        request.urlretrieve(url, file)
        print("houhou")
    try :     
        doc = cif.read_file(PATH_MMCIF+num_pdb+".cif")
        block = doc.sole_block()
        
        cat = block.find_mmcif_category("_entity_poly.") 
        #print(list(cat.tags))
        id = -1
        for r in cat :
            #print(tuple(row))
            num_chaine = r[6].split(",")
            
#                             print("petit rat")
#                             print(num_ch)
            if num_ch in num_chaine :
                id = r[0]
        
        
        cat = block.find_mmcif_category("_entity_src_gen.")         
        
             
        for r in cat :
            if r[0] == id :
                org = r[6]
        
        return org
    except RuntimeError :
        print("probleme de fichier cif : %s"%num_pdb)
        return -1


def test(clusters, graphes):
    with open("Resultats_sequences/nrlist_3.136_all.csv", 'r') as fichier_csv_motifs_atlas:
        csvreader = csv.reader(fichier_csv_motifs_atlas)
        
        compteur = 0
        dico_redondance = {}
        for row in csvreader :
            print(compteur)
            compteur += 1
            if '+' not in row[1] : ## cas de base
                
                repr = (row[1].split("|")[0], row[1].split("|")[2])
                
                org = recup_org_nat(repr[0], repr[1])
                
                if org != None :
                    representes = []
                    autres = row[2].split(",")
                    for elt in autres :
                        org = recup_org_nat(elt.split("|")[0], elt.split("|")[2])
                        if org != None :
                            representes.append((elt.split("|")[0], elt.split("|")[2]))
                    dico_redondance.update({repr:representes})
            else : ## cas où le représentant de la classe est composé de plusieurs structures (je sais pas pourquoi ?)
                #compteur = 0
                repr = row[1]
                bizarre = False
                representes = []
                
                orgs = []
                for elt in repr.split('+') :
                    e = (elt.split("|")[0], elt.split("|")[2])
                    orgs.append(recup_org_nat(e[0], e[1]))
                    
                if None not in orgs :
                    autres = row[2].split(",")
                    for elt in autres :
                        representes_temp = []
                        for elt_2 in elt.split('+') :
                            org = recup_org_nat(elt_2.split("|")[0], elt_2.split("|")[2])
                            if org != None :
                                representes_temp.append((elt_2.split("|")[0], elt_2.split("|")[2]))
                        representes.append(representes_temp)
                    dico_redondance.update({repr:representes})
                    
    
        print(len(dico_redondance))
        print(compteur)
        exit()
    with open("Resultats_sequences/fichier_csv_expr_reg_test.csv", 'w', newline="") as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(["num_cluster", "expr", "cle dico redondance", "occ"])        
        
        compteur = 1
        #exit()
        for cluster in clusters :
            print(compteur)
            if len(cluster) > 2 and compteur == 3:
                print(len(cluster))
                print(cluster)
                
                ## calcule le graphe commun du cluster
                graphe_commun_global, dico_graphes, graphe_commun = recherche_graphe_commun([4], cluster, 0.6, "ramou", "Graphes_communs_mai_2020/extensions", compteur)
                
                #exit()
                #try :
                
                liste_seq_13, pos_aminor_seq_13, liste_seq_24, pos_aminor_seq_24 = etablissement_expr_reg(dico_graphes, graphe_commun)
                for liste_expr in liste_seq_13 :
                                print(liste_expr)
                                print(liste_seq_13[liste_expr])
                                print(len(liste_seq_13[liste_expr]))
                                print(liste_seq_24[liste_expr])
                                print(len(liste_seq_24[liste_expr]))
                                compte_tot = 0
                                for elt in dico_redondance :
                                    if '+' in elt :
                                        print(elt)
                                        elts = elt.split('+')
                                        liste_elts = []
                                        for e in elts :
                                            liste_elts.append((e.split("|")[0],e.split("|")[2]))
                                        i = 0
                                        for e in liste_elts :
                                            compte_13 = 0
                                            compte_24 = 0
                                            dico_trouves = {}
                                            dico_trouves.update({(compteur, liste_expr, e) : [[],[]]})
                                            dico_seq = donne_seq(graphes, [e])
                                            print(dico_seq)
                                            #return
                                            print(liste_seq_13[liste_expr])
                                            print(liste_seq_24[liste_expr])
                                            if e in dico_seq :
                                                match_it_13 = re.finditer(liste_seq_13[liste_expr], dico_seq[e])
                                                match_it_24 = re.finditer(liste_seq_24[liste_expr], dico_seq[e])
                                            
                                                for m in match_it_13 :
                                                    pos = m.span()
                                                    dico_trouves[(compteur, liste_expr,e)][0].append(pos[0])
                                                    print(pos)
                                                    
                                                for m in match_it_24 :
                                                    pos = m.span()
                                                    dico_trouves[(compteur, liste_expr,e)][1].append(pos[0])
                                                    print(pos)   
                                            
                                                csvwriter.writerow([compteur,liste_expr, liste_seq_13[liste_expr],liste_seq_24[liste_expr], e, dico_trouves[(compteur, liste_expr,e)][0], dico_trouves[(compteur, liste_expr,e)][1]  ])
                                            nb_trouves_13 = len(dico_trouves[(compteur, liste_expr,e)][0])
                                            nb_trouves_24 = len(dico_trouves[(compteur, liste_expr,e)][1])
                                            
                                            for elt2 in dico_redondance[elt] :
                                                
                                                    print(elt)
                                                    print(elt2)
                                                    dico_trouves.update({(compteur, liste_expr, elt2[i]) : [[],[]]})
                                                    dico_seq = donne_seq(graphes, [elt2[i]])
                                                    if elt2[i] in dico_seq :
                                                        match_it_13 = re.finditer(liste_seq_13[liste_expr], dico_seq[elt2[i]])
                                                        match_it_24 = re.finditer(liste_seq_24[liste_expr], dico_seq[elt2[i]])
                                                    
                                                        for m in match_it_13 :
                                                            pos = m.span()
                                                            dico_trouves[(compteur, liste_expr, elt2[i])][0].append(pos[0])
                                                            print(pos)
                                                            
                                                        for m in match_it_24 :
                                                            pos = m.span()
                                                            dico_trouves[(compteur, liste_expr, elt2[i])][1].append(pos[0])
                                                            print(pos)   
                                                    
                                                                                                        
                                                        if len(dico_trouves[(compteur, liste_expr,elt2[i])][0]) > nb_trouves_13 or  len(dico_trouves[(compteur, liste_expr,elt2[i])][1]) > nb_trouves_24:
                                                            
                                                            csvwriter.writerow([compteur,liste_expr, liste_seq_13[liste_expr],liste_seq_24[liste_expr], "", elt2[i], dico_trouves[(compteur, liste_expr,elt2[i])][0], dico_trouves[(compteur, liste_expr,elt2[i])][1], "Diff"])
                                                            if len(dico_trouves[(compteur, liste_expr,elt2[i])][0]) > nb_trouves_13 :
                                                                compte_13 += 1
                                                            if len(dico_trouves[(compteur, liste_expr,elt2[i])][1]) > nb_trouves_24 :
                                                                compte_24 += 1
                                                        else :
                                                            csvwriter.writerow([compteur,liste_expr, liste_seq_13[liste_expr],liste_seq_24[liste_expr], "", elt2[i], dico_trouves[(compteur, liste_expr,elt2[i])][0], dico_trouves[(compteur, liste_expr,elt2[i])][1]  ])
    
                                            if compte_13 >= 1 or compte_24 >= 1 : 
                                                compte_tot += 1
                                            i +=  1
                                                
                                    else :    
                                        compte_13 = 0
                                        compte_24 = 0
                                        dico_trouves = {}
                                        dico_trouves.update({(compteur, liste_expr, elt) : [[],[]]})
                                        dico_seq = donne_seq(graphes, [elt])
                                        print(dico_seq)
                                        #return
                                        print(liste_seq_13[liste_expr])
                                        print(liste_seq_24[liste_expr])
                                        if elt in dico_seq :
                                            match_it_13 = re.finditer(liste_seq_13[liste_expr], dico_seq[elt])
                                            match_it_24 = re.finditer(liste_seq_24[liste_expr], dico_seq[elt])
                                        
                                            for m in match_it_13 :
                                                pos = m.span()
                                                dico_trouves[(compteur, liste_expr,elt)][0].append(pos[0])
                                                print(pos)
                                                
                                            for m in match_it_24 :
                                                pos = m.span()
                                                dico_trouves[(compteur, liste_expr,elt)][1].append(pos[0])
                                                print(pos)   
                                        
                                            csvwriter.writerow([compteur,liste_expr, liste_seq_13[liste_expr],liste_seq_24[liste_expr], elt, dico_trouves[(compteur, liste_expr,elt)][0], dico_trouves[(compteur, liste_expr,elt)][1]  ])
                                        nb_trouves_13 = len(dico_trouves[(compteur, liste_expr,elt)][0])
                                        nb_trouves_24 = len(dico_trouves[(compteur, liste_expr,elt)][1])
                                        
                                        for elt2 in dico_redondance[elt] :
                                            dico_trouves.update({(compteur, liste_expr, elt2) : [[],[]]})
                                            dico_seq = donne_seq(graphes, [elt2])
                                            if elt2 in dico_seq :
                                                match_it_13 = re.finditer(liste_seq_13[liste_expr], dico_seq[elt2])
                                                match_it_24 = re.finditer(liste_seq_24[liste_expr], dico_seq[elt2])
                                            
                                                for m in match_it_13 :
                                                    pos = m.span()
                                                    dico_trouves[(compteur, liste_expr, elt2)][0].append(pos[0])
                                                    print(pos)
                                                    
                                                for m in match_it_24 :
                                                    pos = m.span()
                                                    dico_trouves[(compteur, liste_expr, elt2)][1].append(pos[0])
                                                    print(pos)   
                                            
                                                                                                
                                                if len(dico_trouves[(compteur, liste_expr,elt2)][0]) > nb_trouves_13 or  len(dico_trouves[(compteur, liste_expr,elt2)][1]) > nb_trouves_24:
                                                    
                                                    csvwriter.writerow([compteur,liste_expr, liste_seq_13[liste_expr],liste_seq_24[liste_expr], "", elt2, dico_trouves[(compteur, liste_expr,elt2)][0], dico_trouves[(compteur, liste_expr,elt2)][1], "Diff"])
                                                    if len(dico_trouves[(compteur, liste_expr,elt2)][0]) > nb_trouves_13 :
                                                        compte_13 += 1
                                                    if len(dico_trouves[(compteur, liste_expr,elt2)][1]) > nb_trouves_24 :
                                                        compte_24 += 1
                                                else :
                                                    csvwriter.writerow([compteur,liste_expr, liste_seq_13[liste_expr],liste_seq_24[liste_expr], "", elt2, dico_trouves[(compteur, liste_expr,elt2)][0], dico_trouves[(compteur, liste_expr,elt2)][1]  ])

                                        if compte_13 >= 1 or compte_24 >= 1 : 
                                            compte_tot += 1
                                
                                            
                                print(compte_tot)
                                print(len(dico_redondance))
                                        
                                return              
                                        
                
                #except : 
                #    print("expr regu ratee")
            compteur += 1
                
def recherche_expr_reg(clusters, liste_non_redondant, graphes, liste_aminor):
    ''' 
    recherche des expressions regulières correspondant aux représentants des extensions dans les séquences de la PDB
    
    :param clusters: la liste des clusters 
    :param liste_non_redondant: la liste des structures non redondantes de la PDB (selon RNA 3D motif Atlas)
    :param graphes: la liste des graphes de toutes les structures de la PDB
    :param liste_aminor : la liste de toutes les positions des A-minor dans toutes les structures ayant des A-minor dans la PDB

    :type clusters: liste de tuples de taille 2
    :type liste_non_redondant: liste de tuples de taille 2
    :type graphes: dictionnaire de nx.DiGraph
    :type liste_aminor: liste 
'''
    
    ## récupère les séquences de toutes les structures de liste_non_redondant
    liste_non_redondant_tot = []
    for liste in liste_non_redondant :
        liste_non_redondant_tot.extend(liste)
    dico_seq_tot = donne_seq(graphes, liste_non_redondant_tot)
    print(('4V9Q', 'DV') in dico_seq_tot)
    print(('4V9Q', 'DV') in graphes)
    #exit()
    
    
    with open("Resultats_sequences/fichier_csv_expr_reg_new_version.csv", 'a', newline="") as fichier_csv :
        csvwriter = csv.writer(fichier_csv)
        csvwriter.writerow(["Numéro cluster", "Expr reg1", "Expr reg2", "Elts cluster", "Nb seq 1 non A-minor", "Nb seq 2 non A-minor", "Non A-minor", "A-minor", "Nb 3 seq non A-minor", "Nb 3 seq A-minor", "Nb struct non A-minor", "Nb struct A-minor"])

        dico_trouves = {}
        dico_couples = {}
        compteur = 1
        #exit()
        for cluster in clusters :
            print(compteur)
            if len(cluster) > 2 :
                    print(len(cluster))
                    print(cluster)
                
                ## calcule le graphe commun du cluster
                    graphe_commun_global, dico_graphes, graphe_commun = recherche_graphe_commun([4], cluster, 0.6, "ramou", "Graphes_communs_mai_2020/extensions", compteur)
                
                #exit()
                #try :
                    ## calcule les expressions régulières
                    liste_seq_13, pos_aminor_seq_13, liste_seq_24, pos_aminor_seq_24 = etablissement_expr_reg(dico_graphes, graphe_commun)
                    print("ramousnif")
                    print(liste_seq_13)
                    print(liste_seq_24)
                    #exit()
    #                 print(graphe_commun_global.nodes[8])
    #                 print(graphe_commun_global.nodes[2])
    #                 print(graphe_commun_global.nodes[4])
    #                 print(graphe_commun_global.nodes[11])
    #                 print(graphe_commun_global.nodes[6])
                    
                    ## récupère les numéros des noeuds du graphe commun global par chaine 
                    chaines_1 = [-1]*8
                    chaines_2 = [-1]*8
                    
                    for noeud, data in graphe_commun_global.nodes(data=True) :
                        if len(data["position"]) == 1 :
                            if 1 in data["chaine"] :
                                chaines_1[data["position"][0]+4] = noeud
                        
                            elif 3 in data["chaine"] :
                                chaines_1[data["position"][0]-7] = noeud
                            
                            elif 2 in data["chaine"] :
                                chaines_2[data["position"][0]-7] = noeud
                                
                            elif 4 in data["chaine"] :
                                chaines_2[data["position"][0]+4] = noeud
                    
                    
                    print(chaines_1)
                    print(chaines_2)
                    
                    #exit()
                   
                    for liste_expr in liste_seq_13 :
                        #for g in range(len(liste_seq_13[liste_expr])) :
                            #for k in range(len(liste_seq_24[liste_expr])) :
                                ## recherche les occurrences des expressions régulières
                            print(liste_expr)
                            print(liste_seq_13[liste_expr])
                            print(len(liste_seq_13[liste_expr]))
                            print(liste_seq_24[liste_expr])
                            print(len(liste_seq_24[liste_expr]))
                            #exit()
                            
                            c = 0    
                            nb_aminor_1_tot = 0
                            nb_aminor_2_tot = 0
                            nb_non_aminor_1_tot = 0
                            nb_non_aminor_2_tot = 0
                            nb_aminor_tot = 0
                            nb_non_aminor_tot = 0
                            for liste in liste_non_redondant :
                                dico_trouves.update({(compteur, liste_expr, c) : [{},{}]})
                                dico_couples.update({(compteur, liste_expr, c) : {}})
                                for elt in liste :
                                        
                                        print(liste_seq_13[liste_expr])
                                        print(liste_seq_24[liste_expr])
                                        match_it_13 = re.finditer(liste_seq_13[liste_expr], dico_seq_tot[elt])
                                        match_it_24 = re.finditer(liste_seq_24[liste_expr], dico_seq_tot[elt])
                                        
                                        print(match_it_13)
                                        print(match_it_24)
                                        
                                        dico_trouves[(compteur, liste_expr, c)][0].update({elt : []})
                                        dico_trouves[(compteur, liste_expr, c)][1].update({elt : []})
                                        
                                        dico_couples[(compteur, liste_expr, c)].update({elt : []})
                                        
                                        for m in match_it_13 :
                                            pos = m.span()
                                            dico_trouves[(compteur, liste_expr, c)][0][elt].append(pos[0])
                                            print(pos)
                                            
                                        for m in match_it_24 :
                                            pos = m.span()
                                            dico_trouves[(compteur, liste_expr, c)][1][elt].append(pos[0])
                                            print(pos) 
                                        
                                        ## paires d'occurrences dans la même molécule
                                        for occ in dico_trouves[(compteur, liste_expr, c)][0][elt] :
                                            for occ2 in dico_trouves[(compteur, liste_expr, c)][1][elt] :
                                                dico_couples[(compteur, liste_expr, c)][elt].append([occ, occ2])
                                        
                                        dico_trouves[(compteur, liste_expr, c)][0][elt] = []
                                        dico_trouves[(compteur, liste_expr, c)][1][elt] = []
                                        
                                        ##recherche d'A-minor
                                        
                                        for aminor in liste_aminor :
                                            if elt[0] == aminor[0].upper() :
                                                     
                                                if elt[1] == aminor[1][0][0] and elt[1] == aminor[1][1][0]  :
                                                            for occ in dico_couples[(compteur, liste_expr, c)][elt] : 
                                                                print(occ)
                                                                if (occ[0] >= aminor[1][0][1] - 5 and occ[0] <= aminor[1][0][1]) or (occ[1] >= aminor[1][1][1] - 5 and occ[1] <= aminor[1][1][1]) :
                                                                    if (occ[0] >= aminor[1][0][1] - 5 and occ[0] <= aminor[1][0][1]) :
                                                                        if (occ[0], True) not in dico_trouves[(compteur, liste_expr, c)][0][elt] :
                                                                            dico_trouves[(compteur, liste_expr, c)][0][elt].append((occ[0], True))
                                                                    else :
                                                                        if (occ[0], False) not in dico_trouves[(compteur, liste_expr, c)][0][elt] :
                                                                            dico_trouves[(compteur, liste_expr, c)][0][elt].append((occ[0], False))
                                                                    if (occ[1] >= aminor[1][1][1] - 5 and occ[1] <= aminor[1][1][1]) :
                                                                        if (occ[1], True) not in dico_trouves[(compteur, liste_expr, c)][1][elt] :
                                                                            dico_trouves[(compteur, liste_expr, c)][1][elt].append((occ[1], True))
                                                                    else :
                                                                        if (occ[1], False) not in dico_trouves[(compteur, liste_expr, c)][1][elt] :
                                                                            dico_trouves[(compteur, liste_expr, c)][1][elt].append((occ[1], False))
                                                                    if (occ[0] >= aminor[1][0][1] - 5 and occ[0] <= aminor[1][0][1]) and (occ[1] >= aminor[1][1][1] - 5 and occ[1] <= aminor[1][1][1]) :
                                                                        occ.append(True)
                                                                    else :
                                                                        occ.append(False)
                                                                    
        
                                        
                                        
                                        
                                        
                                cle_max = -1
                                nb_tot_max = -1
                                nb_aminor_max = -1
                                
                                for cle in dico_couples[(compteur, liste_expr, c)] :
                                    nb_aminor = 0
                                    print(cle)
                                    for occ in dico_couples[(compteur, liste_expr, c)][cle] :
                                        if len(occ) == 3 and occ[2] == True :
                                            nb_aminor += 1 
                                    if nb_aminor > nb_aminor_max or (nb_aminor == nb_aminor_max and len(dico_couples[(compteur, liste_expr, c)][cle]) > nb_tot_max)  :
                                        nb_tot_max = len(dico_couples[(compteur, liste_expr, c)][cle])
                                        cle_max = cle
                                        nb_aminor_max = nb_aminor
                                
                                nb_aminor_1 = 0
                                nb_aminor_2 = 0
                                for occ in dico_trouves[(compteur, liste_expr, c)][0][cle_max] :
                                    if occ[1] == True :
                                        nb_aminor_1 += 1
                                        
                                for occ in dico_trouves[(compteur, liste_expr, c)][1][cle_max] :
                                    if occ[1] == True :
                                        nb_aminor_2 += 1
                                
                                nb_aminor_1_tot += nb_aminor_1
                                nb_aminor_2_tot += nb_aminor_2
                                nb_non_aminor_1_tot += len( dico_trouves[(compteur, liste_expr, c)][0][cle_max]) - nb_aminor_1
                                nb_non_aminor_2_tot += len( dico_trouves[(compteur, liste_expr, c)][1][cle_max]) - nb_aminor_2
                                nb_aminor_tot += nb_aminor_max
                                nb_non_aminor_tot += nb_tot_max - nb_aminor_max
                                    
                                    
                                
                                
                                c += 1
                            csvwriter.writerow([compteur, liste_seq_13[liste_expr], liste_seq_24[liste_expr], cluster, nb_non_aminor_1_tot, nb_non_aminor_2_tot, nb_aminor_1_tot, nb_aminor_2_tot, nb_non_aminor_tot, nb_aminor_tot])

#                 except :
#                     print("rate pas grave")              
                            
                        #exit()
            compteur += 1

def memes_voisins(graphe, noeud, graphe_commun, noeud_commun, aminor_edges, ch, couple, dico_couple):
    '''
        Donne la liste des aretes communes incidentes aux noeuds noeud/noeud_commun et indique si noeud_commun est impliqué dans une arete inter_brin
        
        :param graphe: un graphe de structure PDB
        :param noeud: un sommet de graphe dont on recherche s'il a des voisins communs avec noeud_commun
        :param graphe_commun: un sous-graphe commun maximum à un cluster
        :param noeud_commun: un sommet de graphe_commun 
        :param aminor_edges: la liste des aretes du motif A-minor (les noeuds du motif ont toujours les mêmes numéros dans tous les sous-graphes communs)
        :param ch: la liste des sommets du brin courant dans l'ordre de la séquence
        :param couple: la liste des aretes communes aux sommets de graphe et graphe_commun pour le brin courant
        :param dico_couple: le dictionnaire des couples de sommets de graphe et graphe_commun déjà associés
        
        :type graphe: nx.DiGraph
        :type noeud: sommet 
        :type graphe_commun: nx.MultiDiGraph
        :type noeud_commun: sommet
        :type aminor_edges: liste
        :type ch: liste
        :type couple: liste de tuples
        :type dico_couple: dict
        
        :return: la liste des aretes communes mise à jour et la liste des arêtes de graphe_commun qui sont inter brin sous forme (noeud_commun, voisin_commun, num de l'arete)
    '''
    
    compte_liaisons_non_b53_commun = 0
    compte_liaisons_non_b53 = 0
    inter_brin = -1
    
    ## juste pour les non apparies
    for voisin in graphe[noeud] : 
        if graphe.edges[noeud,voisin]["label"] != 'B53' :
            compte_liaisons_non_b53 += 1
            
    for voisin_commun in graphe_commun[noeud_commun] :
        for edge_commun in graphe_commun[noeud_commun][voisin_commun] :
            ## juste pour les non apparies
            if graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] != 'B53' : 
                compte_liaisons_non_b53_commun += 1
            if (noeud_commun, voisin_commun, graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"]) not in aminor_edges :
                if graphe_commun.nodes[voisin_commun]["type"] in [None,-1,15] :
                ## cas gnl (noeud de la chaine + noeud None ou artificiel ou 5)
                    for voisin in graphe[noeud] :
                                                
                            #if graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] == graphe.edges[noeud,voisin]["label"] and \

                        if (graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] == graphe.edges[noeud,voisin]["label"] or (graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] == "CWWn" and graphe.edges[noeud,voisin]["label"] == "CWW" and not ((graphe.nodes[noeud]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[noeud]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[noeud]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[noeud]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[noeud]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[noeud]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G')))) and \
                        graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] != "B53" and \
                        (noeud, noeud_commun, voisin, voisin_commun, graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"]) not in couple and \
                        (noeud_commun not in dico_couple or (noeud_commun in dico_couple and dico_couple[noeud_commun] == noeud)) and \
                        (voisin_commun not in dico_couple or (voisin_commun in dico_couple and dico_couple[voisin_commun] == voisin)) :
                            #print(noeud, noeud_commun, voisin, voisin_commun)
                            #print(graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"])
                            #print(graphe.edges[noeud,voisin]["label"])
                            #print((graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] == "CWWn" and graphe.edges[noeud,voisin]["label"] == "CWW" and not ((graphe.nodes[noeud]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[noeud]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[noeud]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[noeud]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[noeud]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[noeud]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G'))))
                            couple.append((noeud, noeud_commun, voisin, voisin_commun, graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"]))
                            if noeud_commun not in dico_couple : 
                                dico_couple.update({noeud_commun : noeud})
                            if voisin_commun not in dico_couple :
                                dico_couple.update({voisin_commun : voisin})
                else :
                   
                    #print("houhou")
                    ch1 = graphe_commun.nodes[noeud_commun]["chaine"]
                    ch2 = graphe_commun.nodes[voisin_commun]["chaine"]
                    
                    #inter = [value for value in ch1 if value in ch2] 
                    ## cas interchaine  (mais dans le même brin)
                    if (3 in ch1 and 3 in ch2) or (1 in ch1 and 1 in ch2) or (2 in ch1 and 2 in ch2) or (4 in ch1 and 4 in ch2) or ((3 in ch1 or 3 in ch2) and (1 in ch1 or 1 in ch2)) or ((2 in ch1 or 2 in ch2) and (4 in ch1 or 4 in ch2)) : ## l'arete relie deux sommets du même brin (soit au sein d'une des chaines soit entre chaines 1 et 3 ou chaines 2 et 4)
#                         print("ahaha")
#                         print(noeud_commun)
#                         print(voisin_commun)
#                         print(ch.index(noeud_commun))
#                         print(ch.index(voisin_commun))
#                         print(noeud)
#                         print(voisin)
                        #print("bouhhh")
                        for voisin in graphe[noeud] :
                            if (graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] == graphe.edges[noeud,voisin]["label"] or (graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] == "CWWn" and graphe.edges[noeud,voisin]["label"] == "CWW" and not ((graphe.nodes[noeud]["nt"] == 'A' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[noeud]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'A') or (graphe.nodes[noeud]["nt"] == 'C' and graphe.nodes[voisin]["nt"] == 'G') or (graphe.nodes[noeud]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'C') or (graphe.nodes[noeud]["nt"] == 'G' and graphe.nodes[voisin]["nt"] == 'U') or (graphe.nodes[noeud]["nt"] == 'U' and graphe.nodes[voisin]["nt"] == 'G')))) and \
                            graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] != "B53" and \
                            (noeud, noeud_commun, voisin, voisin_commun, graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"]) not in couple and \
                            abs(ch.index(noeud_commun) - ch.index(voisin_commun)) == abs(noeud-voisin) and \
                            (voisin, voisin_commun, noeud, noeud_commun, graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"][0]+ graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"][2]+graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"][1]) not in couple and \
                            (noeud_commun not in dico_couple or (noeud_commun in dico_couple and dico_couple[noeud_commun] == noeud)) and \
                            (voisin_commun not in dico_couple or (voisin_commun in dico_couple and dico_couple[voisin_commun] == voisin)) :
                                couple.append((noeud, noeud_commun, voisin, voisin_commun, graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"]))
                                if noeud_commun not in dico_couple : 
                                    dico_couple.update({noeud_commun : noeud})
                                if voisin_commun not in dico_couple :
                                    dico_couple.update({voisin_commun : voisin})
                    ## cas inter brin
                    elif graphe_commun[noeud_commun][voisin_commun][edge_commun]["label"] != "B53" :
                        #print("bahhh")
                        inter_brin = (noeud_commun, voisin_commun, edge_commun)
                        
                    
    ## si les noeuds sont non apparies on ajoute quand même une arete à la liste pour compter une arete artificielle
    if compte_liaisons_non_b53 == 0 and compte_liaisons_non_b53_commun == 0 :
        couple.append((noeud, noeud_commun, '0'))
        if noeud_commun not in dico_couple : 
            dico_couple.update({noeud_commun : noeud})
    return couple, inter_brin

def memes_voisins_inter(graphe, noeud_1, noeud_2, graphe_commun, noeud_commun_1, noeud_commun_2, edge_commun, aretes_communes, dico_couple) :
    '''
        Ajoute à la liste d'aretes communes l'arete entre noeud_1 et noeud_2 si elle est bien de même type que l'arete entre noeud_commun_1 et noeud_commun_2
        et si noeud_1 peut bien être associé à noeud_commun_1 et noeud_2 à noeud_commun_2
        
        :param graphe: le graphe de structure PDB contenant noeud_1 et noeud_2
        :param noeud_1: un sommet de graphe qu'on voudrait associer à noeud_commun_1
        :param noeud_2: un sommet de graphe qu'on voudrait associer à noeud_commun_2
        :param graphe_commun: le graphe commun maximum global contenant noeud_commun_1 et noeud_commun_2
        :param noeud_commun_1: un sommet de graphe_commun qu'on voudrait associer à noeud_1
        :param noeud_commun_2: un sommet de graphe_commun qu'on voudrait associer à noeud_2
        :param edge_commun: le numéro de l'arete entre noeud_commun_1 et noeud_commun_2 qui correspond peut-être à l'arete entre noeud_1 et noeud_2
        :param aretes_communes: la liste des aretes communes déjà trouvées entre les sommets de graphe et graphe_commun
        :param dico_couple: le dictionnaire contenant les couples de sommets de graphe et de graphe_commun déjà associés
        
        :type graphe: nx.DiGraph
        :type noeud_1: sommet
        :type noeud_2: sommet
        :type graphe_commun: nx.MultiDiGraph
        :type noeud_commun_1: sommet
        :type noeud_commun_2: sommet
        :type edge_commun: entier
        :type aretes_communes: liste
        :type dico_couple: dict
        
        :return: aretes_communes
    '''

#     print(noeud_1)
#     print(noeud_2)
#     print(graphe.nodes())
#     print(graphe[543])
#     print(graphe.edges[noeud_1,noeud_2]["label"])
    
    if (noeud_1, noeud_2) in graphe.edges and \
    (graphe_commun[noeud_commun_1][noeud_commun_2][edge_commun]["label"] == graphe.edges[noeud_1,noeud_2]["label"] or (graphe_commun[noeud_commun_1][noeud_commun_2][edge_commun]["label"] == "CWWn" and graphe.edges[noeud_1,noeud_2]["label"] == "CWW" and not ((graphe.nodes[noeud_1]["nt"] == 'A' and graphe.nodes[noeud_2]["nt"] == 'U') or (graphe.nodes[noeud_1]["nt"] == 'U' and graphe.nodes[noeud_2]["nt"] == 'A') or (graphe.nodes[noeud_1]["nt"] == 'C' and graphe.nodes[noeud_2]["nt"] == 'G') or (graphe.nodes[noeud_1]["nt"] == 'G' and graphe.nodes[noeud_2]["nt"] == 'C') or (graphe.nodes[noeud_1]["nt"] == 'G' and graphe.nodes[noeud_2]["nt"] == 'U') or (graphe.nodes[noeud_1]["nt"] == 'U' and graphe.nodes[noeud_2]["nt"] == 'G')))) and \
    (noeud_commun_1 not in dico_couple or (noeud_commun_1 in dico_couple and dico_couple[noeud_commun_1] == noeud_1)) and \
    (noeud_commun_2 not in dico_couple or (noeud_commun_2 in dico_couple and dico_couple[noeud_commun_2] == noeud_2)) and \
    (noeud_1, noeud_commun_1, noeud_2, noeud_commun_2, graphe.edges[noeud_1,noeud_2]["label"]) not in aretes_communes and \
    (noeud_2, noeud_commun_2, noeud_1, noeud_commun_1, graphe.edges[noeud_1,noeud_2]["label"][0]+graphe.edges[noeud_1,noeud_2]["label"][2]+graphe.edges[noeud_1,noeud_2]["label"][1]) not in aretes_communes :
        aretes_communes.append((noeud_1, noeud_commun_1, noeud_2, noeud_commun_2, graphe.edges[noeud_1,noeud_2]["label"]))
    
    return aretes_communes

def tot_arete_non_cov(graphe_commun, aminor_edges):
    ''' 
        Donne le nombre d'aretes du graphe_commun n'appartenant pas au motif A-minor (pondéré par le poids des sommets) 
        
        :param graphe_commun: un sous-graphe commun maximum à un cluster
        :param aminor_edges: la liste des aretes du motif A-minor 
        
        :type graphe_commun: nx.MultiDiGraph
        :type aminor_edges: liste
        
        :return: le nombre d'aretes total du graphe_commun sans les aretes du motif, le nombre d'aretes du brin 1, le nombre d'aretes du brin 2, le nomrbe d'aretes inter brin 
    '''
    nb_aretes = 0
    nb_aretes_ch1 = 0
    nb_aretes_ch2 = 0
    nb_aretes_inter = 0
    for u,v,data in graphe_commun.edges(data=True) :
        if data["label"] != 'B53' and (u,v,data["label"]) not in aminor_edges :
            #if data["near"] != True :
                nb_aretes += graphe_commun.nodes[u]["poids"]
                print(u,v,data)
                ch1 = graphe_commun.nodes[u]["chaine"]
                ch2 = graphe_commun.nodes[v]["chaine"]
                if graphe_commun.nodes[u]["type"] not in [None, -1, 15] and graphe_commun.nodes[v]["type"] not in [None, -1, 15] and not ((3 in ch1 and 3 in ch2) or (1 in ch1 and 1 in ch2) or (2 in ch1 and 2 in ch2) or (4 in ch1 and 4 in ch2) or ((3 in ch1 or 3 in ch2) and (1 in ch1 or 1 in ch2)) or ((2 in ch1 or 2 in ch2) and (4 in ch1 or 4 in ch2))) : ## l'arete relie deux sommets de chaine différente
                    nb_aretes_inter += graphe_commun.nodes[u]["poids"]
                if (graphe_commun.nodes[u]["type"] not in [None, -1, 15] and (1 in ch1 or 3 in ch1)) or  (graphe_commun.nodes[v]["type"] not in [None, -1, 15] and (1 in ch2 or 3 in ch2)):
                    nb_aretes_ch1 += graphe_commun.nodes[u]["poids"]
                else :
                    nb_aretes_ch2 += graphe_commun.nodes[u]["poids"]
    print(nb_aretes)
    return nb_aretes/2, nb_aretes_ch1/2, nb_aretes_ch2/2, nb_aretes_inter/2

def recherche_struct(clusters, liste_non_redondant, graphes, liste_aminor):
    '''
        Algorithme de recherche de la sous-structure max des clusters dans les structures PDB de liste_non_redondant
        
        :param clusters: la liste des clusters
        :param liste_non_redondant: la liste des structures PDB à considérer
        :param graphes: un dictionaire contenant les toutes les structures PDB sous forme de graphes nx
        :param liste_aminor: la liste des positions des motifs A-minor dans toutes les structures de la PDB
        
        :type clusters: liste de listes
        :type liste_non_redondant: liste
        :type graphes: dictionnaire de nx.DiGraph
        :type liste_aminor : liste
    '''
    
    with open("Resultats_sequences/res_csv_struct_2_0.75_suite.csv", 'a', newline="") as fichier_csv :
                csvwriter = csv.writer(fichier_csv)
                csvwriter.writerow(["Structure", "Position chaine 1", "Position chaine 2", "Num cluster", "Cluster", "Positions ch1", "Positions ch2", "Nombre d'aretes trouvees", "Nombre d'aretes trouves dans la chaine1", "Nombre d'aretes trouves dans la chaine 2", "Nombre d'aretes au total"])
        
#         with open("Resultats_sequences/res_csv_struct_aminor_2_0.75_suite.csv", 'a', newline="") as fichier_csv_aminor :
#             csvwriter_aminor = csv.writer(fichier_csv_aminor)
#             csvwriter_aminor.writerow(["Structure", "Position chaine 1", "Position chaine 2", "Num cluster", "Cluster", "Positions ch1", "Positions ch2", "Nombre d'aretes trouvees", "Nombre d'aretes trouves dans la chaine1", "Nombre d'aretes trouves dans la chaine 2", "Nombre d'aretes au total"])
# 
#             with open("Resultats_sequences/res_csv_struct_semi_aminor_2_0.75_suite.csv", 'a', newline="") as fichier_csv_aminor :
#                 csvwriter_semi_aminor = csv.writer(fichier_csv_aminor)
#                 csvwriter_semi_aminor.writerow(["Structure", "Position chaine 1", "Position chaine 2", "Num cluster", "Cluster", "Positions ch1", "Positions ch2", "Nombre d'aretes trouvees", "Nombre d'aretes trouves dans la chaine1", "Nombre d'aretes trouves dans la chaine 2", "Nombre d'aretes au total"])
                
                
#                 with open("Resultats_sequences/dico_recherche_struct.pickle", 'rb') as fichier_pickle_1 : 
#                     mon_depickler = pickle.Unpickler(fichier_pickle_1)
#                     dico_recherche_struct = mon_depickler.load()
#                     
#                 with open("Resultats_sequences/dico_recherche_struct_aminor.pickle", 'rb') as fichier_pickle_2 : 
#                     mon_depickler = pickle.Unpickler(fichier_pickle_2)
#                     dico_recherche_struct_aminor = mon_depickler.load()
#                     
#                 with open("Resultats_sequences/dico_recherche_struct_semi_aminor.pickle", 'rb') as fichier_pickle_3 : 
#                     mon_depickler = pickle.Unpickler(fichier_pickle_3)
#                     dico_recherche_struct_semi_aminor = mon_depickler.load()
                
                dico_recherche_struct = {}
                dico_recherche_struct_aminor = {}
                dico_recherche_struct_semi_aminor = {}
                compteur = 1
                liste = []
                liste_pas_bon = []
                liste_trop_haut = []
                liste_chaine_commune = []
                for cluster in clusters :
                    print(compteur)
                    if len(cluster) > 2  :
                        print(len(cluster))
                        print(cluster)
                        
                        nom = cluster[0]
                        #exit()
                        ## calcule le graphe commun du cluster
                        graphe_commun_global, dico_graphes, graphe_commun = recherche_graphe_commun([4], cluster, 0.6, "ramou", "Graphes_communs_mai_2020/extensions", compteur)
                        
                        #nom = (cluster[0][0], graphe_commun_global.nodes[1]["num_ch"])
                        print(nom)
                        print(graphe_commun_global.nodes[1]["num_seq"][0])
                        
                        ##calcule la position d'un A-minor (pour tester au début)
                        pos1 = graphe_commun_global.nodes[1]["num_seq"][0][1][0] - 4
                        print(pos1)
                        pos2 = graphe_commun_global.nodes[2]["num_seq"][0][1][0] - 3
                        print(pos2)
                        #exit()
                        chaines_1 = []
                        chaines_2 = []
                        chaine_1 = [-1]*8
                        chaine_2 = [-1]*8
                        chaines_1.append(chaine_1)
                        chaines_2.append(chaine_2)
                        
                        for noeud, data in graphe_commun_global.nodes(data=True) :
                            print(noeud, data)
                            
                        for noeud, data in graphe_commun.nodes(data=True) :
                            print(noeud, data)
                        
                        ## établissement des listes de sommets
                        lien_B53 = False
                        for noeud, data in graphe_commun.nodes(data=True) :
                            if len(data["position"]) == 1 :
                                print("est on coince")
                                print(chaines_1)
                                print(chaines_2)
                                print(noeud)
                                if 1 in data["chaine"] :
                                    if min(data["espacement_motif"]) == max(data["espacement_motif"]) :
                                        for chaine in chaines_1 :
                                            for i in range(data["poids"]) :
                                                    j = 0
                                                    noeud_temp = noeud
                                                    while j < i :
                                                        voisin_B53 = -1
                                                        for voisin in graphe_commun_global[noeud_temp] :
                                                            for edge in graphe_commun_global[noeud_temp][voisin] :
                                                                if graphe_commun_global[noeud_temp][voisin][edge]["label"] == "B53" :
                                                                    voisin_B53 = voisin
                                                        j += 1
                                                        noeud_temp = voisin_B53
                                                    if chaine[data["position"][0]+4+i] == -1 :
                                                        chaine[data["position"][0]+4+i] = noeud_temp
                                            
                                    else :
                                        chaines_1_temp = list(chaines_1)  
                                        for chaine in chaines_1_temp :
                                            for elt in set(data["espacement_motif"]) :
                                                new_chaines_1 = list(chaine)
                                                for i in range(data["poids"]) :
                                                    j = 0
                                                    noeud_temp = noeud
                                                    while j < i :
                                                        voisin_B53 = -1
                                                        for voisin in graphe_commun_global[noeud_temp] :
                                                            for edge in graphe_commun_global[noeud_temp][voisin] :
                                                                if graphe_commun_global[noeud_temp][voisin][edge]["label"] == "B53" :
                                                                    voisin_B53 = voisin
                                                        j += 1
                                                        noeud_temp = voisin_B53
                                                    if new_chaines_1[elt+4+i] == -1 :
                                                        new_chaines_1[elt+4+i] = noeud_temp
                                                    
                                                chaines_1.append(new_chaines_1)
                                            chaines_1.remove(chaine)
                            
                                if 3 in data["chaine"] :
                                    if min(data["espacement_motif"]) == max(data["espacement_motif"]) :
                                        for chaine in chaines_1 :
                                            for i in range(data["poids"]) :
                                                    j = 0
                                                    noeud_temp = noeud
                                                    while j < i :
                                                        voisin_B53 = -1
                                                        for voisin in graphe_commun_global.predecessors(noeud_temp) :
                                                            for edge in graphe_commun_global[voisin][noeud_temp] :
                                                                if graphe_commun_global[voisin][noeud_temp][edge]["label"] == "B53" :
                                                                    voisin_B53 = voisin
                                                        j += 1
                                                        noeud_temp = voisin_B53
                                                    if chaine[data["position"][0]-7-i] == -1 :
                                                        chaine[data["position"][0]-7-i] = noeud_temp
                                    else :
                                        chaines_1_temp = list(chaines_1)  
                                        for chaine in chaines_1_temp :
                                            for elt in set(data["espacement_motif"]) :
                                                new_chaines_1 = list(chaine)
                                                for i in range(data["poids"]) :
                                                    j = 0
                                                    noeud_temp = noeud
                                                    while j < i :
                                                        voisin_B53 = -2
                                                        for voisin in graphe_commun_global.predecessors(noeud_temp) :
                                                            for edge in graphe_commun_global[voisin][noeud_temp] :
                                                                if graphe_commun_global[voisin][noeud_temp][edge]["label"] == "B53" :
                                                                    voisin_B53 = voisin
                                                        j += 1
                                                        noeud_temp = voisin_B53
            
                                                    if new_chaines_1[3-elt-i] == -1 :
                                                        new_chaines_1[3-elt-i] = noeud_temp
                                                    
                                                chaines_1.append(new_chaines_1)
                                            chaines_1.remove(chaine)
                                
                                if 2 in data["chaine"] :
                                    if min(data["espacement_motif"]) == max(data["espacement_motif"]) :
                                        for chaine in chaines_2 :
                                            for i in range(data["poids"]) :
                                                    j = 0
                                                    noeud_temp = noeud
                                                    while j < i :
                                                        voisin_B53 = -1
                                                        for voisin in graphe_commun_global.predecessors(noeud_temp) :
                                                            for edge in graphe_commun_global[voisin][noeud_temp] :
                                                                if graphe_commun_global[voisin][noeud_temp][edge]["label"] == "B53" :
                                                                    voisin_B53 = voisin
                                                        j += 1
                                                        noeud_temp = voisin_B53
                                                    print(noeud_temp)
                                                    if chaine[data["position"][0]-7-i] == -1 :
                                                        chaine[data["position"][0]-7-i] = noeud_temp
                                    else : 
                                        chaines_2_temp = list(chaines_2)  
                                        for chaine in chaines_2_temp :
                                            for elt in set(data["espacement_motif"]) :
                                                    new_chaines_2 = list(chaine)
                                                    for i in range(data["poids"]) :
                                                        j = 0
                                                        noeud_temp = noeud
                                                        while j < i :
                                                            voisin_B53 = -2
                                                            for voisin in graphe_commun_global.predecessors(noeud_temp) :
                                                                for edge in graphe_commun_global[voisin][noeud_temp] :
                                                                    if graphe_commun_global[voisin][noeud_temp][edge]["label"] == "B53" :
                                                                        voisin_B53 = voisin
                                                            j += 1
                                                            noeud_temp = voisin_B53
                                                        if new_chaines_2[3-elt-i] == -1 :
                                                            new_chaines_2[3-elt-i] = noeud_temp
                                                            
                                                    chaines_2.append(new_chaines_2)
                                                    print(chaine)
                                            chaines_2.remove(chaine)
                                                
            
                                        print(chaines_2)
                                        #exit()
            
                                            
                                    
                                if 4 in data["chaine"] :
                                    if min(data["espacement_motif"]) == max(data["espacement_motif"]) :
                                        for chaine in chaines_2 :
                                            for i in range(data["poids"]) :
                                                    j = 0
                                                    noeud_temp = noeud
                                                    while j < i :
                                                        voisin_B53 = -1
                                                        for voisin in graphe_commun_global[noeud_temp] :
                                                            for edge in graphe_commun_global[noeud_temp][voisin] :
                                                                if graphe_commun_global[noeud_temp][voisin][edge]["label"] == "B53" :
                                                                    voisin_B53 = voisin
                                                        j += 1
                                                        noeud_temp = voisin_B53
                                                    if chaine[data["position"][0]+4+i] == -1 :
                                                        chaine[data["position"][0]+4+i] = noeud_temp
                                                        print(noeud)
                                                        print(data["position"][0]+4+i)
                                    
                                    else :   
                                        chaines_2_temp = list(chaines_2)  
                                        for chaine in chaines_2_temp :
                                            for elt in set(data["espacement_motif"]) :
                                                new_chaines_2 = list(chaine)
                                                for i in range(data["poids"]) :
                                                    j = 0
                                                    noeud_temp = noeud
                                                    while j < i :
                                                        voisin_B53 = -1
                                                        for voisin in graphe_commun_global[noeud_temp] :
                                                            for edge in graphe_commun_global[noeud_temp][voisin] :
                                                                if graphe_commun_global[noeud_temp][voisin][edge]["label"] == "B53" :
                                                                    voisin_B53 = voisin
                                                        j += 1
                                                        noeud_temp = voisin_B53
                                                        
                                                    if new_chaines_2[elt+4+i] == -1 :
                                                        new_chaines_2[elt+4+i] = noeud_temp
                                                    
                                                chaines_2.append(new_chaines_2)
                                            chaines_2.remove(chaine)
                            elif len(data["position"]) > 1 :
                                lien_B53 = True
            #             c = 0
            #             for c in range(len(chaines_1)) :
            #                 if chaines_1[c] != chaines_1[c+1]
                        if lien_B53 : 
                            liste_chaine_commune.append(compteur)
                        else :
                            aminor_edges = [(1,2,"CSS"),(2,1,"CSS"),(3,1, "B53"),(1,5, "TSS"),(5,1, "TSS"),(2,4, "B53"),(3,4,"CSS"), (4,3, "CSS"), (2,5,"CWW"),(5,2,"CWW")]
                            
                            
                            
                            
                            
                            print(chaines_1)
                            print(chaines_2)
                            print(lien_B53)
                            
                            ## calcule le nombre d'aretes du graphe commun maximum et le nombre d'aretes par brin du sous-graphe commun maximum
                            nb_aretes, nb_aretes_ch1, nb_aretes_ch2, nb_aretes_inter = tot_arete_non_cov(graphe_commun, aminor_edges)
                            print(nb_aretes, nb_aretes_ch1, nb_aretes_ch2, nb_aretes_inter)
                            
                            ## indique si on doit considérer le brin 2 avant le brin 1 en boucle extérieure (on le fait s'il y a plus d'aretes dans le brin 2 que dans le brin 1)
                            change_sens = False
                            if nb_aretes_ch1 < nb_aretes_ch2 :
                                change_sens = True
                            
            #                 for u,v,data in graphe_commun_global.edges(data=True) :
            #                     print(u,v,data)
                            #exit()
                            ok = False
                            compteur_elt = 0
                            for elt in liste_non_redondant :
                                if graphes[elt].number_of_nodes() >= 16 : 
                                
                                    
#                                         liste_pos_aminor_1 = []
#                                         liste_pos_aminor_2 = []
#                                         
#                                         for aminor in liste_aminor :
#                                             if elt[0] == aminor[0].upper() and elt[1] == aminor[1][0][0] :
#                                                 liste_pos_aminor_1.append(aminor[1][0][1] - 3)
#                                                 liste_pos_aminor_2.append(aminor[1][1][1] - 3)
#                                         
#                                         print(liste_pos_aminor_1)
#                                         print(liste_pos_aminor_2)
                                        
                                        dico_struct_1 = {}
                                        dico_struct_2 = {}
                                        
                                        if change_sens :
                                            chaine_new_1 = chaines_2
                                            chaine_new_2 = chaines_1
                                            nb_aretes_ch1_new = nb_aretes_ch2
                                            nb_aretes_ch2_new = nb_aretes_ch1
                                            #liste_pos_aminor_1_new = liste_pos_aminor_2
                                            #liste_pos_aminor_2_new = liste_pos_aminor_1
                                        else :
                                            chaine_new_1 = chaines_1
                                            chaine_new_2 = chaines_2
                                            nb_aretes_ch1_new = nb_aretes_ch1
                                            nb_aretes_ch2_new = nb_aretes_ch2
                                            #liste_pos_aminor_1_new = liste_pos_aminor_1
                                            #liste_pos_aminor_2_new = liste_pos_aminor_2
                                        
                                        for i in range(1, graphes[elt].number_of_nodes()-8) :
                                            #if compteur == 4  and elt == ('1VQP', '0'):
                                            #if i == 474 or i == 693 :
                                                c1 = 0
                                                for ch1 in chaine_new_1 :
                                                    print(ch1)
                                                #print(graphes[elt].nodes())
                                                
                        #                         for noeud  in graphes[elt].nodes() :
                        #                             print(noeud, graphes[elt][noeud])
                                                
                                
                                                    print("i=%d"%i)
                                                    noeuds_inter_brin_1 = []
                                                    aretes_communes_1 = []
                                                    dico_couple_1 = {}
                                                    for j in range(len(ch1)) :
                                                            if ch1[j] != -1 :
                                                                aretes_communes_1, inter_brin = memes_voisins(graphes[elt], i+j, graphe_commun_global, ch1[j], aminor_edges, ch1, aretes_communes_1, dico_couple_1)
                                                                if inter_brin != -1 : 
                                                                    noeuds_inter_brin_1.append(inter_brin)
                                                    dico_struct_1.update({(i, c1) : (dico_couple_1, aretes_communes_1, noeuds_inter_brin_1)})
                                                    c1 += 1
                                                
                                                c2 = 0
                                                for ch2 in chaine_new_2 :
                                                    print(ch2)
        
                                                    aretes_communes_2 = []
                                                    noeuds_inter_brin_2 = []
                                                    dico_couple_2 = {}
                                                    for m in range(len(ch2)) :
                                                        if ch2[m] != -1 :
                                                            aretes_communes_2, inter_brin = memes_voisins(graphes[elt], i+m, graphe_commun_global, ch2[m], aminor_edges, ch2, aretes_communes_2, dico_couple_2)
                                                            if inter_brin != -1 :
                                                                noeuds_inter_brin_2.append(inter_brin)
                                                    dico_struct_2.update({(i, c2) : (dico_couple_2, aretes_communes_2, noeuds_inter_brin_2)})
                                                    c2 += 1
                                        
                                        #print(dico_struct_1)
                                        #print(dico_struct_2)
                                        #return
                                        
                                        #exit()
                            #elt = (nom[0].upper(), nom[1].upper())
                            #print(elt)
                                        print(elt)
                                    #if elt[0] == nom[0].upper() and elt[1] == nom[1].upper() : 
                                        print(graphes[elt].number_of_nodes())
                                            
                                        
                                        
                                        c1 = 0
                                        for ch1 in chaine_new_1 :
                                            print(ch1)
                                            #print(graphes[elt].nodes())
                                            
                    #                         for noeud  in graphes[elt].nodes() :
                    #                             print(noeud, graphes[elt][noeud])
                                            
                                            for i in range(1, graphes[elt].number_of_nodes()-8) :
                                                print("i=%d"%i)
                                                print(0.75*nb_aretes - nb_aretes_ch2_new)
                                                ## on ne teste le i que s'il est possible d'avoir une proportion totale supérieure à 0.75
                                                if len(dico_struct_1[(i, c1)][1]) >= 1.0*nb_aretes - nb_aretes_ch2_new - nb_aretes_inter :
                                                    c2 = 0
                                                    for ch2 in chaine_new_2 :
                                                            print(ch2)
                                                        
                                                            for k in range(1,graphes[elt].number_of_nodes()-8) :
                                                                    #print("i=%d"%i)
                                                                #if k == 1 and c2 == 1 :
                                                    
                                                             
                                                                    print(elt)
                                                                    print(compteur)
                                                                    print("i=%d"%i)
                                                                    print("k=%d"%k)
                                                                    #print(i in liste_pos_aminor_1)
                                                                    #print(k in liste_pos_aminor_2)
                                                                    print(chaines_1)
                                                                    
                                                                #if k == pos2 :
                                                                    ## on ne teste le k que s'il est possible d'avoir une proportion totale supérieure à 0.75
                                                                    if (k > i+8 or k < i-8) and len(dico_struct_2[(k, c2)][1]) >= 1.0*nb_aretes - nb_aretes_ch1_new - nb_aretes_inter :
                                                                        print(0.75*nb_aretes - nb_aretes_ch1_new)
                                                                        
                                                                        aretes_communes = list(dico_struct_1[(i, c1)][1])
                                                                        aretes_communes.extend(dico_struct_2[(k, c2)][1])
                                                                        aretes_communes_1 = list(dico_struct_1[(i, c1)][1])
                                                                        aretes_communes_2 = list(dico_struct_2[(k, c2)][1])
                                                                        
                                                                        dico_couple = dico_struct_1[(i,c1)][0].copy()
                                                                        dico_couple.update(dico_struct_2[(k,c2)][0].copy())
                                                                        
                                                                        noeuds_inter_brin_1 = list(dico_struct_1[(i,c1)][2])
                                                                        noeuds_inter_brin_2 = list(dico_struct_2[(k,c2)][2])
                                                                        #print(dico_couple)
                                                                        #print(noeuds_inter_brin_1)
                                                                        #print(noeuds_inter_brin_2)
                                                                        ## Cas des aretes inter_brin à prendre en compte
                                                                        if len(noeuds_inter_brin_1) > 0 and len(noeuds_inter_brin_2) > 0 :
                                                                                for triplet in noeuds_inter_brin_1 :
                                                                                    print("on est laaaa")
                                                                                    print(ch1.index(triplet[0]))
                                                                                    print(ch2.index(triplet[1]))
                                                                                    aretes_communes = memes_voisins_inter(graphes[elt], i+ch1.index(triplet[0]), k+ch2.index(triplet[1]), graphe_commun_global, triplet[0], triplet[1], triplet[2], aretes_communes, dico_couple)
                                                                        
                                                                        #print(aretes_communes)
                                                                        nb_tot = len(aretes_communes)
                                                                        if i == 693 and k == 1 :
                                                                            print(nb_tot)
                                                                            #return
                                                                        print(nb_tot)
                                                                        if nb_tot >= 1.0*nb_aretes : 
#                                                                             if i in liste_pos_aminor_1_new and k in liste_pos_aminor_2_new :
#                                                                                 #print("est_la")
#                                                                                 #return
#                                                                                 csvwriter_aminor.writerow([elt, i, k, compteur, cluster, ch1, ch2, nb_tot, len(aretes_communes_1), len(aretes_communes_2), nb_aretes, change_sens])
#                                                                                 dico_recherche_struct_aminor.update({(elt, i, k, compteur, str(ch1), str(ch2)) : (aretes_communes, dico_couple)})
#                                                                             elif i in liste_pos_aminor_1_new or k in liste_pos_aminor_2_new :
#                                                                                 csvwriter_semi_aminor.writerow([elt, i, k, compteur, cluster, ch1, ch2, nb_tot, len(aretes_communes_1), len(aretes_communes_2), nb_aretes, change_sens])
#                                                                                 dico_recherche_struct_semi_aminor.update({(elt, i, k, compteur, str(ch1), str(ch2)) : (aretes_communes, dico_couple)})
#     
#                                                                             else :
                                                                                csvwriter.writerow([elt, i, k, compteur, cluster, ch1, ch2, nb_tot, len(aretes_communes_1), len(aretes_communes_2), nb_aretes, change_sens])
                                                                                dico_recherche_struct.update({(elt, i, k, compteur, str(ch1), str(ch2)) : (aretes_communes, dico_couple)})
                    #                                                     if nb_tot == nb_aretes : 
                    #                                                         ok = True
        #                                                             if i == 474 and k == 693 : 
        #                                                                 print(nb_tot)
        #                                                                 print(nb_aretes)
        #                                                                 print("rapoulou")
        #                                                                 return
        
        
                                                            c2 += 1
                                            c1 += 1        
                                compteur_elt += 1
                                        #return
                                            #exit()           
            #                 if not ok :
            #                     if nb_tot < nb_aretes :
            #                         liste_pas_bon.append(compteur)
            #                     else :
            #                         liste_trop_haut.append(compteur)
            #                     print("diff")
            #                     print(elt)
            #                     print(cluster)
            #                     print(compteur)
            #                     print(nb_tot)
            #                     print(nb_aretes)
            #                     #exit()     
            #                 else :
            #                     liste.append(compteur)           
                                
                            
                            
                        
                        #exit()
                        
                        #break
                    compteur += 1   
                 
                
    with open("Resultats_sequences/dico_recherche_struct.pickle", 'wb') as fichier_pickle_1 : 
        mon_pickler = pickle.Pickler(fichier_pickle_1)
        mon_pickler.dump(dico_recherche_struct)
        
    with open("Resultats_sequences/dico_recherche_struct_aminor.pickle", 'wb') as fichier_pickle_2 : 
        mon_pickler = pickle.Pickler(fichier_pickle_2)
        mon_pickler.dump(dico_recherche_struct_aminor)
        
    with open("Resultats_sequences/dico_recherche_struct_semi_aminor.pickle", 'wb') as fichier_pickle_3 : 
        mon_pickler = pickle.Pickler(fichier_pickle_3)
        mon_pickler.dump(dico_recherche_struct_semi_aminor)
#     print(liste)    
#     print(liste_pas_bon)  
#     print(liste_trop_haut)
#     print(liste_chaine_commune)     

def analyse_resultats_struct(liste_aminor, liste_nom_graphes_aminor):
    '''
        Creation de dictionnaires pour stocker les résultats de la recherche des sous-structures cluster par cluster
        
        :param liste_aminor: liste des positions de motifs A-minor dans les structures PDB
        :param liste_nom_graphes_aminor: liste des noms de structures non redondantes de notre jeu de données
        
        :type liste_aminor: liste 
        :type liste_nom_graphes_aminor: liste
        
        :return: stocke dans des fichiers pickle les dictionnaires 
    '''
    
    
    with open("Resultats_sequences/res_csv_struct_2_0.75_suite.csv", 'r', newline="") as fichier_csv :
        csvreader = csv.reader(fichier_csv)
    
        
        dico = {}
        
        compteur = 0
        new_row_3 = -1
        for row in csvreader :
            #if row[0] != "Structure" and (row[0].split(",")[0][2:len(row[0].split(",")[0])-1], row[0].split(",")[1][2:len(row[0].split(",")[1])-2]) in liste_nom_graphes_aminor :
            if row[0] != "Structure" :
                ## ici on ne récupère que les occurrences qui ont une proportion de 1
                if int(row[7])/float(row[10]) == 1.0 :
                    print(row[3])
                    if row[3] not in dico :
                        if new_row_3 != -1 :
                            if "dico_cluster_%s_que_1_2.pickle"%new_row_3 in os.listdir("Resultats_sequences") :
                                with open("recup_data/Resultats_sequences/dico_cluster_%s_que_1_2.pickle"%new_row_3, 'rb') as fichier_pickle :
                                    mon_depickler = pickle.Unpickler(fichier_pickle)
                                    dico_cluster = mon_depickler.load()
                                    
                                    dico[new_row_3].extend(dico_cluster)
                                
                            with open("Resultats_sequences/dico_cluster_%s_que_1_2.pickle"%new_row_3, 'wb') as fichier_pickle :
                                mon_pickler = pickle.Pickler(fichier_pickle)
                                mon_pickler.dump(dico[new_row_3])
                        dico.update({row[3] : []})
                        new_row_3 = row[3]
                        
                    vrai_aminor = False
                    aminor_1 = False
                    aminor_2 = False
                    ## dans le cas où on peut avoir des A-minor (struct A-minor), on regarde si les deux positions appartiennent au même A-minor
                    for cle in liste_aminor :
                        #print(row[0].split(",")[0][2:len(row[0].split(",")[0])-1])
                        #print(row[0].split(",")[1][2:len(row[0].split(",")[1])-2])
                        #exit()
                        
                        if row[0].split(",")[0][2:len(row[0].split(",")[0])-1] == cle.upper() :
                            #print(dico_aminor[cle])
                            #print(cle)
                            #print(row[0])
                            for elt in liste_aminor[cle] :
                                #print("ahaa")
                                #print(row[1])
                                #print(elt[1][0][1] - 3)
                                #print(elt)
                                #print(elt[0][0])
                                #print(elt[0][1])
                                #print(elt[1][1])
                                if row[11] == 'False' :
                                    if row[0].split(",")[1][2:len(row[0].split(",")[1])-2] == elt[0][0] :
                                        if int(row[1]) == elt[0][1] - 3 or int(row[2]) == elt[1][1] -3  :
                                            if int(row[1]) == elt[0][1] - 3 and int(row[2]) == elt[1][1] -3 :
                                                vrai_aminor = True
                                                aminor_1 = True
                                                aminor_2 = True
                                                break
                                            elif int(row[1]) == elt[0][1] - 3 :
                                                aminor_1 = True
                                            elif int(row[2]) == elt[1][1] -3 :
                                                aminor_2 = True
                                            #break
                                else :
                                    if row[0].split(",")[1][2:len(row[0].split(",")[1])-2] == elt[0][0] :
                                        if int(row[1]) == elt[1][1] - 3 or int(row[2]) == elt[0][1] -3 :
                                            if int(row[1]) == elt[1][1] - 3 and int(row[2]) == elt[0][1] -3 :
                                                vrai_aminor = True
                                                aminor_1 = True
                                                aminor_2 = True
                                                break
                                            elif int(row[1]) == elt[1][1] - 3 :
                                                aminor_1 = True
                                            elif int(row[2]) == elt[0][1] -3 :
                                                aminor_2 = True
                                        #vrai_aminor = True
                                        #break
                                if aminor_1 and aminor_2 :
                                    break
                                #exit()
                    ## si on a échangé les brins 1 et 2
                    if row[11] == 'False' :       
                        item = (row[0], row[1], row[2], int(row[7])/float(row[10]), vrai_aminor, aminor_1, aminor_2)
                    else :
                        item = (row[0], row[2], row[1], int(row[7])/float(row[10]), vrai_aminor, aminor_2, aminor_1)
                    
                    ## si on a stocké plusieurs fois la même occurrence
                    existe_deja = False
                    for elt in dico[row[3]] :
                        if elt[0] == item[0] and elt[1] == item[1] and elt[2] == item[2] and elt[3] == item[3] :
                            existe_deja = True
                            break
                    if not existe_deja :
                        dico[row[3]].append(item)
                
            compteur += 1
            print(compteur)
            #if compteur == 100 :
            #    print(dico)
            #    break
        #print(dico['3'])
        
        if "dico_cluster_%s_que_1_2.pickle"%new_row_3 in os.listdir("Resultats_sequences") :
                with open("recup_data/Resultats_sequences/dico_cluster_%s_que_1_2.pickle"%new_row_3, 'rb') as fichier_pickle :
                    mon_depickler = pickle.Unpickler(fichier_pickle)
                    dico_cluster = mon_depickler.load()
                    
                    dico[new_row_3].extend(dico_cluster)
                
        with open("Resultats_sequences/dico_cluster_%s_que_1_2.pickle"%new_row_3, 'wb') as fichier_pickle :
                mon_pickler = pickle.Pickler(fichier_pickle)
                mon_pickler.dump(dico[new_row_3])
        


def analyse_nombre_que_1_par_cluster(num_cluster, csvwriter, cluster):
    '''
        Création d'un fichier csv pour les résultats statistiques cluster par cluster (nombre de structures, moyenne du nombre d'occurrence par structure)
        
        :param num_cluster: le numéro du cluster qu'on considère
        :param csvwriter: le flux du fichier csv dans lequel on veut écrire
        :param cluster: la liste des éléments du cluster
        
        :type num_cluster: entier
        :type csvwriter: un object writer
        :type cluster: liste
        
        :return: crée un fichier csv avec les informations (un cluster par ligne)
    '''
    with open("recup_data/Resultats_sequences/dico_cluster_%s_que_1.pickle"%num_cluster, 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        dico_cluster = mon_depickler.load()
        
    if "dico_cluster_%s.pickle"%num_cluster in os.listdir("recup_data/Resultats_sequences/") :
        with open("recup_data/Resultats_sequences/dico_cluster_%s.pickle"%num_cluster, 'rb') as fichier_pickle_2 :
            mon_depickler = pickle.Unpickler(fichier_pickle_2)
            dico_cluster_2 = mon_depickler.load()
        
            dico_cluster.extend(dico_cluster_2)
    
        
    dico_nb_struct = {}
    dico_nb_struct_1_aminor = {}
    dico_nb_struct_0_aminor = {}
    liste_hors_aminor = []
    liste_hors_tot_aminor = []
    for elt in dico_cluster :
        if elt[3] == 1.0 :
            print(elt)
            if elt[0] not in dico_nb_struct :
                dico_nb_struct.update({elt[0] : 1})
            else :
                dico_nb_struct[elt[0]] += 1
            if not elt[4] and not elt[5] and not elt[6] :
                if elt[0] not in dico_nb_struct_0_aminor :
                    dico_nb_struct_0_aminor.update({elt[0] : 1})
                else :
                    dico_nb_struct_0_aminor[elt[0]] += 1
                liste_hors_aminor.append(elt)
            elif not elt[4] :
                if elt[0] not in dico_nb_struct_1_aminor :
                    dico_nb_struct_1_aminor.update({elt[0] : 1})
                else :
                    dico_nb_struct_1_aminor[elt[0]] += 1
                liste_hors_tot_aminor.append(elt)
    
    somme_1_aminor = 0            
    for elt in dico_nb_struct_1_aminor :
        somme_1_aminor += dico_nb_struct_1_aminor[elt]
    moy_1_aminor = somme_1_aminor/max(1, len(dico_nb_struct_1_aminor))
   
    somme_0_aminor = 0            
    for elt in dico_nb_struct_0_aminor :
        somme_0_aminor += dico_nb_struct_0_aminor[elt]
    moy_0_aminor = somme_0_aminor/max(1,len(dico_nb_struct_0_aminor))
    
    #csvwriter.writerow([num_cluster, len(cluster), len(dico_nb_struct), moy_1_aminor, moy_0_aminor])
    if num_cluster == 60 :
        print(dico_nb_struct)
        print(dico_nb_struct_0_aminor)
        print(dico_nb_struct_1_aminor)
        print(liste_hors_aminor[0:10])
        #exit()
#         print(liste_hors_aminor)
#         print(len(liste_hors_aminor))
#         
#         print(liste_hors_tot_aminor)
#         print(len(liste_hors_tot_aminor))
        
def sequence_apres_structure(cluster, dico_struct_cluster, graphes):
    
    liste_graphes  = []
    for elt in dico_struct_cluster :
        if (elt[0].split(",")[0][2:len(elt[0].split(",")[0])-1], elt[0].split(",")[1][2:len(elt[0].split(",")[1])-2]) not in liste_graphes :
            liste_graphes.append((elt[0].split(",")[0][2:len(elt[0].split(",")[0])-1], elt[0].split(",")[1][2:len(elt[0].split(",")[1])-2]))
    
    print(liste_graphes)
    #exit()
    dico_seq = donne_seq(graphes, liste_graphes)

    
    graphe_commun_global, dico_graphes, graphe_commun = recherche_graphe_commun([4], cluster, 0.6, "ramou", "Graphes_communs_mai_2020/extensions", 60)

    liste_seq_13, pos_aminor_seq_13, liste_seq_24, pos_aminor_seq_24 = etablissement_expr_reg(dico_graphes, graphe_commun)
    
    expr_reg_13 = []
    for cle in liste_seq_13 :
        for elt in liste_seq_13[cle] :
            if elt not in expr_reg_13 :
                expr_reg_13.append(elt)
    print(liste_seq_13)
    
    expr_reg_24 = []
    for cle in liste_seq_24 :
        for elt in liste_seq_24[cle] :
            if elt not in expr_reg_24 :
                expr_reg_24.append(elt)
    print(liste_seq_24)
    
    print(expr_reg_13)
    print(expr_reg_24)
    
    liste_trouve_1 = []
    compteur_1 = 0
    compteur_2 = 0
    compteur = 0
    dico_par_elt = {}
    for elt in dico_struct_cluster :
        print(elt)
        
        compteur += 1
        compteur_1 = False
        compteur_2 = False
        
        nom_graphe = (elt[0].split(",")[0][2:len(elt[0].split(",")[0])-1], elt[0].split(",")[1][2:len(elt[0].split(",")[1])-2])
        #exit()
        
        for expr in expr_reg_13 :
            match_it_13 = re.finditer(expr, dico_seq[nom_graphe][int(elt[1])-1:int(elt[1])+7])
            print(match_it_13)
            
            for m in match_it_13 :
                pos = m.span()
                print(pos)

                compteur_1 = True
                
        #exit()
        
        for expr in expr_reg_24 :
            match_it_24 = re.finditer(expr, dico_seq[nom_graphe][int(elt[2])-1:int(elt[2])+7])
            
            for m in match_it_24 :
                pos = m.span()
                print(pos)
                
                compteur_2 = True
                
        if compteur_1 and compteur_2 :
            if nom_graphe not in dico_par_elt :
                dico_par_elt.update({nom_graphe : 1})
            else :
                dico_par_elt[nom_graphe] += 1
        
    print(compteur_1)
    print(compteur_2)
    print(compteur)
    print(dico_par_elt)
        #exit()

 
if __name__ == '__main__' :
    
#     with open("Resultats_sequences/data_carnaval2_withnear_v3.137.pickle", 'rb') as fichier_pickle_1 :
#         mon_depickler = pickle.Unpickler(fichier_pickle_1)
#         graphes = mon_depickler.load()
#         print(len(graphes))
         
#     with open("Resultats_sequences/dico_organismes_2.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         dico_org =  mon_depickler.load()
#  
#         compteur = 0
#         for elt in dico_org :
#             if dico_org[elt] == None : 
#                 print(elt)
#                 compteur += 1
#         print(compteur)
#         exit()
#         
#         compteur = 0
#         for elt in dico_org :
#             if dico_org[elt] == None :
#                 dico_org[elt] = recup_org_gen(elt[0], elt[1])
#             print(compteur)
#             compteur += 1
#              
#             
#         with open("Resultats_sequences/dico_organismes_2.pickle", 'wb') as fichier_pickle :
#             mon_pickler = pickle.Pickler(fichier_pickle)
#             mon_pickler.dump(dico_org)
#            
#             exit()
               
#         exit()



            
    
#     liste_aminor = recup_liste_pos_aminor()
#     
#     dico_aminor = {}
#     for elt in liste_aminor :
#         if elt[0] not in dico_aminor :
#             dico_aminor.update({elt[0] : [elt[1]]})
#         else :
#             dico_aminor[elt[0]].append(elt[1])
#             
#     for elt in dico_aminor :
#         print(elt, dico_aminor[elt])
#         
#     #exit()
#             
#             
#     analyse_resultats_struct(dico_aminor, [])
#     exit()
    
#     with open("/media/coline/Maxtor/Resultats/clustering_perez_tot_new_data_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
#         mon_depickler = pickle.Unpickler(fichier_sortie)
#         clusters = mon_depickler.load()
#          
#     with open("Resultats_sequences/data_carnaval2_withnear_v3.137.pickle", 'rb') as fichier_pickle_1 :
#         mon_depickler = pickle.Unpickler(fichier_pickle_1)
#         graphes = mon_depickler.load()
#          
#         test(clusters, graphes)
#         exit()
    

    
#     with open("Resultats_sequences/data_carnaval2_withnear_v3.137.pickle", 'rb') as fichier_pickle_1 :
#         mon_depickler = pickle.Unpickler(fichier_pickle_1)
#         graphes = mon_depickler.load()
#     
#     with open("/media/coline/Maxtor/Resultats/clustering_perez_tot_new_data_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
#         mon_depickler = pickle.Unpickler(fichier_sortie)
#         clusters = mon_depickler.load()
#         
#     with open("Resultats_sequences/dico_cluster_60_que_1.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         dico_cluster = mon_depickler.load()
#     
#     sequence_apres_structure(clusters[59], dico_cluster, graphes)
#     
#     exit()
    
#     with open("Resultats_sequences/dico_recherche_struct_semi_aminor.pickle", 'rb') as fichier_pickle_3 : 
#         mon_pickler = pickle.Unpickler(fichier_pickle_3)
#         dico_recherche_struct_semi_aminor = mon_pickler.load()
#         print(list(dico_recherche_struct_semi_aminor.keys())[len(dico_recherche_struct_semi_aminor)-1])
#         print(dico_recherche_struct_semi_aminor[list(dico_recherche_struct_semi_aminor.keys())[len(dico_recherche_struct_semi_aminor)-1]])
#         print((('1VQP', '0'), 2357, 2293, 4, '[35, 29, 26, 3, 1, 21, 24, 25]', '[-1, 7, 6, 2, 4, 13, 15, 17]') in dico_recherche_struct_semi_aminor)
        #exit()
    liste_num_ARN = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]

    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/motif_A-minor/liste_representant_%s.pickle"%(elt), 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] != None and resolutions[element[0]] <= 3 and element[0] != '6hrm': # and element not in liste_pbs :
                            if (elt, element) not in liste_tout :
                                liste_tout.append((elt, element))
    print(liste_tout)
    liste_non_redondant = non_redondantes_seq_RNA3Dmotifatlas_new(liste_tout)
    #exit()
    liste_aminor = recup_liste_pos_aminor()
    
    with open("Resultats_sequences/data_carnaval2_withnear_v3.137.pickle", 'rb') as fichier_pickle_1 :
        mon_depickler = pickle.Unpickler(fichier_pickle_1)
        graphes = mon_depickler.load()
        
        print(('4V9Q', 'DV') in graphes)
        #exit()
    
    with open("/media/coline/Maxtor/Resultats/clustering_perez_tot_new_data_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        clusters = mon_depickler.load()
     
        recherche_expr_reg(clusters, liste_non_redondant, graphes, liste_aminor)
        exit()
#     exit()
#     liste_non_redondant = non_redondantes_seq_RNA3Dmotifatlas(liste_tout)
#     print(liste_non_redondant)
    
#     with open("Resultats_sequences/liste_non_redondant_struct.pickle", 'wb') as fichier_redondance :
#         mon_pickler = pickle.Pickler(fichier_redondance)
#         mon_pickler.dump(liste_non_redondant)
    #exit()
    
#     with open("Graphs/4w2g.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         graphe = mon_depickler.load()
#         print(graphe[('BA',2535)])
#         exit()
#         
#     with open("Nouvelles_donnees/fichier_5dm6_8_10.pickle", 'rb') as fichier_pickle :
#         mon_depickler = pickle.Unpickler(fichier_pickle)
#         extension = mon_depickler.load()
#                      
#         for noeud, data in extension.nodes(data=True) :
#             print(noeud, data)
#                        
#         for u,v,data in extension.edges(data=True) :
#             print(u,v,data)
#                          
#         exit()
            
    
#     
#     compteur = 0
#     for elt in liste_aminor :
#         if elt[0] == '4y4o' :
#             compteur += 1
#             if compteur == 21 :
#                 print(elt)
#              
#     exit()
    
#     with open("Resultats_sequences/liste_non_redondant.pickle", 'rb') as fichier_redondance :
#         mon_depickler = pickle.Unpickler(fichier_redondance)
#         liste_non_redondant = mon_depickler.load()
        
    with open("Resultats_sequences/data_carnaval2_withnear_v3.137.pickle", 'rb') as fichier_pickle_1 :
        mon_depickler = pickle.Unpickler(fichier_pickle_1)
        graphes = mon_depickler.load()
            
        
        liste_nom_graphes_aminor = []
        
        
        ## recherche des structures PDB qui font partie de notre jeu de données d'A-minor
        for cle in liste_tout :
            with open("Nouvelles_donnees/fichier_%s_%s_10.pickle"%(cle[1][0], cle[1][1]), "rb") as fichier_graphe :
                mon_depickler = pickle.Unpickler(fichier_graphe)
                extension = mon_depickler.load()
                
                if cle[1][0] == '4v67' :
                    print(extension.nodes[1]["num_ch"])
                    print(extension.nodes[3]["num_ch"])
                    print(extension.nodes[5]["num_ch"])
                
                if (cle[1][0].upper(), extension.nodes[1]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[1]["num_ch"]))
                
                if (cle[1][0].upper(), extension.nodes[3]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[3]["num_ch"]))
                
                if (cle[1][0].upper(), extension.nodes[5]["num_ch"]) not in liste_nom_graphes_aminor :
                    liste_nom_graphes_aminor.append((cle[1][0].upper(), extension.nodes[5]["num_ch"]))
        
#         print(graphes[('1VQ8', '0')][1106])
#         exit()
    print(liste_nom_graphes_aminor)
    print(len(liste_nom_graphes_aminor))
    
    compteur = 0
    for elt in liste_nom_graphes_aminor :
        compteur += 1
        if elt == ('4V51', 'BA') :
            print(compteur)
    print(compteur)
    #exit()
    
    compteur = 0
    ## on en enlève les doublons
    dico_seq = donne_seq(graphes, liste_nom_graphes_aminor)
    liste_en_double = [] 
    for i in range(len(liste_nom_graphes_aminor)) :
        #print(elt)
        elt = liste_nom_graphes_aminor[i]
        
        for j in range(i+1, len(liste_nom_graphes_aminor)) :
            elt2 = liste_nom_graphes_aminor[j]
            
            if dico_seq[elt] == dico_seq[elt2] :
                print(elt, elt2)
                if elt2 not in liste_en_double :
                    liste_en_double.append(elt2)
                compteur += 1
    print(compteur)
    print(liste_en_double)
    print(len(liste_en_double))
    
    for elt in liste_en_double :
        liste_nom_graphes_aminor.remove(elt)
    print(len(liste_nom_graphes_aminor))
    #exit()
    
    
    liste_aminor = recup_liste_pos_aminor()

        #exit()
    #exit()  
    with open("/media/coline/Maxtor/Resultats/clustering_perez_tot_new_data_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        clusters = mon_depickler.load()
        recherche_struct(clusters, liste_non_redondant, graphes, liste_aminor)
        
        exit()
#     liste_aminor = recup_liste_pos_aminor()
#     with open("Resultats_sequences/graphs_all_2020_07_with_SSEs.pickle", 'rb') as fichier_pickle_1 :
#         mon_depickler = pickle.Unpickler(fichier_pickle_1)
#         graphes = mon_depickler.load()
#             
#         print(graphes[('5DM6', 'X')].nodes[1368])
#         print(graphes[('5DM6', 'X')].nodes[1369])
#         print(graphes[('5DM6', 'X')].nodes[1370])
#         print(graphes[('5DM6', 'X')].nodes[1371])
#         print(graphes[('5DM6', 'X')].nodes[1372])
#         print(graphes[('5DM6', 'X')].nodes[1373])
#         print(graphes[('5DM6', 'X')].nodes[1374])
#         print(graphes[('5DM6', 'X')].nodes[1375])
#         exit()
# #         print(graphes[('1C2W', 'B')].nodes[794])
#         print(graphes[('4IOA', 'X')][475])
#         print(graphes[('4IOA', 'X')][474])
# #         print(graphes[('4IOA', 'X')][794])
#         #exit()
# # #  
#         for elt in liste_aminor :
#             if elt[0] == '4ioa' :
#                 print(elt)
# # #       
# # #           
#         exit()
# # #                 
# # # #     
    with open("Resultats_sequences/communs_aminor.pickle", 'rb') as fichier_aminor :
            mon_depickler = pickle.Unpickler(fichier_aminor)
            communs_aminor = mon_depickler.load()
    with open("Resultats_sequences/communs_non_aminor.pickle", 'rb') as fichier_non_aminor :
            mon_depickler = pickle.Unpickler(fichier_non_aminor)
            communs_non_aminor = mon_depickler.load() 
#             
#     compte = 0        
#     compteur = 0
#     for elt in communs_non_aminor :
#         print(elt)
#         for elt_2 in elt :
#             if len(elt_2) != 4 :
#                 print(elt_2)
#                 print(compteur)
#                 compte += 1
#         compteur += 1
#     print(compte)
#     #exit()
#     
    with open("Resultats_sequences/par_sequence_non_aminor_struct_second.csv", 'w', newline="") as fichier_struct :
        csvwriter = csv.writer(fichier_struct)
        csvwriter.writerow(["Nb struct second", "Elements"])
        for elt in communs_non_aminor :
            compter = 0
            for elt2 in elt :
                if len(elt2) == 4 :
                    compter += 1
            csvwriter.writerow([compter, elt])
    #exit()
#     compte = 0        
#     for elt in communs_aminor :
#         print(elt)
#         compte += len(elt)
#     print(compte)
#     exit()
        
#     liste_aminor = recup_liste_pos_aminor()
#       
#     for elt in liste_aminor :
#         if elt[0] == '4v4b' :
#             print(elt)
#       
#     exit()
                
    liste_num_ARN = ["23S", "18S", "16S","Ribozyme", "Riboswitch", "SRP", "28S", "25S", "Intron", "arnt_16S_arnm", "arnt_16S"]

    liste_tout = []
    for elt in liste_num_ARN :
        with open("resolutions.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            resolutions = mon_depickler.load()
             
            with open("Nouvelles_donnees/motif_A-minor/liste_representant_%s.pickle"%(elt), 'rb') as fichier_sortie :
                    mon_depickler = pickle.Unpickler(fichier_sortie)
                    liste_a_garder_noms = mon_depickler.load()    
                       
                    for element in liste_a_garder_noms :
                        if resolutions[element[0]] != None and resolutions[element[0]] <= 3 and element[0] != '6hrm': # and element not in liste_pbs :
                            if (elt, element) not in liste_tout :
                                liste_tout.append((elt, element))
    print(liste_tout)


    
    
    
#     liste_non_redondant = non_redondantes_seq_RNA3Dmotifatlas(liste_tout)
#     if ('5MDV', '1') in liste_non_redondant :
#         print("gros tas")
#     exit()
#     
#     liste_aminor = recup_liste_pos_aminor()
#     for cle in liste_aminor :
#         if cle[0] == '5j7l' and cle[1][0][0] == 'CA' :
#             print(cle)
#     exit()        
    #diff_occ_aminor("Resultats_sequences/fichier_csv_nb_seq_aminor_par_sequence_4.csv", "Resultats_sequences/fichier_csv_nb_seq_aminor_par_sequence_2.csv")
    #exit()
    
    
                     
    
    
#     with open("Graphs/4ybb.pickle", 'rb') as fichier_graphe :
#         mon_depickler = pickle.Unpickler(fichier_graphe)
#         graphe = mon_depickler.load()
#         
#         for noeud, data in graphe.nodes(data=True) :
#             if noeud[0] == 'CA' :
#                 print(noeud, data)
#                 
#     with open("Nouvelles_donnees/fichier_4ybb_14_10.pickle", 'rb') as fichier_graphe1 :
#             mon_depickler = pickle.Unpickler(fichier_graphe1)
#             extension = mon_depickler.load()
# 
#             seqs = recup_sequences_autour_motif("4ybb", extension, 30, [[3,1], [2,4], [5,5]])
#             print(seqs)
#     exit()

    #stocker_type_arn()
    #exit()
#     with open("Resultats/clustering_perez_tot_new_data_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
#         mon_depickler = pickle.Unpickler(fichier_sortie)
#         clusters = mon_depickler.load()
#         
#         compteur = 1
#         for cluster in clusters :
#     
#             print(compteur)
#             print(cluster)
#             
#             compteur += 1
#     exit()   
#     for i in range(1, 20) :
#         with open("Nouvelles_donnees/fichier_5j7l_%s_10.pickle"%i, 'rb') as fichier_ext :
#             mon_depickler = pickle.Unpickler(fichier_ext)
#             extension = mon_depickler.load()
#             print(extension.nodes[1])
#             print(recup_sequences_autour_motif('5j7l', extension, 3, [[3,1], [2,4], [5, 5]]))
#         
#     exit()
    
    ### cherche les occurrences d'une séquence aléatoire de 8 nts
    with open("Resultats_sequences/graphs_all_2020_07_with_SSEs.pickle", 'rb') as fichier_pickle_1 :
        mon_depickler = pickle.Unpickler(fichier_pickle_1)
        graphes = mon_depickler.load()
         
    with open("type_arn_tot.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        type_arn = mon_depickler.load()
#         
#         ### cherche les structures non redondantes de la PDB selon la méthode et les résultats de RNA 3D motif Atlas
#         liste_non_redondant = non_redondantes_seq_RNA3Dmotifatlas()
#         
#         
#         with open("Resultats_sequences/sequence_aleatoire.csv", 'w', newline="") as fichier_csv_alea :
#             csvwriter = csv.writer(fichier_csv_alea)
#             csvwriter.writerow(["Seq alea", "Arn seq alea", "Nombre d'occurrences"])
#             
#             for _ in range(100) :
#                 
#                 seq_alea, cle_alea, nb_al = sequence_aleatoire(graphes, liste_non_redondant)
#                 nb_occ = recherche_seq_alea(graphes, liste_non_redondant,seq_alea, cle_alea, nb_al)
#                 
#                 csvwriter.writerow([seq_alea, cle_alea, type_arn[cle_alea] ,nb_occ ])
#         exit()
        
    #liste_non_redondant = non_redondantes_seq_RNA3Dmotifatlas(liste_tout)
    with open("Resultats_sequences/liste_non_redondant.pickle", 'rb') as fichier_redondance :
        mon_depickler = pickle.Unpickler(fichier_redondance)
        liste_non_redondant = mon_depickler.load()
    
    #exit()
    
    ### cherche les positions de tous les motifs A-minor
    liste_aminor = recup_liste_pos_aminor()
    
    with open("Resultats/clustering_perez_tot_new_data_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
        mon_depickler = pickle.Unpickler(fichier_sortie)
        clusters = mon_depickler.load()
    
    
    recherche_expr_reg(clusters, liste_non_redondant, graphes, liste_aminor)
    #recherche_tot(graphes, liste_non_redondant, clusters, liste_aminor)
    exit()
    
    with open("Resultats_sequences/nrlist_3.136_all.csv", 'r') as fichier_csv_motifs_atlas:
        csvreader = csv.reader(fichier_csv_motifs_atlas)
        dico_redondance = {}
        for row in csvreader :
            if '+' not in row[1] : ## cas de base
                
                repr = (row[1].split("|")[0], row[1].split("|")[2])
                
                representes = []
                autres = row[2].split(",")
                for elt in autres :
                    representes.append((elt.split("|")[0], elt.split("|")[2]))
                dico_redondance.update({repr:representes})
            else : ## cas où le représentant de la classe est composé de plusieurs structures (je sais pas pourquoi ?)
                compteur = 0
                repr = tuple(row[1])
                bizarre = False
                representes = []
                autres = row[2].split(",")
                for elt in autres :
                    for elt_2 in elt.split('+') :
                        representes.append((elt_2.split("|")[0], elt_2.split("|")[2]))

                dico_redondance.update({repr:representes})
    
    with open("Resultats_sequences/resultats.pickle", 'rb') as fichier_pickle :
        mon_depickler = pickle.Unpickler(fichier_pickle)
        res = mon_depickler.load()
        
        print(list(res.keys())[0])
        print(res[list(res.keys())[0]])
        
        dico_nb_occ = {}
        dico_struct_occ = {}
        aminors = {}
        for elt in res.keys() :
                if elt[0] == "cluster_55_nb_elts_3" and elt[1][0] == '4Y4O' :
                    print(elt)
                    print(res[elt])
                    #exit()
                #if elt[0] == "cluster_56_nb_elts_6" :
                #    exit()
                ## elt = une clé (cluster, (pdb_id, chain))
                
                liste_a_chercher = []
                trouve = False
                for cle in dico_redondance :
                    if elt[1] == cle or elt[1] in cle :
                        liste_a_chercher.append(dico_redondance[cle])
                        trouve = True
                    elif elt[1] in dico_redondance[cle] :
                        liste_a_chercher.append(dico_redondance[cle])
                        trouve = True
                if not trouve :
                    liste_a_chercher.append([elt[1]])
                
                compteur = 0
                for occ_seq in res[elt] :
                    ## occ_seq = un doublet ((a_minor_seq1, a_minor_seq2, a_minor_seq3), [[start1, .....], [start2, ...], [start3, ...]])
                        nb_trouve = [0,0,0]
                        
                        

                        
                        if (elt[0], occ_seq[0]) not in dico_nb_occ.keys() :
                            dico_nb_occ.update({(elt[0], occ_seq[0]) : [0,0,0]})
                            #dico_nb_occ.update({(elt[0], compteur) : [[],[],[]]})
                            dico_struct_occ.update({(elt[0], occ_seq[0]) : [[], [], []]})
                            aminors.update({(elt[0], occ_seq[0]) : [[], [], []]})
                        print((elt[0], compteur))
                        
                        
                        if elt[1] in liste_non_redondant :
#                                 print("entree")
#                                 print(occ_seq)
    #                     print(occ_seq)
                            ### recherche s'il y a un motif A-minor à l'endroit de occ_seq
    #                             if i == 0 :
    #                                 print(compte_occ)
    #                                 for cle in dico_seqs :
    #                                     if len(dico_seqs[cle]) > 0 :
    #                                         print(cle)
    #                                         print(dico_seqs[cle])
                                    #print(dico_seqs)
                                    
                                
                                for i in range(3) :
                                    for liste in liste_a_chercher : 
                                        nb_trouve_temp = 0
                                        aminors_temp = []
                                        existe = False
                                        for l in liste :
                                            if elt[1] != l : 
                                                if l in graphes.keys() :
                                                    seq = ""
                                                    for j in range(1, graphes[l].number_of_nodes()+1) :
                                                        seq += graphes[l].nodes[j]["nt"]
                                                    
                                                    #print(seq)
                                                    
                                                    occ_temp = []
                                                    start = -1
                                                    while start != 0 :
                                                        start = seq.find(occ_seq[0][i], max(start,0))
                                                        if start != -1 :
                                                            occ_temp.append(start)
                                                        start += 1
                                                    if len(occ_temp) > 0 :
                                                    
                                                
                                                        for aminor in liste_aminor :
                            #                                 print(aminor)
                            #                                 print(occ)
                                                            if l[0] == aminor[0].upper() :
                                                                     
                                                                    if l[1] == aminor[1][i][0] :
              
                                                                        for occ2 in occ_temp :
                                                                            if occ2 >= aminor[1][i][1] - 5 and occ2 <= aminor[1][i][1] :
                                   
                                                                                nb_trouve_temp += 1
                                                                                aminors_temp.append((l, occ2, occ_seq[0][i]))
                                                                                existe = True
                                            else :
                                                    for aminor in liste_aminor :
                        #                                 print(aminor)
                        #                                 print(occ)
                                                        if elt[1][0] == aminor[0].upper() :
                                                                 
                                                                if elt[1][1] == aminor[1][i][0] :
                        
                                                                    for occ2 in occ_seq[1][i] :
#                                                                         print(occ2)
                                                                        if int(occ2) >= aminor[1][i][1] - 4 and int(occ2) <= aminor[1][i][1] :
                                                                           
                                                                            nb_trouve_temp += 1
                                                                            aminors_temp.append((elt[1], occ2, occ_seq[0][i]))
                                                                            existe = True
                                            if existe :
                                                nb_trouve[i] += nb_trouve_temp
                                                aminors[(elt[0], occ_seq[0])][i].extend(aminors_temp)
                                                break
                                    
                               
                        
#                             dico_nb_occ[(elt[0], occ_seq[0])][0] += len(occ_seq[1][0]) - nb_trouve_1
#                             dico_nb_occ[(elt[0], occ_seq[0])][1] += len(occ_seq[1][1]) - nb_trouve_2
#                             dico_nb_occ[(elt[0], occ_seq[0])][2] += len(occ_seq[1][2]) - nb_trouve_3
                            
                                dico_nb_occ[(elt[0], occ_seq[0])][0] += nb_trouve[0]
                                dico_nb_occ[(elt[0], occ_seq[0])][1] += nb_trouve[1]
                                dico_nb_occ[(elt[0], occ_seq[0])][2] += nb_trouve[2]
                                
                                if len(occ_seq[1][0]) - nb_trouve[0] > 0 :
                                    dico_struct_occ[(elt[0], occ_seq[0])][0].append((elt[1], occ_seq[1][0], type_arn[elt[1]]))
                                    
                                if len(occ_seq[1][1]) - nb_trouve[1] > 0 :
                                    dico_struct_occ[(elt[0], occ_seq[0])][1].append((elt[1], occ_seq[1][1], type_arn[elt[1]]))
                                    
                                if len(occ_seq[1][2]) - nb_trouve[2] > 0 :
                                    dico_struct_occ[(elt[0], occ_seq[0])][2].append((elt[1], occ_seq[1][2], type_arn[elt[1]] ))
                            
                            
                        
                            #dico_nb_occ[(elt[0], compteur)][0].append(len(occ_seq[1][0]) - nb_trouve_1)
                            #dico_nb_occ[(elt[0], compteur)][1].append(len(occ_seq[1][1]) - nb_trouve_2)
                            #dico_nb_occ[(elt[0], compteur)][2].append(len(occ_seq[1][2]) - nb_trouve_3)
                        compteur += 1
                
                
                #if elt[1] == ('5DM6', 'X') :
                #    exit()
        
        #for elt in dico_nb_occ :
        #    print(elt)
        #    print(dico_nb_occ[elt])
        #exit()    
        
        
        
        with open("Resultats_sequences/graphs_all_2020_07_with_SSEs.pickle", 'rb') as fichier_graphe2 :
            mon_depickler = pickle.Unpickler(fichier_graphe2)
            graphes = mon_depickler.load()
        
        
        with open("Resultats/clustering_perez_tot_new_data_sim_par_branche_0.65.pickle", 'rb') as fichier_sortie :
            mon_depickler = pickle.Unpickler(fichier_sortie)
            clusters = mon_depickler.load()
        
        ### stocke le nombre d'occurrences par séquence (associée à un cluster indiqué)
        with open("Resultats_sequences/fichier_csv_nb_seq_aminor_par_sequence_2.csv", 'w', newline="") as fichier_csv :
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["Numéro cluster", "Sequences", "Elts cluster et position aminor", "Nb seq 1", "Nb seq 2", "Nb seq 3", "A-minor"])
            
            for elt in dico_nb_occ :
                num_cluster = int(elt[0].split("_")[1])
                
                liste_cluster = []
                for aminor in clusters[num_cluster-1] :
                    with open("Nouvelles_donnees/fichier_%s_%s_10.pickle"%(aminor[0], aminor[1]), "rb") as fichier_ext :
                        mon_depickler = pickle.Unpickler(fichier_ext)
                        extension = mon_depickler.load()
                        
                        liste_cluster.append((aminor, (extension.nodes[1]["num_ch"], extension.nodes[1]["position"][0]), (extension.nodes[2]["num_ch"], extension.nodes[2]["position"][0]), (extension.nodes[5]["num_ch"], extension.nodes[5]["position"][0])))
                    
                
                csvwriter.writerow([elt[0], elt[1], liste_cluster, dico_nb_occ[elt][0], dico_nb_occ[elt][1], dico_nb_occ[elt][2], aminors[elt]])
        exit()      
        ### stocke le nombre d'occurrences par séquence dans des structures non homologues 
        ### et donne l'identifiant des structures dans lesquelles sont retrouvées les occurrences
        with open("Resultats_sequences/fichier_csv_nb_seq_sans_aminor_par_sequences_par_nom_structure.csv", 'w', newline="") as fichier_csv :
            csvwriter = csv.writer(fichier_csv)
            csvwriter.writerow(["Numéro cluster", "Sequences", "Elts cluster et position aminor","NB seq 1 non homol", "NB seq 2 non homol", "NB seq 3 non homol","Seq 1", "Seq 2", "Seq 3"])
            
            
            
            
            for elt in dico_struct_occ :
                num_cluster = int(elt[0].split("_")[1])
                
                liste_cluster = []
                for aminor in clusters[num_cluster-1] :
                    with open("Nouvelles_donnees/fichier_%s_%s_10.pickle"%(aminor[0], aminor[1]), "rb") as fichier_ext :
                            mon_depickler = pickle.Unpickler(fichier_ext)
                            extension = mon_depickler.load()
                            
                            liste_cluster.append((aminor, (extension.nodes[1]["num_ch"], extension.nodes[1]["position"][0]), (extension.nodes[2]["num_ch"], extension.nodes[2]["position"][0]), (extension.nodes[5]["num_ch"], extension.nodes[5]["position"][0])))
                
                
                nb_non_homol_1 = alignement_global(clusters[num_cluster-1], dico_struct_occ[elt][0], 1, graphes)
                nb_non_homol_2 = alignement_global(clusters[num_cluster-1], dico_struct_occ[elt][1], 2, graphes)
                nb_non_homol_3 = alignement_global(clusters[num_cluster-1], dico_struct_occ[elt][2], 3, graphes)
                

           
                
                csvwriter.writerow([elt[0], elt[1], liste_cluster, nb_non_homol_1, nb_non_homol_2, nb_non_homol_3, dico_struct_occ[elt][0], dico_struct_occ[elt][1], dico_struct_occ[elt][2]])
                #exit()
                