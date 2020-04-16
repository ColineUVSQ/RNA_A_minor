'''
Created on 18 oct. 2018

@author: Coline Gi

Stockage dans un fichier des nucleotides en position 2 et 4 du motif
pour chaque occurrence

'''


import os
import pickle

dico_classe = {}

for element in os.listdir('graphes_extension/'):
        if "pickle" in element :
            with open("graphes_extension/"+element, 'rb') as fichier_entree :
                print(element)
                
                mon_depickler = pickle.Unpickler(fichier_entree)
                G = mon_depickler.load()
                
                if (G.nodes[2]["nt"], G.nodes[4]["nt"]) not in dico_classe.keys() :
                    dico_classe.update({(G.nodes[2]["nt"], G.nodes[4]["nt"]) : [element]})
                else :
                    dico_classe[(G.nodes[2]["nt"], G.nodes[4]["nt"])].append(element)

with open("fichier_type_nts.txt", 'w') as fichier :
    fichier.write(str(dico_classe))