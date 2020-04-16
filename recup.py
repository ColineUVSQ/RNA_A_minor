'''
Created on 5 déc. 2018

@author: coline
'''
from audioop import max
with open("fichier_similarite_test_3.txt", 'r') as fichier :
    ligne = " "
    dict = {}
    while len(ligne) != 0 :
        ligne1 = fichier.readline()
        ligne2 = fichier.readline()
        element1 = ligne1[:len(ligne1)-1]
        element2 = ligne2[:len(ligne2)-1]
        print(element1)
        print(element2)
        
        for i in range(2) :
            fichier.readline()
            
        ligne = fichier.readline()
        sim = ligne[:len(ligne)-1].split(":")[1]
        dict.update({(element1, element2) : sim})
        ligne = fichier.readline()
    
    print(dict)
    
    compteur = 0
    max = 0 
    pos_max = 0
    for elt in dict.keys() :
        if elt[0] == 'fichier_1FJG_A_48_8' or elt[1] == 'fichier_1FJG_A_48_8' :
            print(elt)
            print(dict[elt])
            compteur += 1
            if float(dict[elt]) > max and float(dict[elt]) < 1.0:
                max = float(dict[elt])
                pos_max = elt
            
    print(compteur)
    print(pos_max)
    print(max)
        
    
    