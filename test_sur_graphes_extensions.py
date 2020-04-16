'''
Created on 5 nov. 2018

@author: Coline Gi
'''
import pickle
import os

def nb_nts_chaines(fichier_compter, G):
    fichier_compter.write(str(element)+'\n')
    
    liaison_B53 = True
    compteur = 1
    compteur_nts  = 1
    nb_liaison_can = 1  
    chaine = [1]
                          
    while liaison_B53 :
        liaison_B53 = False
        nb_liaison_can = 0
        temp = compteur
        for voisin in G[compteur] :
                #print(compteur)
                #print(voisin)
                #print(G.edges.data())
            if voisin not in chaine and G[compteur][voisin]["label"] == 'B53' :
                liaison_B53 = True
                compteur_nts += G.nodes[voisin]["poids"]
                temp = voisin
                nb_liaison_can += 1
                chaine.append(voisin)
        #fichier_compter.write(str(nb_liaison_can)+'\n')
        if nb_liaison_can > 1 and compteur != 1:
            print("probleme")
            fichier_compter.write(str(compteur) + " " + str(nb_liaison_can) +" probleme"+'\n')
        compteur = temp
        
    fichier_compter.write("total 1 : "+str(compteur_nts-1) + "\n")
    print(compteur_nts)
    
    liaison_B53 = True
    compteur = 2
    compteur_nts  = 1
    nb_liaison_can = 1  
    chaine = [2]
                          
    while liaison_B53 :
        liaison_B53 = False
        nb_liaison_can = 0
        temp = compteur
        for voisin in G[compteur] :
                #print(compteur)
                #print(voisin)
                #print(G.edges.data())
            if voisin not in chaine and G[compteur][voisin]["label"] == 'B53' :
                liaison_B53 = True
                compteur_nts += G.nodes[voisin]["poids"]
                temp = voisin
                nb_liaison_can += 1
                chaine.append(voisin)
        #fichier_compter.write(str(nb_liaison_can)+'\n')
        if nb_liaison_can > 1 and compteur != 2:
            print("probleme")
            fichier_compter.write(str(compteur) + " " + str(nb_liaison_can) +" probleme"+'\n')
        compteur = temp
        
    fichier_compter.write("total 2 : "+str(compteur_nts-1) + "\n")
    print(compteur_nts)
    
    liaison_B53 = True
    compteur = 3
    compteur_nts  = 1
    nb_liaison_can = 1  
    chaine = [3]
                          
    while liaison_B53 :
        liaison_B53 = False
        nb_liaison_can = 0
        temp = compteur
        for voisin in G[compteur] :
                #print(compteur)
                #print(voisin)
                #print(G.edges.data())
            if voisin not in chaine and G[compteur][voisin]["label"] == 'B53' :
                liaison_B53 = True
                compteur_nts += G.nodes[voisin]["poids"]
                temp = voisin
                nb_liaison_can += 1
                chaine.append(voisin)
        #fichier_compter.write(str(nb_liaison_can)+'\n')
        if nb_liaison_can > 1 and compteur != 3:
            print("probleme")
            fichier_compter.write(str(compteur) + " " + str(nb_liaison_can) +" probleme"+'\n')
        compteur = temp
        
    fichier_compter.write("total 3 : "+str(compteur_nts-1) + "\n")
    print(compteur_nts)
    
    liaison_B53 = True
    compteur = 4
    compteur_nts  = 1
    nb_liaison_can = 1  
    chaine = [4]
                          
    while liaison_B53 :
        liaison_B53 = False
        nb_liaison_can = 0
        temp = compteur
        for voisin in G[compteur] :
                #print(compteur)
                #print(voisin)
                #print(G.edges.data())
            if voisin not in chaine and G[compteur][voisin]["label"] == 'B53' :
                liaison_B53 = True
                compteur_nts += G.nodes[voisin]["poids"]
                temp = voisin
                nb_liaison_can += 1
                chaine.append(voisin)
        #fichier_compter.write(str(nb_liaison_can)+'\n')
        if nb_liaison_can > 1 and compteur != 4:
            print("probleme")
            fichier_compter.write(str(compteur) + " " + str(nb_liaison_can) +" probleme"+'\n')
        compteur = temp
        
    fichier_compter.write("total 4 : "+str(compteur_nts-1) + "\n")
    print(compteur_nts)

def recup_chaines(G):
    chaines = [[1]]
    for i in range(1,5) :
        compteur = i
        if i != 1 : chaines.append([i])
        liaison_B53 = True
        while liaison_B53 :
            liaison_B53 = False
            temp = compteur
            for voisin in G[compteur] :
                    #print(compteur)
                    #print(voisin)
                    #print(G.edges.data())
                if voisin not in [1,2,3,4] and voisin not in chaines[len(chaines)-1] and G[compteur][voisin]["label"] == 'B53' :
                    liaison_B53 = True
                    temp = voisin
                    chaines[len(chaines)-1].append(voisin)
            #fichier_compter.write(str(nb_liaison_can)+'\n')
            compteur = temp
    print(chaines)
    
    
with open("fichier_nb_sommets.txt", 'w') as fichier_compter :
        for element in os.listdir('graphes_extension/'):
            if "pickle" in element :
                with open("graphes_extension/"+element, 'rb') as fichier_entree :
                    print(element)
                    
                    mon_depickler = pickle.Unpickler(fichier_entree)
                    G = mon_depickler.load()
                    #print(G.nodes.data())

                    #nb_nts_chaines(fichier_compter, G)
                    
                    recup_chaines(G)
                    
                        
                    

                    