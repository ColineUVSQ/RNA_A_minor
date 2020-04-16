'''
Created on 12 déc. 2018

@author: coline
'''

import pickle
import os
from recup_data.constantes import EXTENSION_PATH_TAILLE, GROUPE_ARICH,\
    GROUPE_GNRA, NEW_EXTENSION_PATH_TAILLE


def lien_chaines_1_3(graphe):
    for noeud, data in graphe.nodes(data=True) :
        if noeud not in [1, 2, 3, 4, 5] :
            if 1 in data["chaine"] :
                for voisin in graphe[noeud] :
                    if 3 in graphe.nodes[voisin]["chaine"] :
                        return True
            if 3 in data["chaine"] :
                for voisin in graphe[noeud] :
                    if 1 in graphe.nodes[voisin]["chaine"] :
                        return True
    return False

def test():
    compteur = 0
    compteur_tot = 0
    for fic in os.listdir(EXTENSION_PATH_TAILLE%10) :
        if "pickle" in fic and len(fic.split("_")) == 5 and "graphe_comp" not in fic :
            
            with open(EXTENSION_PATH_TAILLE%10+fic, 'rb') as fichier_graphe :
                mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
                graphe = mon_depickler_graphe.load()
                print(fic)
                print(lien_chaines_1_3(graphe))
                if lien_chaines_1_3(graphe) :
                    compteur += 1
                compteur_tot += 1
    print(compteur)
    print(compteur_tot)

''' 12/09/19 '''
def version_new_data():
    
#     with open("graphs_2.92.pickle", 'rb') as fichier_tot :
#     mon_depickler = pickle.Unpickler(fichier_tot)
#     graphe_tot = mon_depickler.load()
    
    
    gnra = []
    arich = []
    gbulge = []
    compteur = 0
    for fic in os.listdir(NEW_EXTENSION_PATH_TAILLE) :
        if "pickle" in fic and "groupe" not in fic :
            
            with open("Graphs/"+fic.split("_")[1]+".pickle", 'rb') as fichier_tot :
                mon_depickler = pickle.Unpickler(fichier_tot)
                graphe_tot = mon_depickler.load()
            
            with open(NEW_EXTENSION_PATH_TAILLE+fic, 'rb') as fichier_graphe :
                mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
                graphe = mon_depickler_graphe.load()
                occ = fic.split("_")[1]
        
                
                #if int(num_motif) in [44,58,62,197,251,271,137,50,328,118,134,74] :
                compteur += 1
                print(compteur)
                print(fic) 
           
            
                pos_1 = graphe.nodes[1]["position"][0]
                num_ch_1 = graphe.nodes[1]["num_ch"]
                pos_3 = graphe.nodes[3]["position"][0]
                num_ch_3 = graphe.nodes[3]["num_ch"]
                
                if graphe_tot.nodes[(num_ch_1, pos_1)]["nt"] == 'A' and (graphe_tot.nodes[(num_ch_3,pos_3)]["nt"] == 'A' or graphe_tot.nodes[(num_ch_3,pos_3)]["nt"] == 'G') and pos_3-2 > 0 and graphe_tot.nodes[(num_ch_3, pos_3-2)]["nt"] == 'G'  : #and graphe_tot[num].nodes[pos_3-4]["part"] == "Stem" and graphe_tot[num].nodes[pos_1+2]["part"] == "Stem" :  
                    
                    boucle_de_4 = False
                    if pos_3-3 > 0 and ((num_ch_1, pos_1+1), (num_ch_3, pos_3-3)) in graphe_tot.edges() :
                        if graphe_tot[(num_ch_1, pos_1+1)][(num_ch_3, pos_3-3)]["label"] == 'CWW' :
                            boucle_de_4 = True
                    
                    pos_4_ok = False
                    pos_4 = graphe.nodes[4]["position"][0] 
                    num_ch_4 = graphe.nodes[4]["num_ch"]
                    #if (graphe_tot.nodes[(num_ch_4, pos_4)]["nt"] == 'C' or graphe_tot.nodes[(num_ch_4, pos_4)]["nt"] == 'U') :
                    for voisin in graphe_tot[(num_ch_4, pos_4)] :
                            if graphe_tot[(num_ch_4, pos_4)][voisin]["label"] == 'CWW' and (graphe_tot.nodes[voisin]["nt"] == 'G' or graphe_tot.nodes[voisin]["nt"] == 'A') :
                                pos_4_ok = True
                                
#                     liaison_ok = False
#                     for voisin in graphe_tot[(num_ch_1, pos_1)] :
#                         if pos_3-2 > 0 and voisin == pos_3-2 and graphe_tot[(num_ch_1, pos_1)][voisin]["label"] == 'THS' :
#                             liaison_ok = True
                    if boucle_de_4 and pos_4_ok :
                        gnra.append(fic)
                    #gnra.append(fic[8:len(fic)-7])
                        
#                 if graphe_tot.nodes[(num_ch_1,pos_1)]["nt"] == 'A' and graphe_tot.nodes[(num_ch_3, pos_3)]["nt"] == 'A' and not lien_chaines_1_3(graphe) :
#                     liaison_1_ok = False
#                     a_voisin_plus_1 = False
#                     for voisin in graphe_tot[(num_ch_1, pos_1)] :
#                         if graphe_tot[(num_ch_1, pos_1)][voisin]["label"] == 'THS' :
#                             liaison_1_ok = True
#                             print(graphe_tot.nodes[(voisin[0],voisin[1]+1)]["nt"])
#                             if graphe_tot.nodes[(voisin[0],voisin[1]+1)]["nt"] == 'A' :
#                                 a_voisin_plus_1 = True
# 
# 
#                             
#                     liaison_3_1_ok = False
#                     if pos_3-1 > 0 :
#                         for voisin in graphe_tot[(num_ch_3,pos_3-1)] :
#                             if graphe_tot[(num_ch_3, pos_3-1)][voisin]["label"] == 'TWH' :
#                                 liaison_3_1_ok = True
#                    # liaison_3_1_ok = True
#                             
#                     liaison_stem = False
#                     g_a_voisin_plus_1 = False
#                     for voisin in graphe_tot[(num_ch_1, pos_1+1)] :
#                         if graphe_tot[(num_ch_1, pos_1+1)][voisin]["label"] == 'CWW' :
#                             liaison_stem = True
#                             if graphe_tot.nodes[(voisin[0],voisin[1]+1)]["nt"] == 'G' and graphe_tot.nodes[(voisin[0],voisin[1]+2)]["nt"] == 'A' :
#                                 g_a_voisin_plus_1 = True
#                             
#                     
#                     if liaison_3_1_ok and liaison_stem and g_a_voisin_plus_1 and a_voisin_plus_1 and liaison_1_ok:
#                         arich.append(fic)
#                         
#                         
#                 if graphe_tot.nodes[(num_ch_1,pos_1)]["nt"] == 'A' and graphe_tot.nodes[(num_ch_3,pos_3)]["nt"] == 'A' :
#                     liaison_1_1_ok = False
#                     pos = -1
#                     for voisin in graphe_tot[(num_ch_1, pos_1)] :
#                         if graphe_tot[(num_ch_1, pos_1)][voisin]["label"] == 'THH' :
#                             liaison_1_1_ok = True
#                             pos = voisin
# 
#                     liaison_1_2_ok = False
#                     if pos != -1 :
#                         for voisin in graphe_tot[(pos[0], pos[1]+1)] :
#                             if graphe_tot[(pos[0], pos[1]+1)][voisin]["label"] == 'CSH' :
#                                 liaison_1_2_ok = True
# 
#                     liaison_3_ok = False
#                     for voisin in graphe_tot[(num_ch_3, pos_3)] :
#                         if graphe_tot[(num_ch_3, pos_3)][voisin]["label"] == 'THW' :
#                             liaison_3_ok = True
#                     
#                     liaison_1_plus_1_ok = False
#                     for voisin in graphe_tot[(num_ch_1, pos_1+1)] :
#                         if graphe_tot[(num_ch_1, pos_1+1)][voisin]["label"] == 'THS' :
#                             liaison_1_plus_1_ok = True
#                             pos = voisin
#                     
#                     liaison_3_moins_1_ok = False
#                     if pos_3-1 > 0 :
#                         for voisin in graphe_tot[(num_ch_3, pos_3-1)] :
#                             if graphe_tot[(num_ch_3, pos_3-1)][voisin]["label"] == 'TSH' :
#                                 liaison_3_moins_1_ok = True
#                     
#                     if liaison_1_1_ok and liaison_3_ok and liaison_1_2_ok and liaison_1_plus_1_ok and liaison_3_moins_1_ok :
#                         gbulge.append(fic)
                        
            
                    #print(graphe_tot[num].nodes.data())
    
    print(gnra)
    print(len(gnra))
    print(arich) 
    print(len(arich))
    print(gbulge)
    print(compteur)
    
#     for elt in gnra : 
#         if elt in GROUPE_GNRA:
#             print(elt)   
#     print("\n")
#     for elt in gnra : 
#         if elt not in GROUPE_GNRA:
#             print(elt)                
    
#test()
def retrouver_GNRA():
    with open("graphs_2.92.pickle", 'rb') as fichier_tot :
        mon_depickler = pickle.Unpickler(fichier_tot)
        graphe_tot = mon_depickler.load()
        
        
        gnra = []
        arich = []
        gbulge = []
        compteur = 0
        for fic in os.listdir(EXTENSION_PATH_TAILLE%10) :
            if "pickle" in fic and len(fic.split("_")) == 5 and "graphe_comp" not in fic :
                
                with open(EXTENSION_PATH_TAILLE%10+fic, 'rb') as fichier_graphe :
                    mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
                    graphe = mon_depickler_graphe.load()
                    occ = fic.split("_")[1]
                    chaine = fic.split("_")[2]
                    num_motif = fic.split("_")[3]
            
                    num = (occ, chaine)
                    if num == ('5J5B', 'BA') :
                        print(graphe_tot[num][1173])
                    
                    #if int(num_motif) in [44,58,62,197,251,271,137,50,328,118,134,74] :
                    compteur += 1
                    print(fic) 
               
                
                    pos_1 = graphe.nodes[1]["position"][0]
                    pos_3 = graphe.nodes[3]["position"][0]
                    
                    if graphe_tot[num].nodes[pos_1]["nt"] == 'A' and (graphe_tot[num].nodes[pos_3]["nt"] == 'A' or graphe_tot[num].nodes[pos_3]["nt"] == 'G') and graphe_tot[num].nodes[pos_3-2]["nt"] == 'G'  : #and graphe_tot[num].nodes[pos_3-4]["part"] == "Stem" and graphe_tot[num].nodes[pos_1+2]["part"] == "Stem" :  
                        
                        boucle_de_4 = False
                        if (pos_1+1, pos_3-3) in graphe_tot[num].edges() :
                            if graphe_tot[num][pos_1+1][pos_3-3]["label"] == 'CWW' :
                                boucle_de_4 = True
                        
                        pos_4_ok = False
                        pos_4 = graphe.nodes[4]["position"][0] 
                        if (graphe_tot[num].nodes[pos_4]["nt"] == 'C' or graphe_tot[num].nodes[pos_4]["nt"] == 'U') :
                            for voisin in graphe_tot[num][pos_4] :
                                if graphe_tot[num][pos_4][voisin]["label"] == 'CWW' and (graphe_tot[num].nodes[voisin]["nt"] == 'G' or graphe_tot[num].nodes[voisin]["nt"] == 'A') :
                                    pos_4_ok = True
                                    
                        liaison_ok = False
                        for voisin in graphe_tot[num][pos_1] :
                            if voisin == pos_3-2 and graphe_tot[num][pos_1][voisin]["label"] == 'THS' :
                                liaison_ok = True
                        if boucle_de_4 and pos_4_ok :
                            gnra.append(fic[8:len(fic)-7])
                        #gnra.append(fic[8:len(fic)-7])
                            
                    if graphe_tot[num].nodes[pos_1]["nt"] == 'A' and graphe_tot[num].nodes[pos_3]["nt"] == 'A' and not lien_chaines_1_3(graphe) :
    #                     liaison_1_ok = False
    #                     a_voisin_plus_1 = False
    #                     for voisin in graphe_tot[num][pos_1] :
    #                         if graphe_tot[num][pos_1][voisin]["label"] == 'THS' :
    #                             liaison_1_ok = True
    #                             print(graphe_tot[num].nodes[voisin+1]["nt"])
    #                             if graphe_tot[num].nodes[voisin+1]["nt"] == 'A' :
    #                                 a_voisin_plus_1 = True
    
    
                                
    #                     liaison_3_1_ok = False
    #                     for voisin in graphe_tot[num][pos_3-1] :
    #                         if graphe_tot[num][pos_3-1][voisin]["label"] == 'TWH' :
    #                             liaison_3_1_ok = True
                        liaison_3_1_ok = True
                                
                        liaison_stem = False
                        g_a_voisin_plus_1 = False
                        for voisin in graphe_tot[num][pos_1+1] :
                            if graphe_tot[num][pos_1+1][voisin]["label"] == 'CWW' :
                                liaison_stem = True
                                if graphe_tot[num].nodes[voisin+1]["nt"] == 'G' and graphe_tot[num].nodes[voisin+2]["nt"] == 'A' :
                                    g_a_voisin_plus_1 = True
                                
                        
                        if liaison_3_1_ok and liaison_stem and g_a_voisin_plus_1 :
                            arich.append(fic[8:len(fic)-7])
                            
                            
                    if graphe_tot[num].nodes[pos_1]["nt"] == 'A' and graphe_tot[num].nodes[pos_3]["nt"] == 'A' :
                        liaison_1_1_ok = False
                        pos = -1
                        for voisin in graphe_tot[num][pos_1] :
                            if graphe_tot[num][pos_1][voisin]["label"] == 'THH' :
                                liaison_1_1_ok = True
                                pos = voisin
    
                        liaison_1_2_ok = False
                        if pos != -1 :
                            for voisin in graphe_tot[num][pos+1] :
                                if graphe_tot[num][pos+1][voisin]["label"] == 'CSH' :
                                    liaison_1_2_ok = True
    
                        liaison_3_ok = False
                        for voisin in graphe_tot[num][pos_3] :
                            if graphe_tot[num][pos_3][voisin]["label"] == 'THW' :
                                liaison_3_ok = True
                        
                        liaison_1_plus_1_ok = False
                        for voisin in graphe_tot[num][pos_1+1] :
                            if graphe_tot[num][pos_1+1][voisin]["label"] == 'THS' :
                                liaison_1_plus_1_ok = True
                                pos = voisin
                        
                        liaison_3_moins_1_ok = False
                        for voisin in graphe_tot[num][pos_3-1] :
                            if graphe_tot[num][pos_3-1][voisin]["label"] == 'TSH' :
                                liaison_3_moins_1_ok = True
                        
                        if liaison_1_1_ok and liaison_3_ok and liaison_1_2_ok :
                            gbulge.append(fic[8:len(fic)-7])
                            
                
                        #print(graphe_tot[num].nodes.data())
        
        print(gnra)
        print(len(gnra))
        print(arich) 
        print(len(arich))
        print(gbulge)
        print(compteur)
        
        for elt in gnra : 
            if elt in GROUPE_GNRA:
                print(elt)   
        print("\n")
        for elt in gnra : 
            if elt not in GROUPE_GNRA:
                print(elt)        
                
if __name__ == '__main__':
    version_new_data()  
    