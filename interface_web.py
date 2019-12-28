'''
Created on 18 mars 2019

@author: coline

Creation de tous les scripts python permettant les ajouts a la base de donnees
de l'interface Web Django 
(version CaRNAval)
'''

import os
import pickle
from recup_data.constantes import EXTENSION_PATH_TAILLE,\
    GROUPES_TOUTES_ARETES_MAX_4_10_07, EXTENSION_PATH, GROUPE_GNRA, GROUPE_ARICH,\
    GROUPE_GBULGE, GROUPE_GNRA_ETENDU, GROUPE_ARICH_ETENDU, PATH_MMCIF
from recup_data.calcul_sim import calcul_sim_aretes_avec_coeff
import numpy as np

''' renvoie un dictionnaire contenant les coordonnees des noeuds sur Gephi (a partir du fichier_gexf) 
pas vraiment utile en fait '''
def recup_coordonnees_gephi(fichier_gexf):
    ligne = fichier_gexf.readline()
    dans_noeud = False
    dict_coordonnees = {}
    while ligne != "" :
        if "<node" in ligne and "/node" not in ligne and "nodes" not in ligne :  #debut d un noeud
            dans_noeud = True
            label = ligne.split("<")[1].split(" ")[2].split('"')[1]
        if dans_noeud and "viz:position" in ligne :
            coordonnee_x = float(ligne.split("<")[1].split(" ")[1].split('"')[1])
            coordonnee_y = float(ligne.split("<")[1].split(" ")[2].split('"')[1])
            
            dict_coordonnees.update({label : (coordonnee_x, coordonnee_y)})
        ligne = fichier_gexf.readline()
    
    return dict_coordonnees       
            

''' renvoie l'ensemble des positions dans le graphe global des nucleotides stockes dans le graphe d'extension associe '''
def recup_positions_extensions(graphe_extension) :
    positions = []

    for noeud, data in graphe_extension.nodes(data=True) :
        if data["position"][0] != -1 :
            for i in range(data["position"][0], data["position"][1]+1) :
                positions.append(i)


    return positions

''' renvoie les lignes html permettant d'afficher une fenetre jmol avec la structure en 3D correspondant au graphe d'extension passe en parametre '''
def motif3D_html(cle, num_motif, num_occ, graphe_extension, graphes) :
    head = """

      <input type="button" id="neighborhood" value="Show neighborhood">
      <label><input type="checkbox" id="showNtNums">Nucleotide numbers</label>
      <input type="button" id="stereo" value="Stereo">
      <div>
        <!-- unit ID list -->
        <input type='radio' name="group" id='s1' class='jmolInline' data-coord='%s'><label for='s1'>%s</label>
      </div>
    </body>
    </html>
    """
    
    nt_data = '%s|1|%s|%s|%s' #% (pdb, model=1, chain, nt, pos)
    pdb = cle[0]
    chain = cle[1]

    data_nodes = recup_positions_extensions(graphe_extension)
    data_coord = ','.join(nt_data % (pdb, chain, graphes[cle].nodes[x]['nt'], graphes[cle].nodes[x]['fr3d']) for x in data_nodes)
    print(data_coord)
    
    return head%(data_coord, '%s_%s_%s_%s' % (pdb, chain, num_motif, num_occ))

''' renvoie les coordonnees des nucleotides de la structure en 3D correspondant au graphe d'extension passe en parametre,
pouvant etre utilisee pour creer une fenetre jmol '''
def motif3D_html_seul(cle, graphe_extension, graphes):
    nt_data = '%s|1|%s|%s|%s' #% (pdb, model=1, chain, nt, pos)
    pdb = cle[0]
    chain = cle[1]
    data_nodes = recup_positions_extensions(graphe_extension)
    data_coord = ','.join(nt_data % (pdb, chain, graphes[cle].nodes[x]['nt'], graphes[cle].nodes[x]['fr3d']) for x in data_nodes)
    return data_coord

''' creation d'un fichier pickle dans lequel sont stockees les informations des fichiers PDB contenant des motifs A-minor dans CaRNAval '''
def ajout_infos_molecules_pickle():
    tab_infos = [{'num_PDB' : '2XD0', 'num_ch' : 'V', 'type_mol' : 'ARN non codant (toxine)', 'organisme' : 'Pectobacterium atrosepticum', 'libelle_pdb' : 'A processed non-coding RNA regulates a bacterial antiviral system'},
                 {'num_PDB' : '1FJG', 'num_ch' : 'A', 'type_mol' : 'ARNr 16S', 'organisme' : 'Thermus  thermophilus', 'libelle_pdb' : 'STRUCTURE OF THE THERMUS THERMOPHILUS 30S RIBOSOMAL SUBUNIT IN COMPLEX WITH THE ANTIBIOTICS STREPTOMYCIN, SPECTINOMYCIN, AND PAROMOMYCIN'},
                 {'num_PDB' : '5J5B', 'num_ch' : 'BA', 'type_mol' : 'ARNr 16S', 'organisme' : 'Escherichia coli', 'libelle_pdb' : 'Structure of the WT E coli ribosome bound to tetracycline'},
                 {'num_PDB' : '4V9F', 'num_ch' : '0', 'type_mol' : 'ARNr 23S', 'organisme' : 'Haloarcula marismortui', 'libelle_pdb' : 'The re-refined crystal structure of the Haloarcula marismortui large ribosomal subunit at 2.4 Angstrom resolution: more complete structure of the L7/L12 and L1 stalk, L5 and LX proteins'},
                 {'num_PDB' : '5J7L', 'num_ch' : 'DA', 'type_mol' : 'ARNr 23S', 'organisme' : 'Escherichia coli', 'libelle_pdb' : 'Structure of the 70S E coli ribosome with the U1052G mutation in the 16S rRNA bound to tetracycline'},
                 {'num_PDB' : '5DM6', 'num_ch' : 'X', 'type_mol' : 'ARNr 23S', 'organisme' : 'Deinococcus radiodurans', 'libelle_pdb' : 'Crystal structure of the 50S ribosomal subunit from Deinococcus radiodurans'},
                 {'num_PDB' : '5FDU', 'num_ch' : '1A', 'type_mol' : 'ARNr 23S', 'organisme' : 'Thermus thermophilus', 'libelle_pdb' : 'Crystal structure of the Metalnikowin I antimicrobial peptide bound to the Thermus thermophilus 70S ribosome'},
                 {'num_PDB' : '4V88', 'num_ch' : 'A5', 'type_mol' : 'ARNr 25S', 'organisme' : 'Saccharomyces cerevisiae', 'libelle_pdb' : 'The structure of the eukaryotic ribosome at 3.0 A resolution'},
                 {'num_PDB' : '3JCS', 'num_ch' : '1', 'type_mol' : 'ARNr 26S alpha', 'organisme' : 'Leishmania donovani', 'libelle_pdb' : '2.8 Angstrom cryo-EM structure of the large ribosomal subunit from the eukaryotic parasite Leishmania'},
                 {'num_PDB' : '3JCS', 'num_ch' : '2', 'type_mol' : 'ARNr 26S delta', 'organisme' : 'Leishmania donovani', 'libelle_pdb' : '2.8 Angstrom cryo-EM structure of the large ribosomal subunit from the eukaryotic parasite Leishmania'},
                 {'num_PDB' : '4FAU', 'num_ch' : 'A', 'type_mol' : 'Intron', 'organisme' : 'Oceanobacillus iheyensis', 'libelle_pdb' : "Structure of Oceanobacillus iheyensis group II intron in the presence of Li+, Mg2+ and 5'-exon"},
                 {'num_PDB' : '1U9S', 'num_ch' : 'A', 'type_mol' : 'Ribonucléase P', 'organisme' : 'Inconnu', 'libelle_pdb' : "Crystal structure of the specificity domain of Ribonuclease P of the A-type"},
                 {'num_PDB' : '4YAZ', 'num_ch' : 'R', 'type_mol' : 'Riboswitch', 'organisme' : 'Geobacter', 'libelle_pdb' : "3',3'-cGAMP riboswitch bound with 3',3'-cGAMP"},
                 {'num_PDB' : '3UCZ', 'num_ch' : 'R', 'type_mol' : 'Riboswitch', 'organisme' : 'Homo sapiens', 'libelle_pdb' : "The c-di-GMP-I riboswitch bound to GpG"},
                 {'num_PDB' : '5FJC', 'num_ch' : 'A', 'type_mol' : 'Riboswitch', 'organisme' : 'Caldanaerobacter subterraneurs subsp. Tengcongensis', 'libelle_pdb' : "SAM-I riboswitch bearing the H. marismortui Kt-7 variant C-2bU"},
                 {'num_PDB' : '4L81', 'num_ch' : 'A', 'type_mol' : 'Riboswitch (aptamère, variant)', 'organisme' : 'Inconnu', 'libelle_pdb' : "Structure of the SAM-I/IV riboswitch (env87(deltaU92, deltaG93))"},
                 {'num_PDB' : '4PRF', 'num_ch' : 'B', 'type_mol' : 'Ribozyme', 'organisme' : 'Hepatitis Delta virus', 'libelle_pdb' : "A Second Look at the HDV Ribozyme Structure and Dynamics"},
                 {'num_PDB' : '1GID', 'num_ch' : 'B', 'type_mol' : 'Ribozyme (domaine P4,P6)', 'organisme' : 'synthétique', 'libelle_pdb' : "CRYSTAL STRUCTURE OF A GROUP I RIBOZYME DOMAIN: PRINCIPLES OF RNA PACKING"},
                 {'num_PDB' : '4V88', 'num_ch' : 'A6', 'type_mol' : 'ARNr 18S', 'organisme' : 'Saccharomyces cerevisiae', 'libelle_pdb' : "The structure of the eukaryotic ribosome at 3.0 A resolution"}]
    with open("infos_molecules.pickle", 'wb') as fichier_pickle :
        mon_pickler = pickle.Pickler(fichier_pickle)
        mon_pickler.dump(tab_infos)            
                 
''' dans le script script_python_ajout_base.py, ajout des lignes de code permettant l'ajout d'un graphe d'extension a la base de donnees '''                 
def ajout_base_extension(cle):
    with open("script_python_ajout_base.py", 'a') as fichier_py :
        with open("infos_molecules.pickle", 'rb') as fichier_pickle :
            mon_depickler = pickle.Unpickler(fichier_pickle)
            tab_infos = mon_depickler.load()  
            for elt in tab_infos :
                if elt["num_PDB"] == cle[0] and elt["num_ch"] == cle[1] :
                    nom_cle = cle[0] + "_" + cle[1] + "_" + str(cle[2]) + "_" + str(cle[3])
                    fichier_py.write("e = Extension(nom_molecule='" + nom_cle + "', libelle_PDB='"+elt["libelle_pdb"]+"', type_molecule='"+ elt["type_mol"]+ "', organisme='"+ elt["organisme"]+ "', num_PDB='" + elt["num_PDB"] + "', num_ch='"+elt["num_ch"]+ "', num_motif='"+cle[2]+"', num_occ='"+cle[3]+ "', url_PDB='https://www.rcsb.org/structure/"+elt["num_PDB"].lower()+ "')\n")
                    fichier_py.write("e.save()\n")

''' dans le script script_python_alter_extension.py, ajout des lignes de code permettant l'ajout dans la base de donnees 
d'un attribut au graphe d'extension passe en parametre : numero de cluster dans un certain type de clustering '''                 
def alter_base_extension_clustering(cle, num_groupe, type_clustering):
    with open("script_python_alter_extension.py", 'a') as fichier_py : 
        fichier_py.write("e = Extension.objects.get(nom_molecule='"+cle+ "')\n")
        fichier_py.write("e.%s=%d\n"%(type_clustering, num_groupe))
        fichier_py.write("e.save()\n")

''' dans le script script_python_ajout_liens_par_taille.py, ajout des lignes de code permettant l'ajout du lien entre cle_1 et cle_2
avec toutes les similarites entre ces deux graphes d'extension par taille d'extension'''
def ajout_base_sim_par_taille(cle_1, cle_2, sim, taille):
    with open("script_python_ajout_liens_par_taille.py", 'a') as fichier_py :
        fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+cle_1+ "')\n")
        fichier_py.write("e2 = Extension.objects.get(nom_molecule='"+cle_2+ "')\n")
        
        fichier_py.write("sim = Sim_par_taille(name='%s', taille_extension=%d, sim=%1.2f)\n"%((cle_1 + "_" + cle_2), taille, sim))
        fichier_py.write("sim.save()\n")

        
        fichier_py.write("try :\n\tl = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)\n") 
        fichier_py.write("\tl.sim_par_taille.add(sim)\n") 
        fichier_py.write("except django.core.exceptions.ObjectDoesNotExist :\n")
        fichier_py.write("\t")
        #fichier_py.write("try :\n")
        #fichier_py.write("\t\t")
        fichier_py.write("l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)\n")
        fichier_py.write("\t")
        fichier_py.write("l.sim_par_taille.add(sim)\n") 
        
        fichier_py.write("l.save()\n")


''' dans le script script_python_ajout_liens_par_taille.py, ajout des lignes de code permettant la modification de la similarite associee
a une taille d'extension entre les graphes d'extension cle_1 et cle_2 '''
def alter_base_sim_par_taille(cle_1, cle_2, sim, taille):
    with open("script_python_alter_liens_par_taille.py", 'a') as fichier_py :
        fichier_py.write("sim = Sim_par_taille.objects.get(name='%s', taille_extension=%d)\n"%((cle_1 + "_" + cle_2), taille))
        fichier_py.write("sim.sim = %1.2f\n"%sim)
        fichier_py.write("sim.save()\n")
        
''' dans le script script_python_alter_liens_rmsd.py, ajout des lignes de code permettant l'ajout d'un attribut au lien entre les graphes
cle_1 et cle_2 : rmsd normalisee taile 4 '''
def alter_base_rmsd(cle_1, cle_2, rmsd):
    with open("script_python_alter_liens_rmsd.py", 'a') as fichier_py :
        fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+cle_1+ "')\n")
        fichier_py.write("e2 = Extension.objects.get(nom_molecule='"+cle_2+ "')\n")
        fichier_py.write("try :\n\tl = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)\n") 
        fichier_py.write("\tl.rmsd_normalise_taille_4=%.2f\n"%round(rmsd,2)) 
        fichier_py.write("except django.core.exceptions.ObjectDoesNotExist :\n")
        fichier_py.write("\t")
        #fichier_py.write("try :\n")
        #fichier_py.write("\t\t")
        fichier_py.write("l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)\n")
        fichier_py.write("\t")
        fichier_py.write("l.rmsd_normalise_taille_4=%.2f\n"%round(rmsd,2)) 
        fichier_py.write("l.save()\n")

''' dans le script script_python_ajout_clustering_kmeans.py, ajout des lignes de code permettant l'ajout d'un attribut a tous les graphes
d'extension de la base de donnes : le numero du cluster dans clustering kmeans selon la rmsd'''
def ajout_clustering_perez(dico_extension):       
    with open("script_python_ajout_clustering_kmeans.py", 'w') as fichier_py :
        fichier_py.write("from aminor.models import Extension\n")
       
        for cle in dico_extension.keys() :
            fichier_py.write("e = Extension.objects.get(nom_molecule='"+cle+ "')\n")
            fichier_py.write("e.clustering_kmeans_rmsd = %s\n"%dico_extension[cle])
            fichier_py.write("e.save()\n")
            
''' ne sais plus trop ce que c'est non_can_2_avec_4_6 '''               
def ajout_base_liens(cle_1, cle_2, sim):
    with open("script_python_ajout_liens.py", 'a') as fichier_py :
        #fichier_py.write("from aminor.models import Extension\nfrom aminor.models import Lien\nimport django.core.exceptions\n")
        fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+cle_1+ "')\n")
        fichier_py.write("e2 = Extension.objects.get(nom_molecule='"+cle_2+ "')\n")
        fichier_py.write("try :\n\tl = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)\n") 
        fichier_py.write("\tl.sim_non_can_2_avec_4_6=%1.2f\n"%round(sim[0],2)) 
        fichier_py.write("\t")
        fichier_py.write("l.sim_k_non_can_2_avec_4_6=%d\n"%sim[1])
        fichier_py.write("except django.core.exceptions.ObjectDoesNotExist :\n")
        fichier_py.write("\t")
        #fichier_py.write("try :\n")
        #fichier_py.write("\t\t")
        fichier_py.write("l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)\n")
        fichier_py.write("\t")
        fichier_py.write("l.sim_non_can_2_avec_4_6=%1.2f\n"%round(sim[0],2)) 
        fichier_py.write("\t")
        fichier_py.write("l.sim_k_non_can_2_avec_4_6=%d\n"%sim[1])
        #fichier_py.write("\t")
        #fichier_py.write("except django.core.exceptions.ObjectDoesNotExist :\n")
        #fichier_py.write("\t\t")
        #fichier_py.write("g = Lien(num_ext_1 = e1, num_ext_2 = e2, sim_new="+str(round(sim,2))+")\n")
        #fichier_py.write("\t\t")
        fichier_py.write("l.save()\n")

''' dans le script script_python_ajout_base_graphes_taille.py, ajout des lignes de code permettant l'ajout d'un nouveau graph contenant les liens
avec une sim_non_can_2_avec_4_6 sup a 0.7 '''
def ajout_base_graphes(graphe, i, taille):
        with open("script_python_ajout_base_graphes_taille.py", 'a') as fichier_py :
            #fichier_py.write("from aminor.models import Extension\nfrom aminor.models import Lien\nfrom aminor.models import Graph_version_2\nimport django.core.exceptions\n\n")
#             fichier_py.write("g = Graph_version_2.objects.get(graph_nom='toutes_aretes_coeffn1_a1_c1_taille_max_non_can_2_2_"+str(i)+"', graph_poids="+str(i)+")\n")
#             fichier_py.write("g.delete()\n")
            fichier_py.write("g = Graph_version_2(graph_nom='toutes_aretes_coeffn1_a1_c1_taille_max_non_can_2_2_avec_4_6_"+str(i)+"', graph_poids="+str(i)+")\n")
            fichier_py.write("g.save() \n")
            noeuds = list(graphe.nodes())
            print(noeuds)
            for j in range(len(noeuds)) :
                for k in range(j+1, len(noeuds)) :
                    #if graphe.edges[noeuds[j], noeuds[k]]["poids"] >= i :
                        cle_1 = graphe.nodes[noeuds[j]]["nom"]
                        cle_2 = graphe.nodes[noeuds[k]]["nom"]
                        fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+cle_1+ "')\n")
                        fichier_py.write("e2 = Extension.objects.get(nom_molecule='"+cle_2+ "')\n")
                        fichier_py.write("try : \n\t l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)\n")
                        fichier_py.write("except django.core.exceptions.ObjectDoesNotExist :\n")
                        fichier_py.write("\t")
                        fichier_py.write("l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)\n")
                        fichier_py.write("if l.sim_non_can_2_avec_4_6 > 0.7 :\n")
                        fichier_py.write("\tg.graphe.add(l) \n")


''' dans le script script_python_ajout_base_graphes_voisinage.py, ajout des lignes de code permettant de
creer un graphe voisinage contenant les elements d'un cluster et leur voisinage (avec des aretes sup a un seuil sans doute)'''
def ajout_base_graphes_voisinage(composante_cluster, composante_voisinage, i, compteur):
    print(composante_cluster)
    print(composante_voisinage)
    composante_totale = list(set(composante_cluster).union(composante_voisinage))
    print(composante_totale)
    with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim_coeff_n1_a1_c1.pickle", 'rb') as fichier_graphe :
        mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
        graphe_complet = mon_depickler_graphe.load()
        with open("script_python_ajout_base_graphes_voisinage.py", 'a') as fichier_py :
            fichier_py.write("g = Graph_voisinage(graph_nom='toutes_aretes_coeffn1_a1_c1_"+str(i)+"_"+str(compteur) + "_voisinage', graph_poids= "+str(i)+",graph_num="+str(compteur)+")\n")
            fichier_py.write("g.save() \n")
            
            for elt in composante_cluster :
                cle = graphe_complet.nodes[elt]["nom"]
                fichier_py.write("e = Extension.objects.get(nom_molecule='"+cle+ "')\n")
                fichier_py.write("g.graph_cluster.add(e) \n")
            
            for elt in composante_voisinage :
                cle = graphe_complet.nodes[elt]["nom"]
                fichier_py.write("e = Extension.objects.get(nom_molecule='"+cle+ "')\n")
                fichier_py.write("g.graph_voisinage.add(e) \n")
            
            for j in range(len(composante_totale)) :
                for k in range(j+1, len(composante_totale)) :
                    cle_1 = graphe_complet.nodes[composante_totale[j]]["nom"]
                    cle_2 = graphe_complet.nodes[composante_totale[k]]["nom"]
                    fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+cle_1+ "')\n")
                    fichier_py.write("e2 = Extension.objects.get(nom_molecule='"+cle_2+ "')\n")
                    fichier_py.write("try :\n\t l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)\n")
                    fichier_py.write("except django.core.exceptions.ObjectDoesNotExist :\n")
                    fichier_py.write("\t")
                    fichier_py.write("l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)\n")
                    fichier_py.write("g.graphe.add(l) \n")
    
''' dans le script script_python_ajout_homologues.py , ajout des lignes de code permettant d'ajouter un attribut a tous les graphes d'extension
: groupe d'homologues '''                    
def ajout_homologues():
    homologues = [['4V9F_0_62_12', '5J7L_DA_62_5', '5DM6_X_328_2', '5FDU_1A_62_14'],
                  ['4V9F_0_30_4', '5J7L_DA_30_15', '5FDU_1A_30_17'],
                  ['4V9F_0_48_13', '5J7L_DA_48_20', '5DM6_X_48_28', '5FDU_1A_48_25'],
                  ['4V9F_0_207_3', '5J7L_DA_272_2', '5FDU_1A_272_1'],
                  ['4V9F_0_25_56', '5J7L_DA_25_12', '5FDU_1A_25_68', '4V88_A5_25_30'],
                  ['4V9F_0_137_5', '5J7L_DA_48_1', '5FDU_1A_137_6', '5DM6_X_48_10', '4V88_A5_48_3'],
                  ['4V9F_0_134_5', '5FDU_1A_74_7', '5DM6_X_127_7'],
                  ['4V9F_0_48_21', '5J7L_DA_197_4', '5FDU_1A_197_3', '5DM6_X_48_9'],
                  ['4V9F_0_127_6', '5J7L_DA_134_1', '5FDU_1A_134_3', '5DM6_X_134_2'],
                  ['4V9F_0_287_2', '5DM6_X_25_15'],
                  ['5J7L_DA_25_10', '5DM6_X_25_34', '4V88_A5_25_47', '5FDU_1A_25_78'],
                  ['1FJG_A_48_17', '5J5B_BA_48_14'],
                  ['1FJG_A_48_8', '5J5B_BA_48_23'],
                  ['1FJG_A_294_1', '5J5B_BA_294_2'],
                  ['1FJG_A_58_23', '5J5B_BA_58_3'],
                  ['1FJG_A_138_3', '5J5B_BA_138_2'],
                  ['5J5B_BA_48_7', '4V88_A6_48_12'],
                  ['4YAZ_R_36_25', '3UCZ_R_62_15', '5FJC_A_138_1', '4L81_A_25_77']]
    with open("script_python_ajout_homologues.py", 'a') as fichier_py :
        compteur = 1
        for groupe_homologues in homologues :
            for elt in groupe_homologues :
                fichier_py.write("e = Extension.objects.get(nom_molecule='%s')\n"%(elt)) 
                fichier_py.write("e.groupe_homologues = %d\n"%(compteur)) 
                fichier_py.write("e.save() \n")
            compteur += 1    

''' dans le script script_python_ajout_image_extension.py, ajout des lignes de code permettant d'ajouter un attribut a tous les graphes d'extension de la base de donnees :
l'image png du graphe '''
def ajout_image_extension():
    with open("script_python_ajout_image_extension.py", 'a') as fichier_py :
        for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension_avec_liaisons_b53_qui_manquent") :
            if "pickle" in fic :
                fichier_py.write("e = Extension.objects.get(nom_molecule='%s')\n"%(fic[8:len(fic)-7]))
                fichier_py.write("e.url_image = '%s'\n"%(fic[:len(fic)-7]+'.png'))
                fichier_py.write("e.save() \n")
                
''' dans le script script_ajout_base_graphe_version_3.py, ajout des lignes de code permettant de creer une instance de graph_version_3
contenant certains liens, je ne sais pas trop lesquels... '''  
def ajout_base_graphe_version_3(dico_sim, seuil, taille_ext, type_sim ):
    with open("script_ajout_base_graphe_version_3.py", 'a') as fichier_py : 
        #fichier_py.write('from aminor.models import Extension \nfrom aminor.models import Lien\nfrom aminor.models import Graph_version_3\nimport django.core.exceptions\n')
        fichier_py.write('g = Graph_version_3(graph_nom="%s", graph_poids=%1.2f, graph_taille_ext="%s", graph_type_sim="%s")\n'%('graphe_toutes_aretes_%s_%s_%s'%(type_sim, seuil, taille_ext), seuil, taille_ext, type_sim))
        fichier_py.write('g.save()\n')
        
        for cle in dico_sim.keys() :
            print(dico_sim[cle])
            if dico_sim[cle] != None and ((dico_sim[cle] > seuil and type_sim == 'sim') or ( dico_sim[cle] < seuil and type_sim == 'rmsd')):
                fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+cle[0]+ "')\n")
                fichier_py.write("e2 = Extension.objects.get(nom_molecule='"+cle[1]+ "')\n")
                fichier_py.write("try : \n\t l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)\n")
                fichier_py.write("except django.core.exceptions.ObjectDoesNotExist :\n")
                fichier_py.write("\t")
                fichier_py.write("l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)\n")
                fichier_py.write("g.graphe.add(l) \n")
                 
        
if __name__ == '__main__':
    #ajout_infos_molecules_pickle()
#     with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
#         mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
#         graphes = mon_depickler_graphes.load()
#         for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension") :
#             if "pickle" in fic :
#                 with open("Extensions/Metrique_toutes_aretes/graphes_extension/"+fic, 'rb') as fichier_pickle :
#                     mon_depickler = pickle.Unpickler(fichier_pickle)
#                     graphe = mon_depickler.load()
#                       
#                     print(fic)
#                     cle = (fic.split("_")[1], fic.split("_")[2])
#                     print(cle)
#                     debut_head = """<html>
#                             <head>
#                               <title>Radio buttons</title>
#                               <script src="https://ajax.googleapis.com/ajax/libs/jquery/1.11.2/jquery.min.js"></script>
#                               
#                               {% load compress %}
#                               {% load static %}
#                               
#                               {% compress js %}
#                               <script src="{% static 'jmolTools-master/docs/jsmol/JSmol.min.nojq.js' %}"></script>
#                               {% endcompress %}
#                               {% compress js %}
#                               <script src="{% static 'jmolTools-master/jquery.jmolTools.js' %}"></script>
#                               {% endcompress %}
#                             
#                             </head>
#                             
#                             <body>
#                             
#                               <div>
#                                 {% compress js %}
#                                 <script src="{% static 'jmolTools-master/jsmol_initializer.js' %}"></script>
#                                 {% endcompress %}
#                               </div>"""
#                     commandes = motif3D_html(cle, fic.split("_")[3], fic.split("_")[4][:len(fic.split("_")[4])-7], graphe, graphes)
#                     with open('jmol_'+fic[:len(fic)-7]+'.html', 'w') as fichier :
#                         fichier.write(debut_head+commandes)
#                         
#                 pos = recup_positions_extensions(graphe)
#                 print(fic)
#                 print(pos)
                

#         for fic in os.listdir("Extensions/Metrique_toutes_aretes/graphes_extension") :
#             if "pickle" in fic :
#                     print(fic)
#                     cle = (fic.split("_")[1], fic.split("_")[2], fic.split("_")[3], fic.split("_")[4][:len(fic.split("_")[4])-7])
#                     print(cle)
#                     
#                     ajout_base(cle)

            
#             with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"graphe_complet_pondere_sim_toutes_aretes_coeff_all1_max_taille_taille_max.pickle", 'rb') as fichier_graphe :
#                 mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#                 graphe_complet = mon_depickler_graphe.load()  
#                 couples_deja_mis = []
#                 i = 0.7
#                 with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+ "composantes_connexes_extensions_toutes_aretes_coeff_all1_max_taille_taille_max_"+str(i)+".pickle", 'rb') as fichier_pickle :
#                         mon_depickler = pickle.Unpickler(fichier_pickle)
#                         composantes = mon_depickler.load()
#                         for composante in composantes :
#                             if len(composante) > 10:
#                                 for j in range(len(composante)) :
#                                     for k in range(j+1, len(composante)) :
#                                         cle_1 = graphe_complet.nodes[composante[j]]["nom"]
#                                         cle_2 = graphe_complet.nodes[composante[k]]["nom"]
#                                         print(cle_1)
#                                         print(cle_2)
#                                         if graphe_complet.edges[composante[j], composante[k]]["poids"] > i and (cle_1, cle_2) not in couples_deja_mis and (cle_2, cle_1) not in couples_deja_mis :
#                                             ajout_base_liens(cle_1, cle_2, graphe_complet.edges[composante[j], composante[k]]["poids"])
#                                             couples_deja_mis.append((cle_1, cle_2))
#                             print(composante)
# ## ajout de tous les liens et de la sim associees (il faut penser à changer le nom de la sim)
#                 with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"sim_extensions_toutes_aretes_coeff_all1_max_taille_taille_max_avec_val_k.pickle", 'rb') as fichier_graphe :
#                     mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#                     dico_sim = mon_depickler_graphe.load()  
#                     couples_deja_mis = []
#                     for cle in dico_sim.keys() :
#                         cle_1 = cle[0].split("_")[0] + "_" + cle[0].split("_")[1] + "_" + cle[0].split("_")[2] + "_" + cle[0].split("_")[3]
#                         cle_2 = cle[1].split("_")[0] + "_" + cle[1].split("_")[1] + "_" + cle[1].split("_")[2] + "_" + cle[1].split("_")[3]
#                         print(cle_1)
#                         print(cle_2)
#                         
#                         with open(EXTENSION_PATH_TAILLE%dico_sim[cle][1]+"fichier_"+cle_1+".pickle", 'rb') as fichier_1 :
#                             mon_depickler_1 = pickle.Unpickler(fichier_1)
#                             graphe1 = mon_depickler_1.load()
#                             
#                             with open(EXTENSION_PATH_TAILLE%dico_sim[cle][1]+"fichier_"+cle_2+".pickle", 'rb') as fichier_2 :
#                                 mon_depickler_2 = pickle.Unpickler(fichier_2)
#                                 graphe2 = mon_depickler_2.load()
#                                 
#                                 with open(EXTENSION_PATH%dico_sim[cle][1]+"dico_comp_complet_metrique_toutes_aretes_coeff_all1_taille_%d.pickle"%dico_sim[cle][1], 'rb') as fichier_graphe_commun :
#                                     mon_depickler_commun = pickle.Unpickler(fichier_graphe_commun)
#                                     graphe_commun = mon_depickler_commun.load()
#                                     
#                                     for elt in graphe_commun.keys() :
#                                         if cle[0] == elt[0][8:] and cle[1] == elt[1][8:] :
#                                             print(cle)
#                                             print(elt)
#                                             sim = calcul_sim_aretes_avec_coeff(graphe1, graphe2, graphe_commun[elt], cle, 1, 1, 1)
#                                             sim_avec_k = (sim, dico_sim[cle][1])
#                                             if (cle_1, cle_2) not in couples_deja_mis and (cle_2, cle_1) not in couples_deja_mis :
#                                                 ajout_base_liens(cle_1, cle_2, sim_avec_k)
#                                                 couples_deja_mis.append((cle_1, cle_2))
    
#                 for i in range(1,5) :
#                     with open(EXTENSION_PATH%i+"sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%i, 'rb') as fichier_graphe :  
#                         mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#                         dico_sim = mon_depickler_graphe.load()  
#                         couples_deja_mis = []
#                         for cle in dico_sim.keys() :
#                             cle_1 = cle[0].split("_")[0] + "_" + cle[0].split("_")[1] + "_" + cle[0].split("_")[2] + "_" + cle[0].split("_")[3]
#                             cle_2 = cle[1].split("_")[0] + "_" + cle[1].split("_")[1] + "_" + cle[1].split("_")[2] + "_" + cle[1].split("_")[3]
#                             print(cle_1)
#                             print(cle_2)
#                             if (cle_1, cle_2) not in couples_deja_mis and (cle_2, cle_1) not in couples_deja_mis :
#                                 #ajout_base_sim_par_taille(cle_1, cle_2, dico_sim[cle], i)
#                                 alter_base_sim_par_taille(cle_1, cle_2, dico_sim[cle], i)
#                                 couples_deja_mis.append((cle_1, cle_2))
#             with open(EXTENSION_PATH%"taille_max/result_k_max_4_10_toutes_aretes"+"graphe_complet_pondere_sim_toutes_aretes_coeff_all1_max_taille_taille_max.pickle", 'rb') as fichier_graphe :
# #             #for i in range(1,11) : 
#             #    with open(EXTENSION_PATH%i+"graphe_complet_pondere_sim_toutes_aretes_coeff_all1_taille_%d.pickle"%i, 'rb') as fichier_graphe : 
#                     mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#                     graphe_complet = mon_depickler_graphe.load()  
#                     ajout_base_graphes(graphe_complet, 0.7, 0)
                #couples_deja_mis = []
#                 i = 0.7
#                 #while i < 1.0 :
#                 with open(EXTENSION_PATH%"taille_max/result_k_max_4_8_toutes_aretes"+ "composantes_connexes_extensions_toutes_aretes_coeff_all1_max_taille_taille_max_"+str(i)+".pickle", 'rb') as fichier_pickle :
#                         mon_depickler = pickle.Unpickler(fichier_pickle)
#                         composantes = mon_depickler.load()
#                         compteur = 1
#                         composantes_regroupees = []
#                         for composante in composantes :
#                             if len(composante) > 10 :
#                                 composantes_regroupees.extend(composante)
#                                 print(composante)
#                                 compteur += 1
#                         print(composantes_regroupees)
                        #print(composante)
                        #i = i+0.1            
                
        #ajout_homologues()
        #ajout_image_extension()
        
       ## ajout de tous les liens 
#         with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim_coeff_n1_a1_c1.pickle", 'rb') as fichier_graphe :
#             mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#             graphe_complet = mon_depickler_graphe.load()     
#             
#             for u,v, sim in graphe_complet.edges(data="poids") :
#                 ajout_base_liens(graphe_complet.nodes[u]["nom"], graphe_complet.nodes[v]["nom"], sim)
        
        
        ##voisinage du cluster 0.6 2 metrique toutes aretes coeff all1
        
#         with open("Extensions/Metrique_toutes_aretes/graphe_complet_pondere_sim_coeff_n1_a1_c1.pickle", 'rb') as fichier_graphe :
#             mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#             graphe_complet = mon_depickler_graphe.load()  
#             
#             liste_voisinage = []
#             liste_cluster = []
#             for noeud, data in graphe_complet.nodes(data=True) :
#                 if data["nom"] in ['3JCS_1_48_18', '4V88_A5_48_3', '5DM6_X_48_10', '5J7L_DA_48_1', '5FDU_1A_137_6'] :
#                     liste_voisinage.append(noeud)
#                 if data["nom"] in ['1FJG_A_271_1', '1FJG_A_109_6', '4V9F_0_48_16', '4V9F_0_30_23'] :
#                     liste_cluster.append(noeud)
#                     
#                     
#             ajout_base_graphes_voisinage(liste_cluster, liste_voisinage, 0.6, 2)
            
            
#     with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
#         mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
#         graphes = mon_depickler_graphes.load()
#         
#         with open(EXTENSION_PATH_TAILLE%4+"fichier_5DM6_X_25_15.pickle", 'rb') as fichier_graphe :
#             mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#             graphe = mon_depickler_graphe.load()
#             
#             pos = recup_positions_extensions(graphe)
#             #pos_motif = pos[:5]
#             data_script = "select ext_5dm6, 5dm6 and chain X and (resi "
#             num = "%s"
#             data_coord = " or resi ".join(num%(graphes[('5DM6','X')].nodes[x]['fr3d']) for x in pos)
#             
#             print(data_script + data_coord + ')')
            
        
#         for fic in os.listdir("/home/coline/Gephi/fichiers") :
#             if ".gexf" in fic :
#                 with open("/home/coline/Gephi/fichiers/"+fic, 'r') as fichier_gexf :
#                     dict_coordonnees = recup_coordonnees_gephi(fichier_gexf)
#                      
#                     with open("script_modif_extension_pour_coordonnees.py", 'w') as fichier_py :
#                         for cle in dict_coordonnees.keys() :
#                             fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+cle+ "')\n")
#                             fichier_py.write("e1.coordonnee_x_groupe = %f\n"%(float(dict_coordonnees[cle][0])))
#                             fichier_py.write("e1.coordonnee_y_groupe = %f\n"%(float(dict_coordonnees[cle][1])))
#                             fichier_py.write("e1.save()\n")
#     with open("script_ajout_attribut_noeud.py", 'w') as fichier_py :
#         with open("liste_extensions.pickle", 'rb') as fichier_pickle :
#             mon_depickler = pickle.Unpickler(fichier_pickle)
#             liste_ext = mon_depickler.load()
#             compteur = 0
#             fichier_py.write('from aminor.models import Extension\n\n')
#             for elt in liste_ext :
#                 print("petit rat")
#                 if elt in GROUPE_GNRA :
#                     fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+elt+ "')\n")
#                     fichier_py.write("e1.groupe = 'GNRA' \n")
#                     fichier_py.write("e1.save() \n")
#                     compteur += 1
#                     print("a")
#                 if elt in GROUPE_ARICH : 
#                     fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+elt+ "')\n")
#                     fichier_py.write("e1.groupe = 'ARICH' \n")
#                     fichier_py.write("e1.save() \n")
#                     compteur += 1
#                     print('b')
#                 if elt in GROUPE_GBULGE :
#                     fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+elt+ "')\n")
#                     fichier_py.write("e1.groupe = 'GBULGE' \n") 
#                     fichier_py.write("e1.save() \n") 
#                     compteur += 1 
#                     print("c")    
#                 if elt in GROUPE_GNRA_ETENDU :
#                     fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+elt+ "')\n")
#                     fichier_py.write("e1.groupe = 'GNRA-p' \n")
#                     fichier_py.write("e1.save() \n")   
#                     compteur += 1
#                     print('d')
#                 if elt in GROUPE_ARICH_ETENDU :
#                     fichier_py.write("e1 = Extension.objects.get(nom_molecule='"+elt+ "')\n")
#                     fichier_py.write("e1.groupe = 'ARICH-p' \n")  
#                     fichier_py.write("e1.save() \n") 
#                     compteur += 1
#                     print("e")
#                     
#             for elt in GROUPE_ARICH_ETENDU :
#                 if elt not in liste_ext :
#                     print(elt)
#             #liste_ext.remove('4V9F_0_48_30')
#     
# #     with open("liste_extensions.pickle", 'wb') as fichier_pickle :
# #         mon_pickler = pickle.Pickler(fichier_pickle)
# #         mon_pickler.dump(liste_ext)       
#     print(len(liste_ext))    
#     print(compteur)
#     
#     print(len(GROUPE_ARICH))
#     print(len(GROUPE_GNRA))
#     print(len(GROUPE_GBULGE))
#     print(len(GROUPE_ARICH_ETENDU))
#     print(len(GROUPE_GNRA_ETENDU))

    
#     with open("fichier_python_data_coord_extension.py", 'w') as fichier_py :
#         fichier_py.write("from aminor.models import Extension\n")
#         fichier_py.write("from aminor.models import Jmol\n\n")
#         with open("graphs_2.92.pickle", 'rb') as fichier_graphes :
#             mon_depickler_graphes = pickle.Unpickler(fichier_graphes)
#             graphes = mon_depickler_graphes.load()
#             for i in range(1,11) :
#                 for fic in os.listdir(EXTENSION_PATH_TAILLE%i) :
#                     if "pickle" in fic and "couples_possibles" not in fic and len(fic.split("_")) == 5 :
#                         with open(EXTENSION_PATH_TAILLE%i+fic, 'rb') as fichier_pickle :
#                             mon_depickler = pickle.Unpickler(fichier_pickle)
#                             graphe = mon_depickler.load()
#                             cle = (fic.split("_")[1], fic.split("_")[2])
#                             data_coord = motif3D_html_seul(cle, graphe, graphes)
#                             
#                             name_ext = fic.split("_")[1] + "_" + fic.split("_")[2] + "_" + fic.split("_")[3] + "_" + fic.split("_")[4][:len(fic.split("_")[4])-7]
#                             
#                             fichier_py.write("e = Extension.objects.get(nom_molecule='%s')\n"%name_ext)
#                             
#                             fichier_py.write("j = Jmol.objects.get(name_extension='%s', taille_extension=%s)\n"%(name_ext,i))
#                             fichier_py.write("e.url_jmol.add(j)\n")
#                             fichier_py.write("e.save()\n")
                        
    
    
    ## ajout des graphes version3
    
#     for j in range(9, 11) :
#         with open(EXTENSION_PATH%j+"sim_extensions_toutes_aretes_coeff_all1_taille_%d.pickle"%j, 'rb') as fichier_graphe :
#             mon_depickler_graphe = pickle.Unpickler(fichier_graphe)
#             dico_sim = mon_depickler_graphe.load()  
#             
#             for i in np.arange(0.0, 1.0, 0.1) :
#                 ajout_base_graphe_version_3(dico_sim, i, j, 'sim')
                
    
    with open(PATH_MMCIF+'fichiers_rmsd_taille_4_que_carbone_1_normalise.pickle', 'rb') as fichier_rmsd :
        mon_depickler_rmsd = pickle.Unpickler(fichier_rmsd)
        dico_rmsd = mon_depickler_rmsd.load() 
        
        dico_rmsd_cles = {}
        for cle in dico_rmsd.keys() :
            cle_0 = cle[0][8:len(cle[0])-13]
            cle_1 = cle[1][8:len(cle[1])-13]
            dico_rmsd_cles.update({(cle_0, cle_1) : dico_rmsd[cle]})
            
            
            
        for i in np.arange(0.0, 1.0, 0.1) :
            ajout_base_graphe_version_3(dico_rmsd_cles, i, 4, 'rmsd')
        
        
        
    


                    
        
                   
            
        
        