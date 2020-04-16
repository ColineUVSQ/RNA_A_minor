from aminor.models import Extension
from aminor.models import Lien
from aminor.models import Graph_version_2
import django.core.exceptions

g = Graph_version_2(graph_nom='toutes_aretes_coeffn1_a1_c1_taille_max_non_can_2_2_0.7', graph_poids=0.7)
g.save() 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_9')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_197_4')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_9')
e2 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_9')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_197_3')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_9')
e2 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_9')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_21')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J5B_BA_48_14')
e2 = Extension.objects.get(nom_molecule='1FJG_A_48_17')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_62_12')
e2 = Extension.objects.get(nom_molecule='3JCS_1_25_16')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_62_12')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_62_14')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_62_12')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_62_5')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_62_12')
e2 = Extension.objects.get(nom_molecule='5DM6_X_328_2')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='3JCS_1_25_16')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_62_14')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='3JCS_1_25_16')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_62_5')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='3JCS_1_25_16')
e2 = Extension.objects.get(nom_molecule='5DM6_X_328_2')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='4V9F_0_137_5')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_272_1')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_272_2')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_272_1')
e2 = Extension.objects.get(nom_molecule='4V9F_0_207_3')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_30_15')
e2 = Extension.objects.get(nom_molecule='4V9F_0_30_4')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_30_15')
e2 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_30_15')
e2 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_30_15')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_30_17')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
e2 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_62_14')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
e2 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
e2 = Extension.objects.get(nom_molecule='4V9F_0_137_5')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_25_10')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_25_78')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_25_10')
e2 = Extension.objects.get(nom_molecule='4V88_A5_25_47')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_25_10')
e2 = Extension.objects.get(nom_molecule='5DM6_X_25_34')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_25_10')
e2 = Extension.objects.get(nom_molecule='3JCS_1_25_46')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_134_3')
e2 = Extension.objects.get(nom_molecule='5DM6_X_134_2')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_134_3')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_134_1')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_134_3')
e2 = Extension.objects.get(nom_molecule='4V9F_0_127_6')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_30_4')
e2 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_30_4')
e2 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_30_4')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_30_17')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A6_48_12')
e2 = Extension.objects.get(nom_molecule='5J5B_BA_48_7')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_197_4')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_197_3')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_197_4')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_21')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A6_17_55')
e2 = Extension.objects.get(nom_molecule='1FJG_A_58_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_48_13')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_48_25')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_48_13')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_20')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_48_13')
e2 = Extension.objects.get(nom_molecule='5DM6_X_48_28')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='3UCZ_R_62_15')
e2 = Extension.objects.get(nom_molecule='1U9S_A_58_11')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
e2 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_62_14')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_62_5')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_62_14')
e2 = Extension.objects.get(nom_molecule='5DM6_X_328_2')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_134_2')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_134_1')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_134_2')
e2 = Extension.objects.get(nom_molecule='4V9F_0_127_6')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_227_2')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_30')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_25_78')
e2 = Extension.objects.get(nom_molecule='4V88_A5_25_47')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_25_78')
e2 = Extension.objects.get(nom_molecule='5DM6_X_25_34')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_25_78')
e2 = Extension.objects.get(nom_molecule='3JCS_1_25_46')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_25_68')
e2 = Extension.objects.get(nom_molecule='4V88_A5_25_30')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_25_68')
e2 = Extension.objects.get(nom_molecule='4V9F_0_25_56')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_25_68')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_25_12')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_48_25')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_20')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_48_25')
e2 = Extension.objects.get(nom_molecule='5DM6_X_48_28')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_15')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
e2 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_48_19')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
e2 = Extension.objects.get(nom_molecule='1U9S_A_58_11')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_48_8')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_30_17')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J5B_BA_138_2')
e2 = Extension.objects.get(nom_molecule='1FJG_A_138_3')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='4V9F_0_137_5')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_197_3')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_21')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_44_3')
e2 = Extension.objects.get(nom_molecule='1FJG_A_58_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_25_30')
e2 = Extension.objects.get(nom_molecule='4V9F_0_25_56')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_25_30')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_25_12')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_134_1')
e2 = Extension.objects.get(nom_molecule='4V9F_0_127_6')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_48_15')
e2 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_48_15')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_48_19')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_48_15')
e2 = Extension.objects.get(nom_molecule='1U9S_A_58_11')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_127_7')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_74_7')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_48_19')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
e2 = Extension.objects.get(nom_molecule='1U9S_A_58_11')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J5B_BA_48_23')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_30_17')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_294_1')
e2 = Extension.objects.get(nom_molecule='5J5B_BA_294_2')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_62_5')
e2 = Extension.objects.get(nom_molecule='5DM6_X_328_2')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_48_19')
e2 = Extension.objects.get(nom_molecule='1U9S_A_58_11')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_48_20')
e2 = Extension.objects.get(nom_molecule='5DM6_X_48_28')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_25_56')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_25_12')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J5B_BA_58_3')
e2 = Extension.objects.get(nom_molecule='1FJG_A_58_23')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_137_5')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_134_5')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_74_7')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_25_47')
e2 = Extension.objects.get(nom_molecule='5DM6_X_25_34')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_25_47')
e2 = Extension.objects.get(nom_molecule='3JCS_1_25_46')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_25_34')
e2 = Extension.objects.get(nom_molecule='3JCS_1_25_46')
try : 
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
