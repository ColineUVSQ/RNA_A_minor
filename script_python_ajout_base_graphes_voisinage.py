from aminor.models import Extension
from aminor.models import Lien
from aminor.models import Graph_voisinage

import django.core.exceptions

g = Graph_voisinage(graph_nom='toutes_aretes_coeffn1_a1_c1_0.6_2_voisinage', graph_poids= 0.6,graph_num=2)
g.save() 
e = Extension.objects.get(nom_molecule='1FJG_A_271_1')
g.graph_cluster.add(e) 
e = Extension.objects.get(nom_molecule='4V9F_0_30_23')
g.graph_cluster.add(e) 
e = Extension.objects.get(nom_molecule='4V9F_0_48_16')
g.graph_cluster.add(e) 
e = Extension.objects.get(nom_molecule='1FJG_A_109_6')
g.graph_cluster.add(e) 
e = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
g.graph_voisinage.add(e) 
e = Extension.objects.get(nom_molecule='5DM6_X_48_10')
g.graph_voisinage.add(e) 
e = Extension.objects.get(nom_molecule='3JCS_1_48_18')
g.graph_voisinage.add(e) 
e = Extension.objects.get(nom_molecule='4V88_A5_48_3')
g.graph_voisinage.add(e) 
e = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
g.graph_voisinage.add(e) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='4V88_A5_48_3')

try :
	l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
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
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
try :
	l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
try :
	l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
try :
	l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_271_1')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try :
	l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
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
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V88_A5_48_3')
e2 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
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
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
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
e1 = Extension.objects.get(nom_molecule='5FDU_1A_137_6')
e2 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
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
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
e2 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
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
e1 = Extension.objects.get(nom_molecule='5DM6_X_48_10')
e2 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
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
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
e2 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
e2 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
e2 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='4V9F_0_30_23')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
e2 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
e2 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='5J7L_DA_48_1')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
e2 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='1FJG_A_109_6')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
e1 = Extension.objects.get(nom_molecule='3JCS_1_48_18')
e2 = Extension.objects.get(nom_molecule='4V9F_0_48_16')
try :
	 l = Lien.objects.get(num_ext_1=e1, num_ext_2=e2)
except django.core.exceptions.ObjectDoesNotExist :
	l = Lien.objects.get(num_ext_1=e2, num_ext_2=e1)
g.graphe.add(l) 
