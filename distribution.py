'''
Created on 28 nov. 2018

@author: coline
'''

import pickle
import seaborn as sns
import matplotlib.pyplot as plt
import csv

with open("fichier_max_test_3.pickle", 'rb') as fichier :
    mon_depickler = pickle.Unpickler(fichier)
    sim_max = mon_depickler.load()
    print(len(sim_max))
    print(sim_max)
    print(min(sim_max))
    print(max(sim_max))
    
    
    with open("fichier_sim_csv_6.csv", 'w', newline='') as fichier_csv :
        csvwriter = csv.writer(fichier_csv, delimiter=',')
        csvwriter.writerow(sim_max)
    
    sns.distplot(sim_max, bins=20)
    plt.show()

    #plt.hist(sim_max,bins=100,color="red",alpha=0.8)
    #plt.show()