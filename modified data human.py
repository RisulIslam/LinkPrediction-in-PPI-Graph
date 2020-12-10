import sys, codecs
from unidecode import unidecode
import re
import matplotlib.pyplot as plt

import networkx as nx
from statistics import median
from collections import Counter
import numpy as np

prefix_human_protein_name="9606.ENSP"
prefix_mouse_protein_name="10090.ENSMUSP"
suffix_number_of_digit=11

f=codecs.open("human.txt","r",encoding='utf-8')
lines=f.readlines()

fw=open("modified human lower500.txt","w")

def processline(lines):
    j=1
    for line in lines:
        line=line.replace('\n','')
        line=line.replace('\r','')
        

        t=line.split(' ')
        
        if(len(t)==3 and len(t[0])==(len(prefix_human_protein_name)+11) and len(t[1])==(len(prefix_human_protein_name)+11) and t[2].isdigit()==True and int(t[2])<=550):
            n1=t[0][-11:]
            n2=t[1][-11:]
            #print("In")
            fw.write(n1+" "+n2+"\n")
    
                

        j+=1
        if j%10000==0:
            print(j)
    fw.close()

processline(lines)
