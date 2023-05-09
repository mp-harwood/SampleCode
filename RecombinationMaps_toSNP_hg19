#Overlap recombination maps onto list of SNPs
#Separate recombination maps for African ancestry (YRI) and European ancestry (CEU)

import pandas as pd 

#####################
# 	Load African Recombination Map	#
#####################
data = open('YRI_recombination_allchr_CS_HRR_4Dec2019_final_sorted', 'r')
snps = data.readlines()
Afrec = []
for i in snps:
        within=[]
        new=i.split("\t")
        new[-1]= new[-1][:-1]
        Afrec.append(new)
Afrec=Afrec[1:]

#####################
# 	Load European Recombination Map	#
#####################
data = open('CEU_recombination_allchr_CS_HRR_4Dec2019_final_sorted', 'r')
snps = data.readlines()
Eurec = []
for i in snps:
        within=[]
        new=i.split("\t")
        new[-1]= new[-1][:-1]
        Eurec.append(new)
Eurec=Eurec[1:]

#####################
# 	Load SNP list of interest	#
#####################
data = open('output_liftover_hg38_to_hg19_ancestry.bed', 'r')
snps = data.readlines()
ASElist = []
for i in snps:
        within=[]
        new=i.split("\t")
        new[-1]= new[-1][:-1]
        ASElist.append(new)

#####################
# 	Map recombination class to list of SNPs	#
#####################

for i in range(1, len(ASElist)):
        if ASElist[i][4] =="3":
                k=0
                foundEU="no"
                while foundEU=="no" and k<=len(Eurec)-1:
                        if Eurec[k][0]==ASElist[i][0]: #check chr   
                                if float(Eurec[k][1]) <= float(ASElist[i][1]) <= float(Eurec[k][2]): #if within- must be CS/HRR
                                        ASElist[i].append(Eurec[k][3])
                                        foundEU="yes"
                                elif k < len(Eurec)-2:
                                        if Eurec[k][0] == Eurec[k+1][0]: #check chr
                                                if float(Eurec[k][2]) <= float(ASElist[i][1]) <= float(Eurec[k+1][1]): #if not within- then normal 
                                                        ASElist[i].append("Normal")
                                                        foundEU="yes"
                        k=k+1
        elif ASElist[i][4]=="2":
                l=0
                foundAf="no"
                while foundAf=="no" and l<= len(Afrec)-2:
                        if Afrec[l][0]==ASElist[i][0]: #check chr
                                if float(Afrec[l][1]) <= float(ASElist[i][1]) <= float(Afrec[l][2]): #if within- must be CS/HRR
                                        ASElist[i].append(Afrec[l][3])
                                        foundAf="yes"
                        elif l < len(Afrec)-1:
                                if Afrec[l][0] == Afrec[l+1][0]: #check chr
                                        if float(Afrec[l][2]) <= float(ASElist[i][1]) <= float(Afrec[l+1][1]): #if not within- then normal
                                                ASElist[i].append("Normal")
                                                foundAf="yes"
                        l=l+1
        

newfile_df= pd.DataFrame(ASElist)
newfile_df.to_csv("GTEX_WB_hg19_recombination", na_rep='NaN',  sep='\t', encoding='utf-8', index=False)



