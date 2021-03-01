from __future__ import division
import ast
import os
import csv
import gzip

##update files
file_name= "variant_file_path_without_surfix(.vcf)" #This file include variants with >=500/1000 quality for substitutions/indels, and variants in regions with >=20x coverge in >=90% of samples
                                                    #The file format is: chr\tpos\tref\talt\tgene\t-\t-\tN\t-\t-\tUK\tN\t-\t-\t-\t-\t-\t-\t genotype genotype genotype...
                                                    #The file surfix should be ".vcf"

hcdiff= "annotated_file"                            #In-house annotation pipeline output, data should include- gene_name,variant_component,synonymous(T/F),gnomAD-G_AF,EXAC_AC_HOM

LOFTEE_output= "LOFTEE_output"                      #The file format is: chr pos ref alt HC/LC/" "

ADAR_genes= "list_of_324_ADAR_genes"                #attached

AR_genes= "list_of_1605_AR_genes"                   #attached

manualvars= "manual_curated_vars"                   #variants that were curated manually with their final classification. attached.

CLINVAR= "clinvar_vcf"                              #clinvar_20200127.vcf.gz

HGMD= "HGMD_database"                               #hgmd_pro_2018.3_hg19.vcf

HGMD_transcripts= "HGMD_with_Ens_transcriits"       #attached

VKGL= "VKGL_databse"                                #VKGL- Dutch database 

INTERVAR= "intervar_output"                         #Intervar output

samples_num= "number_of_samples_in_the_cohort (int)"    


def make_annotated_list(input_file):

    with open (input_file, 'rt') as file,  open(ADAR_genes,'rt') as extra_genes, open(manualvars,'rt') as manual, open (HGMD_transcripts,'rt')as l, gzip.open (CLINVAR, 'rt') as clinvar, open (HGMD, 'rt') as hgmd, open (VKGL,'rt') as classification_file, open (hcdiff, 'rt') as gene_file, open (INTERVAR,'rt') as inter,open (AR_genes,'rt') as gene_names, open(LOFTEE_output,'rt') as LOFTEE, open (file_name+'_annotated_leftovers.vcf', 'w') as filter5, open (file_name+'_annotated_AR.vcf', 'w') as final, open (file_name+'_annotated_ADAR.vcf', 'w') as final1:
	vcf = csv.reader(file, delimiter='\t')
	clin = csv.reader(clinvar, delimiter='\t')
	hgmd_data = csv.reader(hgmd, delimiter='\t')
	classification = csv.reader(classification_file, delimiter=',')
	annot = csv.reader(gene_file, delimiter='\t')
        intervar = csv.reader(inter, delimiter='\t')
	loftee = csv.reader(LOFTEE, delimiter='\t')
	genes = csv.reader(gene_names, delimiter='\t')
	ADAR_Genes = csv.reader(extra_genes, delimiter='\t')
	hgmd_trans = csv.reader(l, delimiter='\t')
	man = csv.reader(manual, delimiter='\t')
	writerAR = csv.writer(final, delimiter='\t')
        writerADAR = csv.writer(final1, delimiter='\t')
	l_writer = csv.writer(filter5, delimiter='\t')

	next(annot)
	next(classification)
	next(intervar)

	d_class= {}
	d_annot= {}
	d_var= {}
	d_clin= {}
	d_inter= {}
	d_hgmd= []
	d_loftee={}
	d={}
	d_man={}
	d_ens={}
	d_names={"HJV":"HFE2","PJVK":"DFNB59"}
        geneslistAR=[]
        geneslistADAR=[]

        for g in genes:                                                             #create 2 lists of all analysed genes 
	    geneslistAR.append(g[0])
	for g in ADAR_Genes:
            geneslistADAR.append(g[0])


        for w in hgmd_trans:                                                        #create a dictionary for HGMD canonical transcripts
            if w[0] not in d_ens.keys():
                if w[0] in d_names.keys():
                    d_ens[d_names[w[0]]]=ast.literal_eval(w[2])
                else:
                    d_ens[w[0]]=ast.literal_eval(w[2])
            else:
                if w[0] in d_names.keys():
                    d_ens[d_names[w[0]]].extend(ast.literal_eval(w[2]))
                else:
                    d_ens[w[0]].extend(ast.literal_eval(w[2]))


        for p in man:                                                               #create a dictionary for manually curated variants {var:calssification}
            if len(p[2])==len(p[3]):
		kvcf= p[0] + "_" + p[1] + "_" + p[2] + "_" + p[3] 
	    elif len(p[2])>len(p[3]):
		kvcf= p[0] + "_" + str(int(p[1])+1) + "_" + p[2][1:] + "_."
	    elif len(p[2])<len(p[3]):
		kvcf= p[0] + "_" + p[1] + "_._" + p[3][1:]
            d_man[kvcf]=p[4]

            	
	for r in classification:                                                    #create a dictionary for VKGL variants {var:calssification}. CHANGE INDEX IF NEEDED.
	    if len(r[4])>1 and len(r[5])>1: 
		continue
	    elif len(r[4])==1 and len(r[5])==1:
                var=  "chr" + r[1] +  "_" + r[2] +"_" + r[4] + "_" + r[5]
            elif len(r[4])> len(r[5])==1:    
		var=  "chr" + r[1] + "_" + str(int(r[2])+1) + "_" + r[4][1:] + "_."
	    elif len(r[5])> len(r[4]):
		var=  "chr" +  r[1] + "_" + r[2] +  "_._" + r[5][1:]                        
		
	    if var!="":
                d_class[var]= r[10]


	for c in clin:                                                              #create a dictionary for ClinVar data {var:[classification, status_review]}. CHANGE INDEX IF NEEDED.
	    if len(c[0])==1 or len(c[0])==2:
	    	if len(c[3])>1 and len(c[4])>1:
                    continue	
                elif len(c[3])==len(c[4]):
                    var= "chr" +  c[0] + "_" + c[1] + "_" + c[3] + "_" + c[4]
                elif len(c[3])>len(c[4]):
                    var= "chr" +  c[0] + "_" + str(int(c[1])+1) + "_" + c[3][1:] + "_."
                elif len(c[3])<len(c[4]):
                    var= "chr" +  c[0] + "_" + c[1] + "_._" + c[4][1:]
	
                inf= c[7].split(';')
		if 'CLNSIGINCL' in inf:
		    d_clin[var]=["-","-"]
		    break
	        else:
		    clin_clas="-"
		    clin_rev="-"
		    for i in inf:
			if "CLNSIG=" in i:
			    clin_clas= i.split('=')[1]
			elif "CLNREVSTAT=" in i:
			    clin_rev= i.split('=')[1]
		    d_clin[var]=[clin_clas, clin_rev]

	
	for h in hgmd_data:                                                         #create a list of variants with a "DM" flag in HGMD. CHANGE INDEX IF NEEDED.
            if len(h[0])==1 or len(h[0])==2:
                if len(h[3])>1 and len(h[4])>1:
		    continue 
 		elif len(h[3])==1 and len(h[4])==1:
                    var= "chr" + h[0] + "_" + h[1] + "_" + h[3] + "_" + h[4]
                elif len(h[3])>len(h[4]) :
                    var= "chr" +  h[0] + "_" +str(int(h[1])+1)  +  "_" + h[3][1:] + "_."
		elif len(h[3])<len(h[4]):
		    var= "chr" +  h[0] + "_" + h[1]  +  "_._" + h[4][1:] 

		for inf in h[7].split(";"):
		    if "CLASS" in inf and inf.split("=")[1]=="DM":
		        d_hgmd.append(var)
			break


	for t in annot:                                                             #create a dictionary with annotated data from our In-house annotation tool {var:[gene_name,variant_component,synonymous(T/F),gnomAD-G_AF,EXAC_AC_HOM]}. CHANGE INDEX IF NEEDED.
            if t!=[] and "chr" in t[0]:
                lof=0
                missense=0
                other=0
                varList=[]
                rightV=[]
                
                LOF=["stop_gained" ,"start_lost" , "frameshift_variant" ,"splice_donor_variant" , "splice_acceptor_variant"]
                Missense= ["missense_variant"]
                
                Gene=t[14]
                tmpL=t[25].split(" ")
	    			
		if len(t[3])==1 and len(t[4])==1:
		    var= t[0] + "_" + t[1] + "_" + t[3] + "_" + t[4]
		elif len(t[3])>len(t[4]):
		    var= t[0] + "_" + t[1] + "_" + t[3] + "_."
		elif len(t[3])<len(t[4]):
		    var= t[0] + "_" + t[1] + "_._" + t[4]

		d_annot[var]= [t[14], t[20], t[28], t[140], t[81], "-"]

                if var in d_man.keys():
                        d_annot[var][-1]= d_man[var]
                        continue

                for u in tmpL:
                    if len(u)>1:
                        ens=u.split("(")[0].split(".")[0]
                        vType=u.split("(")[1][:-1]
                        varList.append([ens,vType])
                    
                if Gene in d_ens.keys():   
                    for pair in varList:
                        if pair[0] in d_ens[Gene]:
                            rightV.append(pair[1])
                            
                    if len(set(rightV))==1 and rightV!=[]:        
                        d_annot[var][-1]= rightV[0]
                        continue
                    elif rightV!=[]:
                        for j in rightV:
                            if j in LOF:
                                lof+=1
                            elif j in Missense:
                                missense+=1
                            else:
                                other+=1
                        if (lof>0 and missense==0 and other==0) or (lof==0 and missense>0 and other==0) or (lof==0 and missense==0 and other>0):
                            d_annot[var][-1]= t[24].split(";")[0]
                            continue
                        else:
                            perlof= lof/(lof+missense+other)*100
                            permis= missense/(lof+missense+other)*100
                            perother= other/(lof+missense+other)*100
                            if perlof>50:
                                d_annot[var][-1]=t[24].split(";")[0]
                                continue
                            elif permis>50:
                                d_annot[var][-1]= "missense_variant"
                                continue
                            elif perother>50:
                                d_annot[var][-1]= "other"
                                continue
                            else:
                                d_annot[var][-1]=rightV
                            continue
                    else:
                        for pair in varList:
                            if pair[1] in LOF:
                                lof+=1
                            elif pair[1] in Missense:
                                missense+=1
                            else:
                                other+=1

                        if (lof>0 and missense==0 and other==0) or (lof==0 and missense>0 and other==0) or (lof==0 and missense==0 and other>0):
                            d_annot[var][-1]=t[24].split(";")[0]
                            continue
                        else:
                            perlof= lof/(lof+missense+other)*100
                            permis= missense/(lof+missense+other)*100
                            perother= other/(lof+missense+other)*100
                            if perlof>50:
                                d_annot[var][-1]=t[24].split(";")[0]
                                continue
                            elif permis>50:
                                d_annot[var][-1]= "missense_variant"
                                continue
                            elif perother>50:
                                d_annot[var][-1]="other"
                                continue
                            else:
                                for pair in varList:
                                    rightV.append(pair[1])
                                d_annot[var][-1]= rightV
                                continue

                else:
                    
                    for pair in varList:
                        if pair[1] in LOF:
                            lof+=1
                        elif pair[1] in Missense:
                            missense+=1
                        else:
                            other+=1

                    if (lof>0 and missense==0 and other==0) or (lof==0 and missense>0 and other==0) or (lof==0 and missense==0 and other>0):
                        d_annot[var][-1]=t[24].split(";")[0]
                        continue
                    else:
                        perlof= lof/(lof+missense+other)*100
                        permis= missense/(lof+missense+other)*100
                        perother= other/(lof+missense+other)*100
                        if perlof>50:
                            d_annot[var][-1]=t[24].split(";")[0]
                            continue
                        elif permis>50:
                            d_annot[var][-1]= "missense_variant"
                            continue
                        elif perother>50:
                            d_annot[var][-1]="other"
                            continue
                        else:
                            for pair in varList:
                                rightV.append(pair[1])
                            d_annot[var][-1]= rightV
                            continue

                                                                                                    
	for l in loftee:                                                            #create a dictionary for LOFTEE info {var:Loftee_HC/LC}. CHANGE INDEX IF NEEDED.
            if len(l[2])==1 and len(l[3])==1:
                var= "chr" +l[0] + "_" + l[1] + "_" + l[2] + "_" + l[3]
            elif len(l[2])>len(l[3]):
                var= "chr" +l[0] + "_" + l[1] + "_" + l[2] + "_."
            elif len(l[2])<len(l[3]):
                var= "chr" +l[0] + "_" + l[1] + "_._" + l[3]
            
	    d_loftee[var]=l[4]


	for d in intervar:                                                          #create a dictionary for intervar classification {var:intervar_classification}. CHANGE INDEX IF NEEDED. 
            if len(d[3])>1 and len(d[4])>1:
		continue
	    elif len(d[3])==len(d[4]) and d[3]!="-" and d[4]!="-":
                var=  "chr" + d[0] + "_" + d[1] + "_" + d[3] + "_" + d[4]
	    elif (len(d[3])==1 and len(d[4])==1 and d[4]=="-") or len(d[3])>len(d[4]):
                var=  "chr" +  d[0] + "_" + d[1] +"_" + d[3] + "_."
	    elif (len(d[3])==1 and len(d[4])==1 and d[3]=="-") or len(d[3])<len(d[4]):
                var=  "chr" +  d[0] + "_" + d[1] +"_._" + d[4]
			
	    d_inter[var]= (d[13].split(": ")[1]).split(" PVS1")[0]


	for row in vcf:                                                             #scan the input file one row at a time and add all information
	    nl= row[:]
	    c_het= 0
	    c_hom= 0

	    if (len(row[2])>1 and len(row[3])>1):
		continue
	    elif len(row[2])==len(row[3]):
		kvcf= row[0] + "_" + row[1] + "_" + row[2] + "_" + row[3] 
	    elif len(row[2])>len(row[3]):
		kvcf= row[0] + "_" + str(int(row[1])+1) + "_" + row[2][1:] + "_."
	    elif len(row[2])<len(row[3]):
		kvcf= row[0] + "_" + row[1] + "_._" + row[3][1:]
                        
	    if kvcf in d_annot.keys():
		nl[4]= d_annot[kvcf][0]
		nl[5]= d_annot[kvcf][1]
		nl[6]= d_annot[kvcf][2]
		nl[8]= d_annot[kvcf][3]
		nl[9]= d_annot[kvcf][4]
		nl[14]= d_annot[kvcf][5]
	    else:
		nl[4]= '-'
		nl[5]= '-'
		nl[6]= '-'
		nl[8]= '-'	
		nl[9]= '-'
		nl[14]= '-'

	    if kvcf in d_class.keys():
                nl[10]= d_class[kvcf]

            if kvcf in d_hgmd:
                nl[11]= "Y"

	    if kvcf in d_clin.keys():
		nl[12]= d_clin[kvcf][0]
		if "no_assertion" in d_clin[kvcf][1]:
		    nl[13]= 0
		elif "single" in d_clin[kvcf][1] or "conflicting" in d_clin[kvcf][1]:
		    nl[13]= 1
		elif "multiple" in d_clin[kvcf][1]:
		    nl[13]= 2
		elif "expert" in d_clin[kvcf][1]:
                    nl[13]= 3
		elif "practice" in d_clin[kvcf][1]:
                    nl[13]= 4
		else:
		    nl[13]= "-"
	    else:
		nl[12]= "-"
		nl[13]= "-"
         
	    if kvcf in d_inter.keys():
		if "significance" in d_inter[kvcf]:
		    nl[15]= "Uncertain_significance"
		elif "benign" in d_inter[kvcf]:
                    nl[15]= "Likely_benign"
		elif "Benign" in d_inter[kvcf]:
                    nl[15]= "Benign"
		elif "pathogenic" in d_inter[kvcf]:
                    nl[15]= "Likely_pathogenic"
		elif "Pathogenic" in d_inter[kvcf]:
                    nl[15]= "Pathogenic"
	
	    if kvcf in d_loftee.keys():
		nl[16]= d_loftee[kvcf]
	    else:
		nl[16]= "-"
	    
            for i in range(19,len(row)):
                if row[i]=="0/1" or row[i]=="1/0":
                    c_het+=1

                elif row[i]=="1/1":
                    c_hom+=1
            
	    nl[17]= c_het
            nl[18]= c_hom
            
	    #write only non-frequent variants in the output file
	    if nl[4] in geneslistAR:		
                if c_het<float(samples_num)*0.05 and c_hom<float(samples_num)*0.01 or (kvcf=="chr3_15686693_G_C" or kvcf=="chr6_26093141_G_A" or kvcf=="chr14_94847262_T_A") :
		    writerAR.writerow(nl)
	        elif not (kvcf=="chr3_15686693_G_C" or kvcf=="chr6_26093141_G_A" or kvcf=="chr14_94847262_T_A"):
		    l_writer.writerow(nl)
	    if nl[4] in geneslistADAR:		
                if c_het<float(samples_num)*0.05 and c_hom<float(samples_num)*0.01 or (kvcf=="chr3_15686693_G_C" or kvcf=="chr6_26093141_G_A" or kvcf=="chr14_94847262_T_A") :
		    writerADAR.writerow(nl)
	        elif not (kvcf=="chr3_15686693_G_C" or kvcf=="chr6_26093141_G_A" or kvcf=="chr14_94847262_T_A"):
		    l_writer.writerow(nl)
	    

def make_plp_list1(input_file):
    #make separate files for substitutions and indels AR/ADAR genes
    os.system("cat " +input_file+ '''_annotated_AR.vcf|awk -v FS="\t" '{if (length($4)==1 && length($3)==1) print}'>'''+input_file+ "_annotated_subs_AR.vcf")
    os.system("cat " +input_file+ '''_annotated_AR.vcf|awk -v FS="\t" '{if (!(length($4)==1 && length($3)==1)) print}'>'''+input_file+ "_annotated_indels_AR.vcf")
    os.system("cat " +input_file+ '''_annotated_ADAR.vcf|awk -v FS="\t" '{if (length($4)==1 && length($3)==1) print}'>'''+input_file+ "_annotated_subs_ADAR.vcf")
    os.system("cat " +input_file+ '''_annotated_ADAR.vcf|awk -v FS="\t" '{if (!(length($4)==1 && length($3)==1)) print}'>'''+input_file+ "_annotated_indels_ADAR.vcf")

    #make plp list for the AR/ADAR subs
    os.system ("cat " +input_file+ '''_annotated_subs_AR.vcf |awk -v FS="\t" '{if (($11=="LP" || ($13~"athogenic" && $14>=2)) || (($6~"CANONICAL"||$6=="EXON_REGION") && ($7=="FALSE") && ($9<=0.01) && ($15~"stop"||$15~"start"||$15~"frameshift"||$15~"splice_donor"||$15~"splice_acceptor")) || (($11!="LB" && !($13~"enign" && $14>=2))&&(($16~"athogenic" && ($13~"Pathogenic"||$13=="Likely_pathogenic")) || ($12=="Y" && ($13~"Pathogenic"||$13=="Likely_pathogenic")) || ($12=="Y" && $16~"athogenic"))))print}'>''' +input_file+ "_annotated_subs_presumable_plp.vcf")
    os.system ("cat " +input_file+ '''_annotated_subs_ADAR.vcf |awk -v FS="\t" '{if ((($11=="LP" || ($13~"athogenic" && $14>=2))&& ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor")) || (($6~"CANONICAL"||$6=="EXON_REGION") && ($7=="FALSE") && ($9<=0.01) && ($15~"stop"||$15~"start"||$15~"frameshift"||$15~"splice_donor"||$15~"splice_acceptor")) || (($11!="LB" && !($13~"enign" && $14>=2))&&(($16~"athogenic" && ($13~"Pathogenic"||$13=="Likely_pathogenic")) || ($12=="Y" && ($13~"Pathogenic"||$13=="Likely_pathogenic")) || ($12=="Y" && $16~"athogenic"))&& ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor")))print}'>>''' +input_file+ "_annotated_subs_presumable_plp.vcf")  

    #make plp list for the AR indels - part 1-2  
    os.system("cat " +input_file+ '''_annotated_indels_AR.vcf |awk -v FS="\t" '{if (($13~"athogenic" && $14>=2)||$11=="LP") print }'> '''+input_file+ '''_annotated_indels_path1_AR.vcf''')
    os.system("cat " +input_file+ '''_annotated_indels_AR.vcf |awk -v FS="\t" '{if (!($11=="LP"||($13~"athogenic" && $14>=2)) && (($6~"CANONICAL" || $6=="EXON_REGION") && ($9<=0.01) && ($15~"stop"||$15~"start"||$15~"frameshift"||$15~"splice_donor"||$15~"splice_acceptor")) && !(length($3)>10 || length($4)>10)) print}'>''' +input_file+ '''_annotated_indels_partial_path2_AR.vcf''')

    #make plp list for the ADAR indels - part 1-2
    os.system("cat " +input_file+ '''_annotated_indels_ADAR.vcf |awk -v FS="\t" '{if ((($13~"athogenic" && $14>=2)||$11=="LP")&& ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor")) print }'> '''+input_file+ '''_annotated_indels_path1_ADAR.vcf''')
    os.system("cat " +input_file+ '''_annotated_indels_ADAR.vcf |awk -v FS="\t" '{if (!($11=="LP"||($13~"athogenic" && $14>=2)&& ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor")) && (($6~"CANONICAL" || $6=="EXON_REGION") && ($9<=0.01) && ($15~"stop"||$15~"start"||$15~"frameshift"||$15~"splice_donor"||$15~"splice_acceptor")) && !(length($3)>10 || length($4)>10)) print}'>''' +input_file+ '''_annotated_indels_partial_path2_ADAR.vcf''')

    
def make_plp_list2(input_file):
    ##edit indels_partial_path2_AR/ADAR file and remove manualy adajcent indels <10bp range and save as "file_name + '''annotated_indels_partial_path2_AR/ADAR_.txt
    #make plp list for the AR indels - part 2-3 
    os.system("dos2unix " +input_file+ "_annotated_indels_partial_path2_AR_.txt")    
    os.system("cat " +input_file+ '''_annotated_indels_partial_path2_AR_.txt|awk -v FS="\t" '{if($17!="LC")print}'>'''+input_file+ '''_annotated_indels_path2_AR.vcf''')
    os.system("cat " +input_file+ '''_annotated_indels_AR.vcf |awk -v FS="\t" '{if (!($11=="LP" || ($13~"athogenic" && $14>=2)) && !(($6~"CANONICAL"||$6=="EXON_REGION") && ($9<=0.01) && ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor"))&& !(length($3)>10 || length($4)>10) && (($11!="LB" && !($13~"enign" && $14>=2) && $17!="LC")&&(($16~"athogenic" && ($13~"Pathogenic" || $13=="Likely_pathogenic"))||($12=="Y" && ($13~"Pathogenic" || $13=="Likely_pathogenic"))||($12=="Y" && $16~"athogenic"))))print}'> '''+input_file+ '''_annotated_indels_path3_AR.vcf''')

    #combine all indels AR plps
    os.system("cat " +input_file+ '''_annotated_indels_path1_AR.vcf '''+input_file+ '''_annotated_indels_path2_AR.vcf '''+input_file+ '''_annotated_indels_path3_AR.vcf >''' +input_file+ '''_annotated_indels_presumable_plp.vcf''')

    #make plp list for the ADAR indels - part 2-3
    os.system("dos2unix " +input_file+ "_annotated_indels_partial_path2_ADAR_.txt")
    os.system("cat "+input_file+ '''_annotated_indels_partial_path2_ADAR_.txt|awk -v FS="\t" '{if($17!="LC")print}'>'''+input_file+ '''_annotated_indels_path2_ADAR.vcf''')
    os.system("cat " +input_file+ '''_annotated_indels_ADAR.vcf |awk -v FS="\t" '{if (!($11=="LP" || ($13~"athogenic" && $14>=2)&& ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor")) && !(($6~"CANONICAL"||$6=="EXON_REGION") && ($9<=0.01) && ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor"))&& !(length($3)>10 || length($4)>10) && (($11!="LB" && !($13~"enign" && $14>=2) && $17!="LC")&&(($16~"athogenic" && ($13~"Pathogenic" || $13=="Likely_pathogenic"))||($12=="Y" && ($13~"Pathogenic" || $13=="Likely_pathogenic"))||($12=="Y" && $16~"athogenic"))&& ($15~"stop" || $15~"start"|| $15~"frameshift"|| $15~"splice_donor"||$15~"splice_acceptor")))print}'> '''+input_file+ '''_annotated_indels_path3_ADAR.vcf''')

    #combine all indels ADAR plps
    os.system("cat " +input_file+ '''_annotated_indels_path1_ADAR.vcf '''+input_file+ '''_annotated_indels_path2_ADAR.vcf '''+input_file+ '''_annotated_indels_path3_ADAR.vcf >>''' +input_file+ '''_annotated_indels_presumable_plp.vcf''')

    #make final total plp file
    os.system("cat " +input_file+ '''_annotated_indels_presumable_plp.vcf ''' +input_file+ "_annotated_subs_presumable_plp.vcf > "  +input_file+ "_annotated_total_presumable_plp.vcf")


make_annotated_list(file_name+".vcf")
##make_plp_list1(file_name)
##make_plp_list2(file_name)

