import os
import sys
import subprocess
import shutil
import psutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
#######################################Importing packages ######################################
#######################################  Information  ##########################################
print('''
         \t\t#########################################
         \t\t#    1. Package name: PY-GEMINI         #
         \t\t#    2. Please Cite:                    #
         \t\t#    (Protein Files only)               #
         \t\t#    Thanku.                            #
         \t\t#########################################
      ''')
##################################################  Path Information ############################
print('**Please Input the directory and number of cores for the operation**\n\n\n')
input_path=input('Enter the path of the directory containing all the protein sequences\n\n\n')
################################################# core detection ################################
cr=psutil.cpu_count()
cores=int(input(('Your system has '+str(cr)+' CORES.\nPlease choose any value lower to it.\n')))
if(cores>1):
    nc=cores
else:
    nc=1
################################### Copying input file ##########################################
print('\n\n\nProcess A.\n**Creating a copy of your input file, as blast will generate subsidiary files in your sequence directory**\n\n\n')
try:
    in_path=('/'.join(input_path.split('/')[:-1]))
    shutil.copytree((input_path),os.path.join(in_path,'copy_of_'+(input_path.split('/')[-1])))
except FileExistsError:
    shutil.rmtree(os.path.join(in_path,'copy_of_'+(input_path.split('/')[-1])))
    in_path=('/'.join(input_path.split('/')[:-1]))
    shutil.copytree((input_path),os.path.join(in_path,'copy_of_'+(input_path.split('/')[-1])))
################################## Creating result directory ####################################
print('Process B.\n**Creating a folder named "Result" to save all your results**\n')
try:
    in_path=('/'.join(input_path.split('/')[:-1]))
    os.mkdir(os.path.join(in_path,'Result'))
    out_path=(os.path.join(in_path,'Result'))
except FileExistsError:
    in_path=('/'.join(input_path.split('/')[:-1]))
    out_path=(os.path.join(in_path,'Result'))
    shutil.rmtree(out_path)
    os.mkdir(out_path)
print('Result saved at:'+str(os.path.join(in_path,'Result'))+'\n\n\n')
                        ############   Re-naming  ###################
os.mkdir(os.path.join(out_path,'first_blast'))
items0=os.listdir(input_path)
for eachres in items0:
    resname=eachres.split('.')[0]
    fresname=resname+'.fasta'
    os.rename(os.path.join(input_path,eachres),os.path.join(input_path,fresname))
##################################### Blasting ##################################################
print('Process C.\n**Running BLAST**\n\n\n')
items=os.listdir(input_path)
os.mkdir(os.path.join(out_path,'Blast_log'))
for i in range(len(items)):
    os.mkdir(os.path.join(out_path,'first_blast',''.join(items[i].split('.')[0])))
    seqdb=(os.path.join(input_path,items[i]))
    cmd='makeblastdb -in '+seqdb+' -dbtype prot -title '+seqdb+'title -out '+seqdb+'db > '+os.path.join(out_path,'Blast_log','BLAST-logfile')
    subprocess.call(cmd, shell=True)
    for j in range(len(items)):
        if(items[j]!=items[i]):
            seqqu=(os.path.join(input_path,items[j]))
            cmd='blastp -query '+seqqu+' -db '+seqdb+'db -out '+(os.path.join(out_path,'first_blast',''.join(items[i].split('.')[0]),(''.join(items[i].split('.')[0])+'_'+''.join(items[j].split('.')[0]))))+'.tsv -num_threads '+str(nc)+' -outfmt 6'
            subprocess.call(cmd, shell=True)
##############################################  Blast data extraction  ##########################
print('Process D.\n**Blast data extraction**\n\n\n')
os.mkdir((os.path.join(out_path,'homogenome')))
entity=os.listdir(os.path.join(out_path,'first_blast'))
for i in range(len(entity)):
    line=[]
    member=os.listdir(os.path.join(out_path,'first_blast',entity[i]))
    for each in member:
        with open((os.path.join(out_path,'first_blast',entity[i],each)),'r') as f1:
            rf=f1.readlines()
            for every in rf:
                data=((''.join(every.split('\n'))).split('\t'))
                if(float(data[2])>=85 and int(data[5])<=3 and float(data[-2])<10e-5):
                    line.append(data[1])
    setline=set(line)
    for ele in setline:
        if(line.count(ele)==(len(items)-1)):
            with open((os.path.join(out_path,'homogenome',(entity[i]+'.fasta'))),'a') as f2:
                f2.write(ele+'\n')
################################### Counting the least number of core genes ######################
print('Process E.\n**Counting the Least Number of Genes (Homologs) to be considered for Core Gene interpretation**\n\n\n')
atoms=os.listdir(os.path.join(out_path,'homogenome'))
seqlen={}
for protons in atoms:
    with open(os.path.join(out_path,'homogenome',protons),'r') as file1:
        readf=file1.readlines()
    mem=len(readf)
    seqlen.update({protons:mem})
min_seq=(min(seqlen,key=seqlen.get))
################################### making the sequences #########################################
print('Process F.\n**Creating Fasta Sequences for homologous genes**\n\n\n')
os.mkdir(os.path.join(out_path,'homogenome_seq'))
sequences=[]
for electrons in atoms:
    with open(os.path.join(out_path,'homogenome',electrons),'r') as fl:#electrons are the filesnames
        headers=fl.readlines()#here we have opened the file and sequence names are taken in a list
    for ehed in headers:
        checkhead=(ehed.split('\n')[0])
        for eseq in SeqIO.parse(os.path.join(input_path,electrons),'fasta'):
            if(checkhead==eseq.id):
                sequences.append(SeqRecord(Seq(str(eseq.seq)),id=eseq.id,description=eseq.description))
    SeqIO.write(sequences,os.path.join(out_path,'homogenome_seq',''.join(electrons.split('.')[0]))+'.fasta','fasta')
    sequences=[]
################################  dummy second blast+Reciprocal blast ##############################
print('Process G.\n**Finding Core Genes in all the Fasta files and filtering the FASTA sequences for the same**\n\n\n')
start=[]
end=[]
vs_start=[]
vs_end=[]
os.mkdir(os.path.join(out_path,'Reciprocal'))
os.mkdir(os.path.join(out_path,'Reciprocal','reciprocal_blast_fullGenome'))
os.mkdir(os.path.join(out_path,'Reciprocal','reciprocal_blast_homologousGenome'))
os.mkdir(os.path.join(out_path,'Reciprocal','homologousPID=>85%Genome'))
os.mkdir(os.path.join(out_path,'Reciprocal','reciprocal_renamed'))
os.mkdir(os.path.join(out_path,'Reciprocal','reciprocal_genome'))
cseq=os.listdir(os.path.join(out_path,'homogenome'))
dseq=os.listdir(os.path.join(out_path,'first_blast'))
for echdseq in dseq:
    infile=os.listdir(os.path.join(out_path,'first_blast',echdseq))
    for nums in range(len(infile)):
        if(min_seq.split('.fasta')[0] in infile[nums]):
            shutil.copyfile(os.path.join(out_path,'first_blast',echdseq,infile[nums]),os.path.join(out_path,'Reciprocal','reciprocal_blast_fullGenome',infile[nums]))
rseq=os.listdir(os.path.join(out_path,'Reciprocal','reciprocal_blast_fullGenome'))
for echrseq in rseq:
    with open(os.path.join(out_path,'Reciprocal','reciprocal_blast_fullGenome',echrseq),'r') as dfile1:
        rd_dfile1=dfile1.readlines()
    for echcseq in cseq:
        if(echcseq.split('.fasta')[0] in echrseq):
            with open(os.path.join(out_path,'homogenome',echcseq),'r') as dfile2:
                rd_dfile2=dfile2.readlines()
    for nums in range(len(rd_dfile2)):
        numel=rd_dfile2[nums].split('\n')[0]
        for gums in range(len(rd_dfile1)):
            if(numel in rd_dfile1[gums]):
                with open((os.path.join(out_path,'Reciprocal','reciprocal_blast_homologousGenome',echrseq)),'a') as dfile3:
                    dfile3.write(rd_dfile1[gums])
rcpg=os.listdir(os.path.join(out_path,'Reciprocal','reciprocal_blast_homologousGenome'))
for echrcpg in rcpg:
    rcp_ar1=[]
    rcp_ar2=[]
    if(echrcpg.startswith(min_seq.split('.fasta')[0])==False):
        with open(os.path.join(out_path,'Reciprocal','reciprocal_blast_homologousGenome',echrcpg),'r') as read1:
            rcp_read1=read1.readlines()
        for rcp_each1 in rcp_read1:
            datum1=((''.join(rcp_each1.split('\n'))).split('\t'))
            if(float(datum1[2])>=85 and int(datum1[5])<=3 and float(datum1[-2])<10e-5):
                (rcp_ar1.append(str(datum1[0]+':'+datum1[1])))
        for echrcp_ar1 in rcp_ar1:
            with open(os.path.join(out_path,'Reciprocal','homologousPID=>85%Genome',echrcpg),'a') as dfile4:
                dfile4.write(echrcp_ar1+'\n')
    elif(echrcpg.startswith(min_seq.split('.fasta')[0])==True):
        with open(os.path.join(out_path,'Reciprocal','reciprocal_blast_homologousGenome',echrcpg),'r') as read2:
            rcp_read2=read2.readlines()
        for rcp_each2 in rcp_read2:
            datum2=((''.join(rcp_each2.split('\n'))).split('\t'))
            if(float(datum2[2])>=85 and int(datum2[5])<=3 and float(datum2[-2])<10e-5):
                (rcp_ar2.append(str(datum2[1]+':'+datum2[0])))
        for echrcp_ar2 in rcp_ar2:
            with open(os.path.join(out_path,'Reciprocal','homologousPID=>85%Genome',echrcpg),'a') as dfile5:
                dfile5.write(echrcp_ar2+'\n')
rchgnm=os.listdir(os.path.join(out_path,'Reciprocal','homologousPID=>85%Genome'))
for echrchgnm in rchgnm:
    if(echrchgnm.startswith(min_seq.split('.fasta')[0])):
        with open(os.path.join(out_path,'Reciprocal','homologousPID=>85%Genome',echrchgnm),'r') as ofile1:
            rd_ofile1=ofile1.readlines()
        for echrd_ofile1 in rd_ofile1:
            liness1=(''.join(echrd_ofile1.split('\n'))).split(':')
            with open(os.path.join(out_path,'Reciprocal','reciprocal_renamed','1.'+echrchgnm),'a') as ofile2:
                ofile2.write(str(liness1[1]+':'+liness1[0]+'\n'))
    if(echrchgnm.startswith(min_seq.split('.fasta')[0])==False):
        icespice=echrchgnm.split('_'+min_seq.split('.fasta')[0]+'.tsv')[0]
        nameice='2.'+min_seq.split('.fasta')[0]+'_'+icespice+'.tsv'
        shutil.copyfile(os.path.join(out_path,'Reciprocal','homologousPID=>85%Genome',echrchgnm),os.path.join(out_path,'Reciprocal','reciprocal_renamed',nameice))
sorgnom=os.listdir(os.path.join(out_path,'Reciprocal','reciprocal_renamed'))
for soreach in sorgnom:
    if(soreach.startswith('1.')):
        nameselect1=soreach.split('1.')[1]
        with open(os.path.join(out_path,'Reciprocal','reciprocal_renamed',('1.'+nameselect1)),'r') as srgfile1:
            rd_srgfile1=srgfile1.readlines()
        with open(os.path.join(out_path,'Reciprocal','reciprocal_renamed',('2.'+nameselect1)),'r') as srgfile2:
            rd_srgfile2=srgfile2.readlines()
    elif(soreach.startswith('2.')):
        continue
    for echrcpa in rd_srgfile1:
        cmpdat1=(echrcpa.split('\n')[0]).split(':')[0]
        for echrcpb in rd_srgfile2:
            cmpdat2=(echrcpb.split('\n')[0]).split(':')[1]
            if((cmpdat1)==(cmpdat2)):
                with open((os.path.join(out_path,'Reciprocal','reciprocal_genome',min_seq.split('.fasta')[0]+'_vs'+nameselect1.split(min_seq.split('.fasta')[0])[1].rstrip('.tsv')+'.txt')),'a') as fh1:
                    fh1.write(echrcpb)
                break
#################################### counting numbers #########################################
lisech=os.listdir(os.path.join(out_path,'Reciprocal','reciprocal_genome'))
rra=[]
corra=[]
for echlis in lisech:
    with open(os.path.join(out_path,'Reciprocal','reciprocal_genome',echlis)) as philoe1:
        rdphiloe=philoe1.readlines()
    for slay in rdphiloe:
        escdat=(slay).split(":")[0]
        rra.append(escdat)
srra=set(rra)
for elsrra in srra:
    if(rra.count(elsrra)==len(lisech)):
        corra.append(elsrra.rstrip("\n"))
############################### Set of Core Genes ###############################################
os.mkdir(os.path.join(out_path,'CORE_Genes'))
os.mkdir(os.path.join(out_path,'CORE_Genes','gene_names'))
os.mkdir(os.path.join(out_path,'CORE_Genes','gene_sequences'))
os.mkdir(os.path.join(out_path,'CORE_Genes','msa_seq'))
for echseqlis in lisech:
    echfilname=''.join(''.join(echseqlis.split("_vs_")[1]).split(".txt")[0])
    with open(os.path.join(out_path,'Reciprocal',"reciprocal_genome",echseqlis)) as philoe2:
        rdphiloe2=philoe2.readlines()
        for xe in range(len(corra)):
            for x in range(len(rdphiloe2)):
                echlinx=rdphiloe2[x].split('\n')[0]
                if(corra[xe] in echlinx):
                    spechlinx=echlinx.split(":")[1]
                    if(echseqlis==lisech[-1]):
                        with open(os.path.join(out_path,'CORE_Genes','gene_names',(min_seq.split('.fasta')[0])),'a') as philoe4:
                            philoe4.write(corra[xe]+'\n')
                    with open(os.path.join(out_path,'CORE_Genes','gene_names',echfilname),'a') as philoe3:
                        philoe3.write(spechlinx+'\n')
lseqin=os.listdir(os.path.join(out_path,'CORE_Genes','gene_names'))
for echlseqin in lseqin:
    chequences=[]
    nmechlseqin=echlseqin+'.fasta'
    with open(os.path.join(out_path,'CORE_Genes','gene_names',echlseqin),'r') as philoe4:
        rdphiloe4=philoe4.readlines()
    for eachnos in rdphiloe4:
        philoe4line=eachnos.rstrip('\n')
        for pseq in SeqIO.parse(os.path.join(out_path,'homogenome_seq',nmechlseqin),'fasta'):
            if(philoe4line==pseq.id):
                chequences.append(SeqRecord(Seq(str(pseq.seq)),id=pseq.id,description=pseq.description))
    SeqIO.write(chequences,os.path.join(out_path,'CORE_Genes','gene_sequences',nmechlseqin),'fasta')
#######################              single Sequence file         #########################
print('Process H.\n**Making Sequences for MSA**\n\n\n')
os.mkdir(os.path.join(out_path,'CLUSTALW'))
mseq=os.listdir(os.path.join(out_path,'CORE_Genes','gene_sequences'))
for x in range(len(mseq)):
    xmseq=[]
    join_seq=''
    for xseq in SeqIO.parse(os.path.join(out_path,'CORE_Genes','gene_sequences',mseq[x]),'fasta'):
        join_seq=join_seq+xseq.seq
    xmseq.append('>'+mseq[x].rstrip('.fasta')+'\n'+str(join_seq))
    with open(os.path.join(out_path,'CORE_Genes','msa_seq',mseq[x]),'a') as dic1:
        dic1.write('>'+mseq[x].rstrip('.fasta')+'\n')
        dic1.write(str(join_seq))
    for seqloeach in xmseq:
        with open (os.path.join(out_path,'CLUSTALW','all_seq_for_msa.fasta'),'a') as dic2:
            dic2.write(seqloeach+'\n')
######################        Doing MSA     ###############################################
print('Process I.\n**Running ClustalW**\n\n\n')
align='clustalw -INFILE='+os.path.join(out_path,'CLUSTALW','all_seq_for_msa.fasta')+' -QUICKTREE -OUTFILE='+os.path.join(out_path,'CLUSTALW','all_seq_for_msa.aln.fasta')+' -OUTPUT=FASTA -KTUPLE=2 -TOPDIAGS=5 -WINDOW=5 -PAIRGAP=3 -SCORE=PERCENT -GAPOPEN=8.0 > '+os.path.join(out_path,'CLUSTALW','CLUSTALW-logfile')
os.system(align)
#####################   Making the Newick File for Core Genes  ############################
print('Process J.\n**Making Newick file for the core genes\n\n\n')
os.mkdir(os.path.join(out_path,'FAST-TREE'))
NWK=('fasttree < '+os.path.join(out_path,'CLUSTALW','all_seq_for_msa.aln.fasta')+' > '+os.path.join(out_path,'FAST-TREE','core_genes.nwk'))
os.system(NWK)
#####################   Renaming homogenome files   ##################
prenames=os.listdir(os.path.join(out_path,'homogenome'))
for bleach in prenames:
    renames=bleach.split('.fasta')[0]+'.txt'
    os.rename(os.path.join(out_path,'homogenome',bleach),os.path.join(out_path,'homogenome',renames))
#####################   Non-Recombinant Starts  ######################
print('Process K.\n**Non-Recombinant Gene Analysis starts**\n\n\n')
os.mkdir(os.path.join(out_path,'Non-Recombinant'))
os.mkdir(os.path.join(out_path,'Non-Recombinant','Log-Files'))
os.mkdir(os.path.join(out_path,'Non-Recombinant','Single-Liners_FASTA'))
os.mkdir(os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA'))
funitems=os.listdir(os.path.join(out_path,'CORE_Genes','gene_sequences'))
for funeach in (funitems):
    for funseq in SeqIO.parse(os.path.join(out_path,'CORE_Genes','gene_sequences',funeach),'fasta'):
        with open(os.path.join(out_path,'Non-Recombinant','Single-Liners_FASTA',funeach),'a') as funf:
            funf.write('>'+funseq.id+'\n'+str(funseq.seq)+'\n')
itemsfun=os.listdir(os.path.join(out_path,'Non-Recombinant','Single-Liners_FASTA'))
for eachfun in (itemsfun):
    with open(os.path.join(out_path,'Non-Recombinant','Single-Liners_FASTA',eachfun),'r') as funnf:
        rdfunnf=funnf.readlines()
    for fnum in range(0,len(rdfunnf),2):
        if(fnum!=0):
            nfnum=fnum/2
        else:
            nfnum=fnum
        with open(os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA','rec'+str(int(nfnum))+'.fasta'),'a') as fnumb:
            fnumb.write(rdfunnf[fnum]+rdfunnf[fnum+1])
fintems=os.listdir(os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA'))
os.mkdir(os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA_aligned'))
##########################################################################
print('Process L.\n**Running clustalw for each set of core genes in all the organisms**\n\n\n')
for fineach in fintems:
    align='clustalw -INFILE='+os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA',fineach)+' -QUICKTREE -OUTFILE='+os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA_aligned',fineach.split('.fasta')[0]+'.aln.fasta')+' -OUTPUT=FASTA -KTUPLE=2 -TOPDIAGS=5 -WINDOW=5 -PAIRGAP=3 -SCORE=PERCENT -GAPOPEN=8.0 > '+os.path.join(out_path,'Non-Recombinant','Log-Files',fineach.split('.fasta')[0]+'.log')
    os.system(align)
##############################
print('Process M.\n**Phipack analysis for recombination\n\n\n')
os.mkdir(os.path.join(out_path,'Non-Recombinant','Phipack'))
philist=os.listdir(os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA_aligned'))
for phieach in philist:
    phicmd='/usr/bin/phipack-phi -f '+os.path.join(out_path,'Non-Recombinant','Gene_set_FASTA_aligned',phieach)+' -t A -w 10 -v -g i > '+os.path.join(out_path,'Non-Recombinant','Phipack',phieach.split('.aln.fasta')[0])
    #os.system(phicmd)
    try:
        os.system(phicmd)
    except ERROR:
        continue
##############################
print('Process N.\n**Creating List of Recombinant and Non-Recombinant genes: Tables,Fasta files')
val_nrc=[]
Val_nrc={}
names_gene=[]
os.mkdir(os.path.join(out_path,'Non-Recombinant','Recombinant-data'))
nrclist=os.listdir(os.path.join(out_path,'Non-Recombinant','Phipack'))
for nreach in range(len(nrclist)):
    with open(os.path.join(out_path,'Non-Recombinant','Phipack',nrclist[nreach]),'r') as nrcf:
        nrcfl=nrcf.readlines()
    for inr in range(len(nrcfl)):
        if('PHI (Normal):' in nrcfl[inr]):
            elnrcfl=nrcfl[inr]
            val_nrc.append(nrclist[nreach]+':'+(elnrcfl.rstrip('\n').split(':'))[1].strip(' '))
for rnum in range(len(val_nrc)):
    nvar='rec'+str(rnum)
    for start in val_nrc:
        if(nvar == start.split(':')[0]):
            try:
                Val_nrc.update({start.split(':')[0]:float(start.split(':')[1])})
            except ValueError:
                Val_nrc.update({start.split(':')[0]:'--'})
            with open(os.path.join(out_path,'Non-Recombinant','Recombinant-data','sorted_List_Phi_Values.txt'),'a') as wrfn:
                wrfn.write(start+'\n')
                break
for key, value in Val_nrc.items():
    if(value != '--' and float(value)<0.05):
        names_gene.append(int(key[3:]))
        with open(os.path.join(out_path,'Non-Recombinant','Recombinant-data','all_core_genes_having_recombination.txt'),'a') as rrfn:
            rrfn.write(key+':'+str(value)+'\n')
#############################################################################################################################################################
os.mkdir(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-data'))
os.mkdir(os.path.join(out_path,'Non-Recombinant','Recombinant-data','All-organisms-Reco-data'))
for rcneach in lseqin:
    with open(os.path.join(out_path,'CORE_Genes','gene_names',rcneach),'r') as rcnf:
        rd_rcnf=rcnf.readlines()
    for dist in range(len(names_gene)):
        del rd_rcnf[dist]
        with open(os.path.join(out_path,'Non-Recombinant','Recombinant-data','All-organisms-Reco-data',rcneach),'a') as wnrdf:
            wnrdf.write(rd_rcnf[dist])
    with open(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-data',rcneach),'a') as wrdf:
        for freach in range(len(rd_rcnf)):
            wrdf.write((rd_rcnf[freach]))
os.mkdir(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-Sequences'))
listseq=os.listdir(os.path.join(out_path,'Non-Recombinant','Single-Liners_FASTA'))
listseqdat=os.listdir(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-data'))
for oseach in range(len(listseq)):
    sweq=[]
    for soseach in range(len(listseqdat)):
        with open(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-data',listseqdat[soseach]),'r') as filesos:
            rd_filesos=filesos.readlines()
        for tin in range(len(rd_filesos)):
            hedseq=rd_filesos[tin].split('\n')[0]
            for esseq in SeqIO.parse(os.path.join(out_path,'Non-Recombinant','Single-Liners_FASTA',listseq[oseach]),'fasta'):
                if(hedseq == esseq.id):
                    sweq.append(SeqRecord(Seq(str(esseq.seq)),id=esseq.id,description=esseq.description))
    SeqIO.write(sweq,os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-Sequences',listseq[oseach]),'fasta')
os.mkdir(os.path.join(out_path,'Non-Recombinant','nR-CLUSTALW'))
msseq=os.listdir(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-Sequences'))
for x in range(len(msseq)):
    xmsseq=[]
    join_sseq=''
    for xsseq in SeqIO.parse(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-Sequences',msseq[x]),'fasta'):
        join_sseq=join_sseq+xsseq.seq
    xmsseq.append('>'+msseq[x].rstrip('.fasta')+'\n'+str(join_sseq))
    with open(os.path.join(out_path,'Non-Recombinant','All-organisms-nReco-Sequences',msseq[x]),'a') as dic1:
        dic1.write('>'+msseq[x].rstrip('.fasta')+'\n')
        dic1.write(str(join_sseq))
    for sseqloeach in xmsseq:
        with open (os.path.join(out_path,'Non-Recombinant','nR-CLUSTALW','nR_all_seq_for_msa.fasta'),'a') as dic2:
            dic2.write(sseqloeach+'\n')
###################################################################################
print('Process O.\n**Running clustalw for Non-Recombinant genes\n\n\n')
align='clustalw -INFILE='+os.path.join(out_path,'Non-Recombinant','nR-CLUSTALW','nR_all_seq_for_msa.fasta')+' -QUICKTREE -OUTFILE='+os.path.join(out_path,'Non-Recombinant','nR-CLUSTALW','nR_all_seq_for_msa.aln.fasta')+' -OUTPUT=FASTA -KTUPLE=2 -TOPDIAGS=5 -WINDOW=5 -PAIRGAP=3 -SCORE=PERCENT -GAPOPEN=8.0 > '+os.path.join(out_path,'Non-Recombinant','nR-CLUSTALW','CLUSTALW-logfile')
os.system(align)
###################################################################################
print('Process P.\n**Making Newick file for the Non-Recombinant genes\n\n\n')
os.mkdir(os.path.join(out_path,'Non-Recombinant','FAST-TREE'))
NWK=('fasttree < '+os.path.join(out_path,'Non-Recombinant','nR-CLUSTALW','nR_all_seq_for_msa.aln.fasta')+' > '+os.path.join(out_path,'Non-Recombinant','FAST-TREE','nReco_genes.nwk'))
os.system(NWK)
