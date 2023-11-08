##############################################################################################
##############################################################################################
####### script for preprocessing reads, requires sickle and sabre for demultiplexing #########
####### useful if reads are on paired end format in two different files   	     #########
####### Import the full restriction sequence, in my case CATG                        ######### 
#######	trim reads that recreate the cutting site in the sequence		     #########
####### remove also the reads that does not have the initial cut site		     #########
#######                                                                              #########
####### Require installation of vcftools, samtools, bwa, R and r-dependencies (ReQON)#########
####### USAGE: python pipeline.py forward_seq rev_seq genome_file                    #########
##############################################################################################
##############################################################################################

from __future__ import print_function
import sys, os, math
from commands import getoutput

def remove_chimera_pe(fwd, rev, site):
	#### open in parallel the fastq file
	#### trim chimera with CATG, remove reads shorter than 30bp,
	#### keep reads starting with CATG after barcode
	#### return name filtered files 
	n=len(site)
	fwd_in=open(fwd)
	fwd_out=open('%s_good.fq' % fwd.split('.')[0], 'w')
	f_file='%s_good.fq' % fwd.split('.')[0]
	rev_in=open(rev)
	rev_out=open('%s_good.fq' % rev.split('.')[0], 'w')
	r_file='%s_good.fq' % rev.split('.')[0]
	while True:
		fname, rname= fwd_in.readline(), rev_in.readline()
		if fname=='':
			break
		fseq, rseq= fwd_in.readline(), rev_in.readline()
		fplus, rplus= fwd_in.readline(), rev_in.readline()
		fqual, rqual= fwd_in.readline(), rev_in.readline()
		fseq, fqual=clip_chimera(fseq, fqual, site)##### first clip for RE site
		rseq, rqual=clip_chimera(rseq, rqual, site)
		fseq, fqual=qual_trim(fseq, fqual)###### then trim for quality
		rseq, rqual=qual_trim(rseq, rqual)
		minlen=min(len(fseq.strip()), len(rseq.strip()))
		if minlen>=30:##### filter reads for lenght
			if fseq[5:5+n]==site and rseq[5:5+n]==site: ##### filter good reads that starts with the RE site
				if fseq.count('N') < 3 and rseq.count('N') < 3:#remove n
					fwd_out.write(fname+fseq+fplus+fqual)
					rev_out.write(rname+rseq+rplus+rqual)
	fwd_in.close()
	fwd_out.close()
	rev_in.close()
	rev_out.close()
	return f_file, r_file
######### 
def qual_trim(seq, qual):
	n_qual=[ord(x)-33 for x in qual[:-1]]
	for i in range(len(n_qual)-4): ###### trim with a sliding window of 5, the 4 at the end is for avoid error
		clip=n_qual[i:i+5]
		m=float(sum(clip)/5)
		if m <= 20.0:
			qual=qual[:i]+'\n'
			seq=seq[:i]+'\n'
			break
	return seq, qual
#############
def clip_chimera(s, qual, site):
	if site in s[len(site)+5:]: #adding also len barcode modify if necessary
		p=s.index(site, len(site)+5)
		s=s[:p]+'\n'
		qual=qual[:p] +'\n'
	return s, qual


####### STARTING ANALYSIS
fwd_seq=sys.argv[1]
rev_seq=sys.argv[2]
genome=sys.argv[3] 
pattern='CATG'
print('Preprocessing files, quality trimming, trimming chimeras chimeras and removing sequence without the initial cutting site,  more than 2 Ns and shorter than 30bp after trimming')
fwd=remove_chimera_pe(fwd_seq, pattern)
###### call the barcode file barcode_sabre.txt and use the name [barcode|bc|lib]XXX_[f|r].fq where X could be the barcode sequence or a number, f and r stands for forward and reverse
###### the prefix is important for the alignment subroutes
###### later

print('starting demultiplexing samples')
dem='sabre pe -f s_1_filter.txt -r s_2_filter.txt -b barcode_sabre.txt -u unknow_f -w unknow_r -c -m 1' ###### BARCODE FILE NAME!!!!!!!!
print(dem)
os.system(dem)

#########################################################################################################
####### part for automate alignment using bwa aligner and the index already demultiplexed
####### specify the prefix of the file name. The file should be named by common_prefix+whatever+[f|r].fq. 
####### in PE sequencing name f for forward and r for reverse
####### PAIRED END ALIGNMENT
############################################################################################################

print('start genome indexing')
os.system('bwa index -a bwtsw %s' % genome) #INSERT GENOME folder

print('start alignment')
prefix='BC' ##### !!!!PREFIX of FILES!!!!!
fwd_files=getoutput('ls %s*f.fq' % prefix).split('\n')
fwd_files.sort()
rev_files=getoutput('ls %s*r.fq' % prefix).split('\n')
rev_files.sort()
for i in range(len(fwd_files)):
	print(fwd_files[i], rev_files[i])
	sai_fwd_out=fwd_files[i].split('.')[0] + '_aln.sai'
	aln1='bwa aln -n 0.1 -o 2 -e 6 -l 20 -k 3 -t 4 %s %s > %s' % (genome, fwd_files[i], sai_fwd_out) ####### GENOME FOLDER
	print(aln1)
	os.system(aln1)
	sai_rev_out=rev_files[i].split('.')[0]+'_aln.sai'
	aln2='bwa aln -n 0.1 -o 2 -e 6 -l 20 -k 3 -t 4 %s %s > %s' % (genome, rev_files[i], sai_rev_out)  ###### GENOME FOLDER
	print(aln2)
	os.system(aln2)
	sam_out='_'.join(sai_rev_out.split('_')[:2]) + '_pe_aln.sam'
	cmd2='bwa sampe %s %s %s %s %s > %s' % (genome, sai_fwd_out, sai_rev_out, fwd_files[i], rev_files[i], sam_out) ###### GENOME FOLDER
	print(cmd2)
	os.system(cmd2)

#################################################################################################################
#### Genome indexing
os.system('samtools index %s' % genome) 

#################################################################################################
##### filter reads with a mapping quality below 10
##### script for post process reads and call snp
##### reads should be in sam format
##### 

aln=getoutput('ls *sam').split('\n')
for i in aln:
	name=i.split('.')[0] + '_sort'
	cmd='samtools view -bShq 10 %s | samtools sort - %s' % (i, name)
	print(cmd)
	os.system(cmd)
######################## creates a bunch of files that finish with sort.bam
####################################################################################################
#########	QUALITY SCORE RECALIBRATION
#########	call an r script (in the same folder)
######### 	remember to install r libraries (ReQON and Rsamtools)
#################################################################

print('quality recalibration')
os.system('R CMD BATCH qual_recal.R')
	
##### subscript for indexing bam files

bam_files=getoutput('ls *recalibrated.bam').split('\n')
for i in bam_files:
	cmd='samtools index %s' % i
	print(cmd)
	os.system(cmd)

###########################################################################
##############
############## SNP calling and filtering subroutine
print('starting snps calling')
os.system('samtools mpileup -uDgf %s *recalibrated.bam|bcftools view -vcg -d 0.1 - > all_snp_raw.vcf' % genome) ### GENOME FOLDER!!!!!!
os.system('vcftools --vcf all_snp_raw.vcf --maf 0.05 --minQ 10 --min-meanDP 5 --max-meanDP 1000 --recode --out filtered_snp') ### 
####### export some basic statistics
#vcftools --vcf filtered_snp.recode.vcf --out filter_snp_stats --depth --site-depth --site-mean-depth --SNPdensity 100000 --hist-indel-len --het
