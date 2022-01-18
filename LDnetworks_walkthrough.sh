# assuming you work on ls5.tacc.utexas.edu ...do you? I will be glad to add you to my account if you want
# all the lines starting with ls5_launcher_creator are specific to TACC, I think; but all you need is to execute the commands written into the file with a name corresponding to the -j argument in the call to ls5_launcher_creator

# replace matz@utexas.edu with your email

#---------- getting the latest version of my scripts

cd ~/bin 
git clone https://github.com/z0on/2bRAD_denovo.git
mv 2bRAD_denovo/* . 
rm -rf 2bRAD_denovo 
cd -

#---------- installing ngsLD

cd 
git clone https://github.com/fgvieira/ngsLD.git
cd ngsLD

nano Makefile
add -I${TACC_GSL_INC}  to CC and CXX macros (CFLAGS= ...);
and -L${TACC_GSL_LIB} to the 'LIB = ...' line.

module load gsl
export PKG_CONFIG_PATH=/opt/apps/intel18/gsl/2.2.1/lib/pkgconfig/
make
cp ngsLD ~/bin

ls *.bam > bams 

# Update March 15, 2020: new angsd commit by 
# nspope https://github.com/nspope/angsd.git --branch hetFilter
# has a bult-in filter -maxHetFreq [float] to -doHWE 
# (ADJUST -minInd argument to 70-80% of the total number of individuals in bams!): 

# filtering sites and getting genotypes (-doGeno 8, for posterior number of derived alleles):
FILTERS="-minMapQ 20 -minQ 20 -dosnpstat 1 -doHWE 1 -maxHetFreq 0.5 -sb_pval 1e-5 -hetbias_pval 1e-5 -skipTriallelic 1 -minInd 150 -snp_pval 1e-3 -minMaf 0.05"
TODO="-doMajorMinor 1 -doMaf 1 -doCounts 1 -makeMatrix 1 -doIBS 1 -doCov 1 -doGeno 8 -doPost 1"
echo "angsd -b bams -GL 1 $FILTERS $TODO -P 12 -out zz4" > zz4
ls5_launcher_creator.py -j zz4 -n zz4 -t 12:00:00 -e matz@utexas.edu -w 1 -a tagmap -q normal
sbatch zz4.slurm

# listing sites
zcat zz4.mafs.gz | cut -f 1,2 | tail -n +2 >allSites

# unzipping genotypes
zcat zz4.geno.gz > zz4.geno

# Thinning sites to 25 kb, 48 random replicates
>thi25
for R in `seq 1 48`; do
echo "thinner.pl vcf=allSites interval=25000 criterion=random >z25k_sites_$R">>thi25;
done
ls5_launcher_creator.py -j thi25 -n thi25 -t 1:00:00 -e matz@utexas.edu -w 12 -a tagmap 
sbatch thi25.slurm

# extracting sites from genotype files and launching ngsLD for each replicate: 
module load gsl
NB=`cat bams | wc -l`
>ld
for R in `seq 1 48`; do
echo "awk 'NR==FNR{a[\$1\$2]; next} \$1\$2 in a{print}' z25k_sites_$R zz4.geno > zz_$R.geno && NS=\`wc -l zz_$R.geno\` && cut -f 1,2 zz_$R.geno > zz_$R.geno.sites && gzip zz_$R.geno && ngsLD --geno zz_$R.geno.gz --probs 1 --n_ind $NB --n_sites \$NS --max_kb_dist 0 --pos zz_$R.geno.sites --out z25k_$R.LD --n_threads 12 ">> ld;
done
module load gsl
NB=`cat bams_clean | wc -l`
ls5_launcher_creator.py -j ld -n ld -t 24:00:00 -e matz@utexas.edu -w 4 -a tagmap -q normal
sbatch ld.slurm

# collecting sites that correlate with each other better than R2 = 0.1 :

>si
for R in `seq 1 48`; do
echo "awk '\$7>0.1' z25k_$R.LD | awk '{a=\$1 \"\\n\" \$2; print a}' >z25k_$R.LD.sites" >>si;
done
ls5_launcher_creator.py -j si -n si -t 1:00:00 -e matz@utexas.edu -w 12 -a tagmap 
sbatch si.slurm

cat *LD.sites | sort | uniq | perl -pe 's/:/\t/' | sort -k 1,1 -k 2,2n > LDsites

#------ Round two! now we only have "interesting" sites

# let's thin them to 10 kb in 48 replicates:
>thiL
for R in `seq 1 48`; do
echo "thinner.pl vcf=LDsites interval=10000 criterion=random >LD_sites_$R">>thiL;
done
ls5_launcher_creator.py -j thiL -n thiL -t 2:00:00 -e matz@utexas.edu -w 4 -a tagmap 
sbatch thiL.slurm

# how many sites we have in each batch?
wc -l LD_sites_1

# if the number is larger than 10000, thin them to larger distance (15000). Ideally aim for 6000-10000 sites.

# rerunning ngsLD
module load gsl
NB=`cat bams | wc -l`
>ld2
for R in `seq 1 48`; do
echo "awk 'NR==FNR{a[\$1\$2]; next} \$1\$2 in a{print}' LD_sites_$R zz4.geno > LD_$R.geno && NS=\`wc -l LD_$R.geno\` && cut -f 1,2 LD_$R.geno > LD_$R.geno.sites && gzip LD_$R.geno && ngsLD --geno LD_$R.geno.gz --probs 1 --n_ind $NB --n_sites \$NS --max_kb_dist 0 --pos LD_$R.geno.sites --out LD_$R.LD --n_threads 12 ">> ld2;
done
ls5_launcher_creator.py -j ld2 -n ld2 -t 24:00:00 -e matz@utexas.edu -w 4 -a tagmap -q normal
sbatch ld2.slurm

# making LD files a bit smaller:
>sl
for R in `seq 1 48`; do
echo "cut -f 1,2,7 LD_$R.LD > LD_$R.LD.slim">>sl;
done
ls5_launcher_creator.py -j sl -n sl -t 0:30:00 -e matz@utexas.edu -w 4 -a tagmap 
sbatch sl.slurm

# geno files for wgcna (-geno 8, posterior probs of all genotypes)
for R in `seq 1 48`; do
awk 'NR==FNR{a[$1$2]; next} $1$2 in a{print}' LD_sites_$R zz4.geno > zz4_$R.geno;
done

# making square matrices of correlations among sites:
>lm
for R in `seq 1 48`; do
echo "Rscript ld2matrix.R LD_$R.LD.slim">>lm;
done
ls5_launcher_creator.py -j lm -n lm -t 24:00:00 -e matz@utexas.edu -w 4 -a tagmap -q normal
sbatch lm.slurm

# --------------- next, WGCNA! I did not write the walkthrough for that yet, let me know when you get to this point!

# serial WGCNA
>ldss2
for I in `seq 1 48`;do echo "Rscript LD_WGCNA_auto.R LD_$I.LD.slim_matrix.RData bams zz4_$I.geno">>ldss2;done
ls5_launcher_creator.py -j ldss2 -n ldss2 -t 1:00:00 -e matz@utexas.edu -w 6 -a tagmap -q normal
sbatch ldss2.slurm


