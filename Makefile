
all:

################################################################
# LD blocks
step1: ldblocks/EUR

ldblocks/EUR: ldblocks
	ln -s ./nygcresearch-ldetect-data-ac125e47bf7f/EUR $@

ldblocks: ac125e47bf7f.zip
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	unzip $< -d $@

ac125e47bf7f.zip:
	wget https://bitbucket.org/nygcresearch/ldetect-data/get/ac125e47bf7f.zip

################################################################
# break down LD block by LD block
PHENOTYPES := $(shell ls -1 UKBB/*.sumstats.gz | xargs -I file basename file .sumstats.gz)
TEMPDIR := /broad/hptmp/ypp/ukbb/tempdata
LD := ldblocks/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)
ldlist := $(shell seq 1 $(NBLOCKS))

step2: jobs/step2.txt.gz

jobs/step2.txt.gz: $(foreach pheno, $(PHENOTYPES), jobs/step2/$(pheno)-jobs)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p log/
	@cat $^ | gzip > $@
	qsub -P compbio_lab -o log/step2.log -binding "linear:1" -cwd -V -l h_vmem=16g -l h_rt=17200 -b y -j y -N UKBB_DATA -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

# % = $(pheno)
jobs/step2/%-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@echo "./make.data-ldblock.R UKBB/$*.sumstats.gz $(LD) $(TEMPDIR)" > $@

################################################################
# run FQTL independently
CHR := $(shell seq 1 22)
STEP3-JOBS := $(foreach chr, $(CHR), $(shell ls -1 $(TEMPDIR)/$(chr)/ | sed 's/\///' 2> /dev/null | awk '{ print "jobs/step3/$(chr)/" $$1 "-jobs" }'))

step3: jobs/step3.txt.gz

jobs/step3.txt.gz: $(STEP3-JOBS)
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@mkdir -p log/
	@cat $^ | gzip > $@
	qsub -P compbio_lab -o log/step3.log -binding "linear:1" -cwd -V -l h_vmem=16g -l h_rt=14400 -b y -j y -N UKBB_RUN -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

jobs/step3-resubmit.txt.gz:
	@zcat jobs/step3.txt.gz | awk 'system("[ ! -f " $$NF ".zscore.gz ]") = 0' | gzip > $@
	[ $$(zcat $@ | wc -l) -gt 0 ] && qsub -P compbio_lab -o log/step3.log -binding "linear:1" -cwd -V -l h_vmem=16g -l h_rt=86400 -b y -j y -N UKBB_RE_RUN -t 1-$$(zcat $@ | wc -l) ./run_rscript.sh $@

# % = $(chr)/$(ld)
jobs/step3/%-jobs:
	@[ -d $(dir $@) ] || mkdir -p $(dir $@)
	@echo ./make.run-fqtl.R $(TEMPDIR)/$* 1KG/plink/chr$(shell echo $* | awk -F'/' '{ print $$1 }') 45 result/fqtl/$*/fqtl > $@
