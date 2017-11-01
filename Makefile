
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
LD := ldblocks/EUR/fourier_ls-all.bed
NBLOCKS := $(shell cat $(LD) 2> /dev/null | tail -n+2 | wc -l)
ldlist := $(shell seq 1 $(NBLOCKS))

step2: $(foreach pheno, $(PHENOTYPES), $(foreach ld, $(ldlist), stat/$(ld)/$(pheno)-stat.txt.gz))
	echo $(NBLOCKS)

# % = $(ld)/$(pheno)
stat/%-stat.txt.gz: stat/%-small.bed $(LD)
	cat $<
	zcat UKBB/$(shell echo $* | awk -F'/' '{ print $$2 }').sumstats.gz | tail -n+2 | awk -F'\t' '{ printf "chr" $$2 FS ($$3 - 1) FS $$3 FS $$1; for(j=4; j<=NF; ++j) printf FS $$j; printf "\n"; }' | bedtools intersect -a stdin -b $< -wa | gzip > $@

stat/%-small.bed:
	[ -d $(dir $@) ] || mkdir -p $(dir $@)
	cat $(LD) | tail -n+2 | head -n $(shell echo $* | awk -F'/' '{ print $$1 }') | tail -n1 | awk -F'\t' '{ gsub(/ /,"",$$0); print $$1 FS $$2 FS $$3 }' > $@

