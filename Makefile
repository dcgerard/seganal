nc = 12
rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

.PHONY: all
all: null_sims

.PHONY: huge_null
huge_null: ./output/null_pvals.RDS

./output/null_pvals.RDS : ./code/null_huge.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

.PHONY: null_sims
null_sims: ./output/nullsims/null_pvalues.RDS

./output/nullsims/null_pvalues.RDS: ./code/null_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout
