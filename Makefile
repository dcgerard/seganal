nc = 12
rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

.PHONY: all
all: null_sims nood_null_sims nood_alt_sims

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

## No outliers null simulations ----
.PHONY: nood_null_sims
nood_null_sims: ./output/nood_nullsims/nood_qq_example.pdf

./output/nood_nullsims/nood_qq_example.pdf: ./code/nood_null_plots.R ./output/nood_nullsims/nood_null_pvalues.RDS
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

./output/nood_nullsims/nood_null_pvalues.RDS: ./code/nood_null_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## No outlier alternative simulations ----

.PHONY: nood_alt_sims
nood_alt_sims: ./output/nood_altsims/nood_alt_pvalues.RDS

./output/nood_altsims/nood_alt_pvalues.RDS: ./code/nood_alt_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout
