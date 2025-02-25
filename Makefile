nc = 12
rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

## null sims output with outliers
nullplots = ./output/nullsims/null_t1e.pdf \
            ./output/nullsims/qq_example.pdf

## null sims output, no outliers
nood_nullplots = ./output/nood_nullsims/nood_null_t1e.pdf \
                 ./output/nood_nullsims/nood_qq_example.pdf

## alt sims output
alt_output = ./output/altsims/alt_tab.txt

.PHONY: all
all: null_sims nood_null_sims alt_sims dr_null_sims

## Null simulations under fewer scenarios but massive sample size
## Just to verify the asymptotics are correct.
## Not part of main manuscript
.PHONY: huge_null
huge_null: ./output/null_pvals.RDS

./output/null_pvals.RDS : ./code/null_huge.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

## Null simulations with outliers ----
.PHONY: null_sims
null_sims: $(nullplots)

$(nullplots): ./code/null_plots.R ./output/nullsims/null_pvalues.RDS
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

./output/nullsims/null_pvalues.RDS: ./code/null_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Null simulations, no outliers ----
.PHONY: nood_null_sims
nood_null_sims: $(nood_nullplots)

$(nood_nullplots): ./code/nood_null_plots.R ./output/nood_nullsims/nood_null_pvalues.RDS
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

./output/nood_nullsims/nood_null_pvalues.RDS: ./code/nood_null_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Double Reduction Null Simulations ----
.PHONY: dr_null_sims
dr_null_sims: ./output/dr_nullsims/dr_null_pvalues.RDS

./output/dr_nullsims/dr_null_pvalues.RDS: ./code/dr_null_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout

## Alternative simulations ----
.PHONY: alt_sims
alt_sims: $(alt_output)

$(alt_output): ./code/alt_plots.R ./output/altsims/alt_pvalues.RDS
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout

./output/altsims/alt_pvalues.RDS: ./code/alt_sims.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) '--args nc=$(nc)' $< $(rout)/$(basename $(<F)).Rout
