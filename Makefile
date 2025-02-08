rexec = R CMD BATCH --no-save --no-restore
rout = ./output/rout

.PHONY: all
all: huge_null

.PHONY: huge_null
huge_null: ./output/null_pvals.RDS

../output/null_pvals.RDS : ./code/null_huge.R
	mkdir -p $(rout)
	mkdir -p $(@D)
	$(rexec) $< $(rout)/$(basename $(<F)).Rout
