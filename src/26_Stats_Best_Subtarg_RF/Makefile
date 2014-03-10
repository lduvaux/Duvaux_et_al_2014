
MYDIR := $(shell basename `pwd`)
OUT := ../$(MYDIR).Rdata

UTILSDEPS := $(shell ls ../utils/*.R)
LOCALDEPS := $(shell ls ./*.R | grep -v script.R)


$(OUT): script.R $(LOCALDEPS) $(UTILSDEPS)
	Rscript $< $@ 
		
clean:
	rm -f *.Rdata *.txt *.pdf *.log result_*.csv TableS*.csv $(OUT)
