GLMs := 12_ 

toclean := $(shell ls [0-9][0-9]_*/Makefile |sed s/Makefile//g)

all : 22_ 24_ 09_ 11_ 12_ 14_ 15_ 


clean:
	rm -f *.Rdata

cleanall: clean
	echo $(toclean)
	for i in $(toclean); do $(MAKE) -C $$i clean; done

15_: 11_
	@$(MAKE) -C 15_Binom_Dup-CpDup_Div-Rel1Mdl

14_: 11_
	@$(MAKE) -C 14_Binom_Polym_Div-Rel-1Mdl

#13_: 11_
	#@$(MAKE) -C 13_Binom2Steps-polym_Diverg-Related

12_: 11_
	@$(MAKE) -C 12_Multinom-polym_Diverg-Related

11_: 03_ 23_
	@$(MAKE) -C 11_GLMTruncFetchData

10_: 03_
	@$(MAKE) -C 10_GeneralStatsA

09_: 08_
	@$(MAKE) -C 09_NJFromRF

08_: 04_
	@$(MAKE) -C 08_RandomForest

04_: 03_
	@$(MAKE) -C 04_Stability-RF

05_: 03_
	@$(MAKE) -C 05_NJ

23_: 03_
	@$(MAKE) -C 23_TargetAndBaitCount

03_: 02_
	@$(MAKE) -C 03_ProcessingSegmentation

02_: 01_
	@$(MAKE) -C 02_PreProcessing
	
01_: 
	@$(MAKE) -C 01_GetData

22_:
	@$(MAKE) -C 22_SequencingStats

24_: 02_
	@$(MAKE) -C 24_ShowPreProcessingEffic


