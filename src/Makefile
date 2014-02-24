GLMs := 12_ 

toclean := $(shell ls [0-9][0-9]_*/Makefile |sed s/Makefile//g)

all : 22_ 24_ 25_ 11_ 12_ 14_ 15_ 16_ 17_ 18_ 07_ 08_ 09_


clean:
	rm -f *.Rdata

cleanall: clean
	echo $(toclean)
	for i in $(toclean); do $(MAKE) -C $$i clean; done

18_: 11_
	@$(MAKE) -C 18_Plot_Dup-CpDup_Div-Rel

17_: 11_
	@$(MAKE) -C 17_Plot_Polym_Div-Rel

16_: 11_
	@$(MAKE) -C 16_Binom_fqcy_Div-Rel

15_: 11_
	@$(MAKE) -C 15_Binom_Dup-CpDup_Div-Rel

14_: 11_
	@$(MAKE) -C 14_Binom_Polym_Div-Rel

#13_: 11_
	#@$(MAKE) -C 13_Binom2Steps-polym_Diverg-Related

12_: 11_
	@$(MAKE) -C 12_Multinom-polym_Diverg-Related

11_: 03_ 23_
	@$(MAKE) -C 11_GLMTruncFetchData

#10_: 03_
	#@$(MAKE) -C 10_GeneralStatsA

09_: 03_ 06_
	@$(MAKE) -C 09_Count_CNV

08_: 06_
	@$(MAKE) -C 08_GeneFamImptceTest

07_: 06_
	@$(MAKE) -C 07_NJFromRF

06_: 05_
	@$(MAKE) -C 06_RandomForest

05_: 03_ 04_
	@$(MAKE) -C 05_Stability-RF

04_: 03_
	@$(MAKE) -C 04_NJ

25_: 03_
	@$(MAKE) -C 25_CNVplot_alongChr

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


