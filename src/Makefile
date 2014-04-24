GLMs:= 12_ 

toclean:= $(shell ls [0-9][0-9]_*/Makefile |sed s/Makefile//g)

all: part1 part2

part1: 22_ 24_ 25_ 27_ 14_ 15_ 16_ 17_ 18_ 19_ 07_ 09_ 28_

part2: 08_ 26_	# long process

clean:
	rm -f *.Rdata

cleanall: clean
	echo $(toclean)
	for i in $(toclean); do $(MAKE) -C $$i clean; done

19_: 16_
	@$(MAKE) -C 19_Plot_fqcy_Div-Rel

18_: 15_
	@$(MAKE) -C 18_Plot_Dup-CpDup_Div-Rel

17_: 14_
	@$(MAKE) -C 17_Plot_Polym_Div-Rel

16_: 11_
	@$(MAKE) -C 16_Binom_fqcy_Div-Rel

15_: 11_
	@$(MAKE) -C 15_Binom_Dup-CpDup_Div-Rel

14_: 11_
	@$(MAKE) -C 14_Binom_Polym_Div-Rel

11_: 03_ 23_
	@$(MAKE) -C 11_GLMTruncFetchData

26_: 06_ 08_
	@$(MAKE) -C 26_Stats_Best_Subtarg_RF

28_: 03_ 06_
	@$(MAKE) -C 28_Compute_Vst

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

27_: 03_
	@$(MAKE) -C 27_Count_Singletons

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

24_: 02_
	@$(MAKE) -C 24_ShowPreProcessingEffic

22_:
	@$(MAKE) -C 22_SequencingStats



