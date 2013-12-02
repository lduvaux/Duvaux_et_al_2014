GLMs := 12_ 

#~ #
#~ 
#~ all : 09_ $(GLMs)
#~ 
#~ 
#~ 12_: 11_
	#~ @$(MAKE) -C 12_GLM00_TEMPLATE
#~ 
#~ 13_: 11_

#all : 21_ 22_09_ 11_ $(GLMs)
all : 22_ 24_ 09_ 11_ 12_

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


clean:
	rm -f *.Rdata
