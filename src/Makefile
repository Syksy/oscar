# Makefile for the proximal bundle method DBDC with COX's PROPORTIONAL HAZARDS MODEL

#FF = f95
FF = gfortran 
#FFLAGS = 
FFLAGS = -fbounds-check -Wall
OPEN = 
#OPEN =  -fopenmp 
#RM = del   #in windows
RM = rm    #in unix/linux/mac

all: tdbdc

tdbdc: constants.o bundle1.o bundle2.o functions.o dbdc.o cox_dbdc.o tdbdc.o plqdf1.o 
	$(FF) -o tdbdc $(FFLAGS) $(OPEN) constants.o bundle1.o bundle2.o functions.o dbdc.o cox_dbdc.o tdbdc.o plqdf1.o 

constants.mod: constants.o constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95
	
constants.o: constants.f95
	$(FF) -c $(FFLAGS) $(OPEN) constants.f95

bundle1.mod: constants.mod bundle1.o bundle1.f95 
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95
	
bundle1.o: constants.mod bundle1.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle1.f95 

bundle2.mod: constants.mod bundle2.o bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 
	
bundle2.o: constants.mod bundle2.f95
	$(FF) -c $(FFLAGS) $(OPEN) bundle2.f95 

functions.mod: constants.mod functions.o functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 
	
functions.o: constants.mod functions.f95
	$(FF) -c $(FFLAGS) $(OPEN) functions.f95 

dbdc.mod: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.o dbdc.f95 
	$(FF) -c $(FFLAGS) $(OPEN) dbdc.f95	 
	
dbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) dbdc.f95 
	
cox_dbdc.mod: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.mod cox_dbdc.o cox_dbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) cox_dbdc.f95 
	
cox_dbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.mod cox_dbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) cox_dbdc.f95 	
	
tdbdc.o: constants.mod bundle1.mod bundle2.mod functions.mod dbdc.mod cox_dbdc.mod tdbdc.f95
	$(FF) -c $(FFLAGS) $(OPEN) tdbdc.f95 
	
plqdf1.o: plqdf1.f
	$(FF) -c $(FFLAGS) $(OPEN) plqdf1.f 
	
clean:	
	$(RM) tdbdc constants.mod constants.o bundle1.mod bundle1.o bundle2.mod bundle2.o functions.mod functions.o dbdc.mod dbdc.o cox_dbdc.mod cox_dbdc.o tdbdc.o plqdf1.o  
	echo Clean done