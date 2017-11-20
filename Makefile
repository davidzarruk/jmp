#---------------------
#Housekeeping        #
#---------------------

CORES := $(shell bash -c 'read -p "Initial or final SS: " pwd; echo $$pwd') 

#1. C++ compiler
NVCC = g++
NVCC2 = nvcc

#2. Flags
CFLAGS=-fopenmp -Wall  

#3. Name of output
OUTPUT=mainDZ

#4. Files to run in main:
SOURCE=CUDA_model.cu

#5. Home directory
HOME=/home/davidza/
HOMEAWS=/home/

#6. Dynamic libraries 
DLIB=-lgsl -lgslcblas 

#7. Set number of threads	
export OMP_NUM_THREADS=32

#-------------------------------------
#3. Include directories:
INCDIR3=$(HOME)install/include/
INCDIR3AWS=$(HOMEAWS)install/include/

#Generating the INC for all the paths
INC=$(INCDIR3)
INC_PARAMS=$(foreach d, $(INC), -I$d)

INCAWS=$(INCDIR3AWS)
INC_PARAMSAWS=$(foreach d, $(INCAWS), -I$d)
#-------------------------------------

#4. Lnlopt library 
LNLOPT=-L$(HOME)install/lib/ -lm -lnlopt 
LNLOPTAWS=-L$(HOMEAWS)install/lib/ -lm -lnlopt 


myprog:
	$(NVCC) $(SOURCE) $(INC_PARAMS) $(LNLOPT)  $(DLIB) $(CFLAGS)  -o $(OUTPUT)
	./$(OUTPUT)

cuda:
	$(NVCC2) $(SOURCE) $(INC_PARAMS) $(LNLOPT)  $(DLIB) -Wno-deprecated-gpu-targets -o $(OUTPUT)
	./$(OUTPUT) $(CORES)

cudaws:
	$(NVCC2) $(SOURCE) $(INC_PARAMSAWS) $(LNLOPTAWS)  $(DLIB) -o $(OUTPUT)
	./$(OUTPUT) $(CORES)

cudasilent:
	$(NVCC2) $(SOURCE) $(INC_PARAMS) $(LNLOPT)  $(DLIB) -o $(OUTPUT)
	nohup ./$(OUTPUT) $(CORES) > cuda_output.out&

nloptinst:
	wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz -P $(HOME)
	tar  -zxvf $(HOME)nlopt-2.4.2.tar.gz -C$(HOME)
	mkdir $(HOME)install
	cd $(HOME)nlopt-2.4.2 && ./configure --prefix=$(HOME)install && make && make install

nloptaws:
	wget http://ab-initio.mit.edu/nlopt/nlopt-2.4.2.tar.gz -P $(HOMEAWS)
	tar  -zxvf $(HOMEAWS)nlopt-2.4.2.tar.gz -C$(HOMEAWS)
	mkdir $(HOMEAWS)install
	cd $(HOMEAWS)nlopt-2.4.2 && ./configure --prefix=$(HOMEAWS)install && make && make install

execute: myprog
	./$(OUTPUT)
	scp davidza@tesla.ssc.upenn.edu:/home/davidza/heterogeneous/Model-V1/CUDA/V2/results.txt "/home/david/Dropbox/Documents/Doctorado/Research/Heterogeneous households/Model-V1/CUDA/V2/"

executesilent: myprog
	nohup ./$(OUTPUT) > console_output.out&


clean: 
	rm $(OUTPUT)

