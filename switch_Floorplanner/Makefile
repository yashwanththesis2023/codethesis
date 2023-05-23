#---------------------------------DIR------------------------------------
GUROBIDIR= /pc2/groups/hpc-prf-lgrsfpga/yashwanth/Gurobi/gurobi1000/linux64
CPLEXDIR      = /pc2/groups/hpc-prf-lgrsfpga/yashwanth/cplex/cplex
CONCERTDIR    = /pc2/groups/hpc-prf-lgrsfpga/yashwanth/cplex/concert


#---------------------------------LIB--------------------------------------

CPLEXLIBDIR   = $(CPLEXDIR)/lib/x86-64_linux/static_pic
CONCERTLIBDIR = $(CONCERTDIR)/lib/x86-64_linux/static_pic
GUROBILIBDIR= $(GUROBIDIR)/lib

#---------------------------------LINK-------------------------------------
CPLEXINCDIR   = $(CPLEXDIR)/include
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINC 	= -I$(CONCERTINCDIR) -I$(CPLEXINCDIR) 
CPLEXDIRS  	= -L$(CPLEXLIBDIR) -L$(CONCERTLIBDIR) 
CPLEXFLAGS 	= -lconcert -lilocplex -lcplex -lm -pthread 
GUROBILIB   = -L$(GUROBILIBDIR) -lgurobi_c++ -lgurobi100 

#-------------------------------LPSOLVE------------------------------------------
lpsolveDIR=-L/usr/lib/lp_solve -llpsolve55

#------------------------------COMPILER-OPTION----------------------------------
CPP      = g++ -m64 -O -fPIC -fno-strict-aliasing -fexceptions -DNDEBUG -DIL_STD 
#--------------------------------------------------------------------------



CODE=switchFloorplanner
CODE_FILE= newfpgaplanner

#--------------------------------------------------------------------------
ifeq ($(ALL_USE), 1)
CPLEX_USE=1
GUROBI_USE=1
LPSOLVE_USE=1
endif

ifeq ($(CPLEX_USE), 1)  
CPLEX =-DCPLEX_USE
else
CPLEX=
CPLEXINC=
CPLEXDIRS=
CPLEXFLAGS= 	
endif

ifeq ($(GUROBI_USE), 1)  
GUROBI =-DGUROBI_USE 
else
GUROBI =
GUROBILIB=
endif

ifeq ($(LPSOLVE_USE), 1)  
LPSOLVE =-DLPSOLVE_USE 
else
LPSOLVE =
lpsolveDIR=
endif


all:	
	$(CPP) $(LPSOLVE) $(GUROBI) $(CPLEX) $(CPLEXINC) $(CPLEXDIRS) -o $(CODE_FILE) $(CODE_FILE).cpp $(CPLEXFLAGS) $(GUROBILIB) -ldl $(lpsolveDIR)

clean:
	rm  $(CODE)