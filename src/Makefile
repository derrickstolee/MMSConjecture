#
# Makefile for MMSConjecture project.
#

TESTS	=	CStarTest.exe

TARGETS = 	mmsconj.exe					


OBJECTS = 	MMSBranchingManager.o  		\
			MMSImplementationCPLEX.o	\
			MMSImplementationGLPK.o	


CC 	= gcc
CXX = gcc

OPT = -O0 -g
# -O7
DEBUG =  
# -g 
# -static 
WARNINGS = -Wall


CFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)
CXXFLAGS 	= $(OPT) $(DEBUG) $(WARNINGS)
LFLAGS 		= $(OPT) $(DEBUG) $(WARNINGS)  -lglpk -lgmp
# Todo: put in -static when compiling for grid computing



TREESEARCHDIR = ../../TreeSearch
UTILSDIR = ../../Utilities/src
NAUTYDIR = ../../nauty
CLIQUERDIR = ../../cliquer
CPLEXDIR = ../../cplex
GLPKDIR = ../../glpk

INCLUDES =  -I$(TREESEARCHDIR) 	\
			-I$(NAUTYDIR)		\
			-I$(CLIQUERDIR)		\
			-I$(UTILSDIR)		\
			-I$(GLPKDIR)


TREESEARCHOBJS	= $(TREESEARCHDIR)/SearchManager.o
 
UTILSOBJS	= $(UTILSDIR)/Set.o				\
				$(UTILSDIR)/translation.o				\
				$(UTILSDIR)/TreeSet.o				 
				  
NAUTYOBJS     	= $(NAUTYDIR)/nauty.o 			\
				  $(NAUTYDIR)/nausparse.o		\
				  $(NAUTYDIR)/gtools.o			\
				  $(NAUTYDIR)/nautil.o			\
				  $(NAUTYDIR)/naugraph.o
				  
CLIQUEROBJS		= $(CLIQUERDIR)/cliquer.o		\
					$(CLIQUERDIR)/graph.o		\
					$(CLIQUERDIR)/reorder.o		
				  
LIBOBJS			= $(TREESEARCHOBJS) $(NAUTYOBJS) $(CLIQUEROBJS) $(UTILSOBJS)
			
			
.SUFFIXES: .c .cpp .o .obj .exe 

all: $(OBJECTS) $(TESTS) $(TARGETS)


# The default object compiler
.c.o: $<
	$(CC) $(CFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@
        
.cpp.o: $<
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(LIBS) -c $< -o $@
        
.cpp.exe: $< $(OBJECTS) 
	$(CXX) $(LFLAGS)			\
        	$(INCLUDES)	$(DEBUG)			\
        	$(LIBOBJS) $(LIBS)				\
        	`cat $@.objs`           		\
            $< -o $@
	cp $@ ../bin/
        
.c.exe: $< $(COBJECTS)
	$(CC) 	$(LFLAGS)			    \
        	$(INCLUDES)				\
        	$(NAUTYOBJS)  $(COBJECTS) $(LIBS)		\
            $< -o $@
	cp $@ ../bin/
        
clean:
	rm $(OBJECTS) $(TARGETS) $(TESTS)
	
cleanexe:
	rm $(TARGETS)

cleantest:
	rm $(TESTS)

clexe:
	rm $(TARGETS)
