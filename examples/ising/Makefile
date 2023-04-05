# Defining variables
FC:=gfortran
EXE:=ising
EXT:=.exe

all:			TYPE :=
debug:    TYPE := _Debug
catch:    TYPE := _Catch
coverage: TYPE := _Coverage

all:      OPT := -O3 -Wall -Werror
debug:    OPT := -O0 -g -fbacktrace -ffpe-trap=zero,overflow,underflow -Wall -Werror
catch:    OPT := -O0 -g -Wall -W -Wsurprising -Werror
coverage: OPT := -O0 -g --coverage -Wall -Werror

#-----------------------------------------------------------------------------------
#Directory Paths
#-----------------------------------------------------------------------------------
SIMANN := ../../src

#-----------------------------------------------------------------------------------
#Add source files as necessary
#If the files are not comiled using the generic rules, add commands under their
#rule declaration. Add these items to FILTER
#-----------------------------------------------------------------------------------

SRC := ising.f90

OBJ := $(SRC:.f90=.o)
MOD := $(OBJ:.o=.mod)

FILTER := \

OBJ_FILTER := $(FILTER:.f90=.o)
MOD_FILTER := $(FILTER:.f90=.mod)

#-----------------------------------------------------------------------------------
#Add simple compile SIMANN files
#-----------------------------------------------------------------------------------
SIMANN_SRC := $(SIMANN)/OpenFSAM.f90

SIMANN_OBJ := $(patsubst %.f90,%.o,$(notdir $(SIMANN_SRC)))

SIMANN_MOD := $(SIMANN_OBJ:.o=.mod)

#-----------------------------------------------------------------------------------
#Phony targets for cleaning and building
#-----------------------------------------------------------------------------------
.PHONY: all debug catch coverage clean reset

print-%  : ; @echo $* = $($*)

all: $(EXE)

debug: $(EXE)

catch: $(EXE)

coverage: $(EXE)

#Intended to clean up compilation artifacts but leave executable & coverage
clean:
	rm -f $(OBJ) $(SIMANN_OBJ) $(LIB_OBJ)
	rm -f $(MOD) $(SIMANN_MOD) $(LIB_MOD)
	rm -f $(COMP_DEP)

#Intended to reset directory to fresh state with no exe or artifacts
reset: clean
	rm -f *.gcno *.gcda *.o *.mod
	rm -f -r $(EXE)*.dSYM
	rm -f *.exe

#-----------------------------------------------------------------------------------
#Generics for source files
#-----------------------------------------------------------------------------------
$(filter-out $(OBJ_FILTER), $(OBJ)): %.o:	%.f90
	$(FC) -c $(OPT) $<

$(filter-out $(MOD_FILTER), $(MOD)):	%.mod:	%.f90
	$(FC) -c $(OPT) $<

$(EXE): $(OBJ) $(SIMANN_OBJ) $(LIB_OBJ)
	$(FC) -o $@$(TYPE)$(EXT) $(OPT) $(OBJ) $(SIMANN_OBJ) $(LIB_OBJ)

#-----------------------------------------------------------------------------------
#Generics for contrib files
#-----------------------------------------------------------------------------------
$(SIMANN_OBJ): %.o: $(filter %.f90, $(SIMANN_SRC))
	$(FC) -c $(OPT) $^
$(SIMANN_MOD):	%.mod:	$(filter %.f90, $(SIMANN_SRC))
	$(FC) -c $(OPT) $^

#-----------------------------------------------------------------------------------
#Dependency List
#Use [gfortran -M -cpp *.f90] repeatedly until clean compile to update rules below
#may need to make first to get OpenFSAM module
#-----------------------------------------------------------------------------------
ising.o: ising.f90 OpenFSAM.mod
