# Please chech the manual for make
# https://www.gnu.org/software/make/manual/make.html

CC    := gfortran
OPTION = -O0 -Wall -Wextra -J module -I module
DEBAG =

TARGET = zdp_chem
SKN_O = src/senkin/cklib.o src/senkin/dasac.o src/senkin/driv.o \
        src/senkin/senkin.o src/senkin/zdplib.o src/senkin/chemkin_m.o
SKN_M = module/chemkin.mod
ZDP_O = src/zdplaskin/dvode_f90_m.o src/zdplaskin/zdplaskin_m.o
ZDP_M = module/dvode_f90_m.mod module/zdplaskin.mod

$(TARGET): $(ZDP_O) $(SKN_O) 
	$(CC) $(OPTION) $(DEBAG) $^ src/zdplaskin/bolsig_x86_64.so -o $@

%.o: %.f
	$(CC) $(OPTION) $(DEBAG) -c $< -o $@

.SUFFIXES: .F90

%.o: %.F90
	$(CC) $(OPTION) $(DEBAG) -c $< -o $@ 

%.mod: %.F90 %.o
    @:

.PHONY: clean
clean :
	rm $(TARGET) $(SKN_O) $(ZDP_O) $(ZDP_M)