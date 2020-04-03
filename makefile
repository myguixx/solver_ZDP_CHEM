# Please chech the manual for make
# https://www.gnu.org/software/make/manual/make.html

CC    := gfortran
OPTION = -O0 -Wall -Wextra
DEBAG =

TARGET = zdplaskine
OBJS = src/zdplaskin/dvode_f90_m.o src/zdplaskin/zdplaskin_m.o src/zdplaskin/calc_rop.o
MODS = dvode_f90_m.mod zdplaskin.mod

.SUFFIXES: .F90

$(TARGET): $(OBJS)
	$(CC) $(OPTION) $(DEBAG) $^ src/zdplaskin/bolsig_x86_64.so -o $@

%.o: %.F90
	$(CC) $(OPTION) $(DEBAG) -o $@ -c $<

%.mod: %.F90 %.o
    @:

.PHONY: clean
clean :
	rm $(TARGET) $(OBJS) $(MODS)