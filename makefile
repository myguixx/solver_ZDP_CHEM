# Please chech the manual for make
# https://www.gnu.org/software/make/manual/make.html

CC    := gfortran
OPTION = -O0 -Wall -Wextra
DEBAG =

TARGET = zdp_chem
SKN_O = src/senkin/cklib.o src/senkin/dasac.o src/senkin/driv.o src/senkin/senkin.o
ZDP_O = src/zdplaskin/dvode_f90_m.o src/zdplaskin/zdplaskin_m.o src/zdplaskin/calc_rop.o
ZDP_M = dvode_f90_m.mod zdplaskin.mod

$(TARGET): $(SKN_O)
	$(CC) $(OPTION) $(DEBAG) $^ -o $@

%.o: %.f
	$(CC) $(OPTION) $(DEBAG) -c $< -o $@

.PHONY: clean
clean :
	rm $(TARGET) $(SKN_O)



# .SUFFIXES: .F90

# $(TARGET): $(OBJS)
# 	$(CC) $(OPTION) $(DEBAG) $^ src/zdplaskin/bolsig_x86_64.so -o $@

# %.o: %.F90
# 	$(CC) $(OPTION) $(DEBAG) -o $@ -c $<

# %.mod: %.F90 %.o
#     @:

# .PHONY: clean
# clean :
# 	rm $(TARGET) $(OBJS) $(MODS)