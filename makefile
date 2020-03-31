# Please chech the manual for make
# https://www.gnu.org/software/make/manual/make.html

OPTION = -O0 -Wall -Wextra
DEBAG = -g

CC    := gfortran
ALL_O := $(patsubst src/zdplaskin/%.F90,src/zdplaskin/%.o,$(wildcard src/zdplaskin/*.F90))

zdplaskine: $(ALL_O)
	$(CC) $(OPTION) $(DEBAG) $^ src/zdplaskin/bolsig_x86_64.so -o $@
%.o: %.F90
	$(CC) $(OPTION) $(DEBAG) -c $< -o $@

# senkine: $(ALL_O)
# 	$(CC) $(OPTION) $(DEBAG) $^ -o $@
# %.o: %.f
# 	$(CC) $(OPTION) $(DEBAG) -c $< -o $@

.PHONY: clean
clean :
	rm senkine $(ALL_O)