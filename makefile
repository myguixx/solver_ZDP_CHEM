# Please chech the manual for make
# https://www.gnu.org/software/make/manual/make.html

OPTION = -O0 -Wall -Wextra
DEBAG = -g

CC    := gfortran
ALL_O := $(patsubst src/senkin/%.f,src/senkin/%.o,$(wildcard src/senkin/*.f))

senkine: $(ALL_O)
	$(CC) $(OPTION) $(DEBAG) $^ -o $@

%.o: %.f
	$(CC) $(OPTION) $(DEBAG) -c $< -o $@

.PHONY: clean
clean :
	rm senkine $(ALL_O)