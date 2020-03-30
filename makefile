COMPILE = gfortran

COMLIB1 =

COMOTT1 = -O0 -Wall -Wextra

COMCONV = 

FILES = \
	senkin.o \
	driv.o \
	cklib.o \
	dasac.o



.SUFFIXES: .f .o

.f.o:
	$(COMPILE) $(COMOTT1) $(COMCONV) -c $<

senkine: $(FILES)
	$(COMPILE) $(COMOTT1) $(COMCONV) $(FILES) -o senkine $(COMLIB1)

senkin.o: senkin.f
driv.o: driv.f
cklib.o: cklib.f
dasac.o: dasac.f