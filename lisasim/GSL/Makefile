PLATFORM=$(shell uname -s)-$(shell uname -r)

OPTIONS= -Wall -O3

BUILD= build/$(PLATFORM)

CFILES= mt.c ranlux.c rng.c

HFILES= gsl_errno.h gsl_rng.h gsl_types.h

OFILES= $(foreach cfile,$(CFILES),$(BUILD)/$(cfile:.c=.o))

libs: $(BUILD) $(OFILES)

$(BUILD)/%.o: %.c $(HFILES)
	gcc $(OPTIONS) -c $< -o $@

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -rf $(BUILD)/*.o
