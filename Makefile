MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

.PHONY: all pbF12 clean

all: pbF12

pbF12:
	@make -s -C pb
	@ln -sf pb/pbF12 pbF12

clean:
	@rm -f pbF12
	@make clean -s -C pb

distclean: clean
