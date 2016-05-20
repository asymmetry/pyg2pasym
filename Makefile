MYOS        := $(shell uname)
ARCH        := $(shell uname -m)
USER        := $(shell whoami)
MYHOST      := $(shell hostname -s)

.PHONY: all clean pb polint

all: pb polint

pb:
	@make -s -C pb
	@ln -sf pb/pbF12 pbF12

polint:
	@make -s -C polint
	@ln -sf polint/polel polel
	@ln -sf polint/polin polin

clean:
	@rm -f pbF12 polel polin
	@make clean -s -C pb
	@make clean -s -C polint

distclean: clean
