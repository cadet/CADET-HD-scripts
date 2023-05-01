SCRIPTS=$(shell find postproc runners utils -type f) 

prefix?=$(HOME)/local/chromahd

all: install

install: 
	mkdir -p $(prefix)
	cp $(SCRIPTS) $(prefix)
	echo installed to $(prefix)
