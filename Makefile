SCRIPTS=$(shell find postproc runners utils -type f) 

PREFIX?=$(HOME)/local/chromahd

all: install

install: 
	@mkdir -p $(PREFIX)
	@for script in $(SCRIPTS); do \
		ln -sfn $(PWD)/$$script $(PREFIX)/; \
	done
	@echo installed to $(PREFIX)
