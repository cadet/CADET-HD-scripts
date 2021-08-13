SCRIPTS=$(shell find . -maxdepth 1 -type f -printf "%P\n") 

prefix?=$(HOME)/local

all: install

install: 
	@for script in $(SCRIPTS); do \
		ln -sfn $(PWD)/$$script $(prefix)/bin/$$script; \
	done

