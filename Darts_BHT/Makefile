install:
	R CMD build Darts
	R CMD INSTALL Darts_0.1.2.tar.gz
	if [ ! -d $(HOME)/.darts ]; then mkdir $(HOME)/.darts/; fi
	cp -f ./bin/Darts_BHT.mp $(HOME)/.darts/Darts_BHT
	chmod +x $(HOME)/.darts/Darts_BHT

uninstall:
	R CMD REMOVE Darts
	rm $(HOME)/.darts/Darts_BHT
