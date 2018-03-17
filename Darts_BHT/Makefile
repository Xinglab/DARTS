install:
	R CMD INSTALL Darts_0.1.1.tar.gz
	if [ ! -d $(HOME)/.darts ]; then mkdir $(HOME)/.darts/; fi
	cp ./bin/Darts_BHT $(HOME)/.darts/
	chmod +x $(HOME)/.darts/Darts_BHT

uninstall:
	R CMD REMOVE Darts
	rm $(HOME)/.darts/Darts_BHT
