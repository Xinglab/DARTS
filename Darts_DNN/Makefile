uninstall :
	rm -rf $(HOME)/.darts/sequenceFeature/*; \
	rm -rf $(HOME)/.darts/trainedParam/*; \
	pip uninstall Darts_DNN

install :
	if [ ! -f "$(HOME)/.darts/sequenceFeature/ENCODE_sequenceFeature_absmax_normalized.h5" ]; then \
		echo "fetching cis features.."; \
		mkdir -p $(HOME)/.darts/sequenceFeature/ $(HOME)/.darts/sequenceFeature; \
		wget -O $(HOME)/.darts/sequenceFeature/ENCODE_sequenceFeature_absmax_normalized.h5 "https://drive.google.com/open?id=1U-uHVBfr6cWFIIKaINnXbWjCh7khVx-I"; \
		fi; \
	if [ ! -f "$(HOME)/.darts/trainedParam/encode_plus_roadmap/best_testing_model.h5" ]; then \
		echo "fetching trained model params.."; \
		mkdir -p $(HOME)/.darts/trainedParam/encode_plus_roadmap; \
		wget -O $(HOME)/.darts/trainedParam/encode_plus_roadmap/best_testing_model.h5 "https://drive.google.com/open?id=12F5QurlD1Mqj4blcKX941kz3KAM5GQu5"; \
		fi; \
	python setup.py install