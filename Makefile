
all: node_modules node-dependencies env

# Python

env: ppqm.git
	mamba env create -f ./environment.yml -p ./env --quiet
	cd ppqm.git && ../env/bin/python -m pip install -e .
	./env/bin/python -m pip install -e .

ppqm.git:
	git clone https://github.com/ppqm/ppqm ppqm.git

# Node

node_modules:
	npm install

node-build: src/regiosqm_api/public
	npm run build

src/regiosqm_api/public:
	cd src/regiosqm_api && ln -s ../../public public

node-dependencies: public/chemdoodleweb public/fontawesome public/rdkit

public/chemdoodleweb:
	bash scripts/setup_chemdoodle.sh

public/fontawesome:
	bash scripts/setup_fontawesome.sh

public/rdkit:
	bash scripts/setup_rdkit.sh

node-dev:
	npm run dev

node-format:
	npm run pretty

# Development

start_flask:
	python -m python_api

# Interaction

start_jupyter:
	jupyter lab
