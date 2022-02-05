
all: node_modules dependencies

node_modules:
	npm install

build:
	npm run build

dependencies: public/chemdoodleweb public/fontawesome public/rdkit

public/chemdoodleweb:
	bash scripts/setup_chemdoodle.sh

public/fontawesome:
	bash scripts/setup_fontawesome.sh

public/rdkit:
	bash scripts/setup_rdkit.sh

dev:
	npm run dev

format:
	npm run pretty
