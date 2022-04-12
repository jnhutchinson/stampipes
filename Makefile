.PHONY: all docker docker-stampipes-dnase docker-stampipes-rna

all: docker

docker: docker-stampipes-dnase docker-stampipes-rna

docker-stampipes-dnase:
	docker build --target stampipes-dnase --tag stampipes-dnase .

docker-stampipes-rna:
	docker build --target stampipes-rna --tag stampipes-rna .
