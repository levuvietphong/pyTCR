all: paper.pdf

paper.pdf: paper.md paper.bib
	pandoc --citeproc --bibliography paper.bib  paper.md \
  	-o paper.pdf --pdf-engine=xelatex

clean:
	rm paper.pdf

.PHONY: clean
