all: jags.pdf openbugs.pdf

clean:
	rm -rf  *.toc *.bbl *.blg *.aux *.out *.log *.spl *tikzDictionary *.md figure/ cache/ *.dvi *.4* *.css *.idv *.lg *.tmp *.xref *.gz my*.bib *.nav *.xml *.vrb *.snm 

empty:
	rm *.tex *.pdf *.html *.docx *.odt

%.pdf: %.tex
	xelatex -interaction=nonstopmode $<;
	xelatex -interaction=nonstopmode $<;
	R --slave -e "tools::compactPDF('$@', gs_quality='ebook')"