object = ShallowIceApproximations

pdf:
	pdflatex $(object).tex

bib:
	pdflatex $(object).tex
	bibtex $(object)
	pdflatex $(object).tex
	pdflatex $(object).tex

clean:
	rm -f $(object).aux $(object).log $(object).out $(object).blg $(object).bbl $(object).toc

open:
	open $(object).pdf

