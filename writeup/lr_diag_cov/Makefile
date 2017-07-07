vem_derivation.pdf: vem_derivation.tex vem_bib.bib
	pdflatex vem_derivation.tex
	pdflatex vem_derivation.tex
	bibtex vem_derivation.aux
	pdflatex vem_derivation.tex
	rm -f *.aux *.bbl *.blg *.log *.out
