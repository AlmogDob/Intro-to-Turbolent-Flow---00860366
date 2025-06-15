file: Turb_hw2.tex
	@pdflatex Turb_hw2.tex > latex_temp.txt
	@makeindex Turb_hw2.nlo -s nomencl.ist -o Turb_hw2.nls > latex_temp.txt
	@pdflatex Turb_hw2.tex > latex_temp.txt
	@pdflatex Turb_hw2.tex > latex_temp.txt
	@cp -r *.log logs/.
	@rm *.aux *.ilg *.lof *.nlo *.nls *.out *.toc latex_temp.txt *.log *.fls *.fdb_latexmk *.lol *.bbl *.blg

clean_after_vscode:
	@test ! $(wildcard *.log) || cp -r *.log logs/.
	@rm *.fdb_latexmk *.fls *.synctex.gz *.aux *.bbl *.log *.blg *.toc *.out *.lof *.nlo *.lol

test:
	@test ! $(wildcard *.pdf) || echo File does exist