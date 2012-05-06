filename=fcsr_and_bs

pdf: ps
	 ps2pdf ${filename}.ps

ps: dvi
	 dvips ${filename}.dvi

dvi:
	 latex ${filename}
	  bibtex ${filename}||true
		 latex ${filename}
		  latex ${filename}

clean:
	 rm -f ${filename}.ps ${filename}.pdf ${filename}.log ${filename}.aux ${filename}.out ${filename}.dvi ${filename}.bbl ${filename}.blg
