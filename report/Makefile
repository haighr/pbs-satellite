SHELL=/bin/bash

figs=figure.pdf ETS_polygons.pdf ETS_polygons_clipped.pdf coast_bc_poly.pdf coast_bc_joinPoly.pdf sst_coast_bc.pdf plot_ts_1.pdf sst_a.pdf sst_b.pdf sst_c.pdf trivial.pdf
rdfiles=$(wildcard ../PBSsatellite/man/*.Rd)

.PHONY: all clean mrproper

all: report.pdf

report.pdf: report.tex man.tex ${figs}
	latexmk -quiet -xelatex report.tex || { cat report.log; exit 1; }

man.tex: mkman ${rdfiles}
	./mkman

%.pdf: %.R
	R --vanilla < $<

clean:
	latexmk -C
	${RM} report.bbl

mrproper: clean
	$(RM) ${figs}
