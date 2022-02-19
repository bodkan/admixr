.PHONY: build vignettes docs

version := $(shell less DESCRIPTION | grep 'Version' | sed 's/Version: \(.*\)$$/\1/')

docs:
	R -e 'devtools::document()'
	R -e 'pkgdown::build_reference(examples = FALSE)'
	R -e 'pkgdown::build_reference_index()'

website:
	rm -rf docs/
	R -e 'devtools::install(upgrade = "never")'
	R -e 'knitr::knit("README.Rmd", output = "README.md")'
	R -e 'devtools::document()'
	R -e 'pkgdown::build_site(examples = FALSE)'
#	git restore docs/pkgdown.yml

build: $(pkg)

check: $(pkg)
	cd build; R CMD check --as-cran $(notdir $<)

winrel: README.md
	R -e 'devtools::check_win_release()'

windev: README.md
	R -e 'devtools::check_win_devel()'

winold: README.md
	R -e 'devtools::check_win_oldrelease()'

rhub: README.md
	R -e 'rhub::check_for_cran()'

$(pkg): README.md
	R -e 'devtools::document()'
	mkdir -p build; cd build; R CMD build --log ../../admixr

README.md: README.Rmd
	R -e 'devtools::install(upgrade = "never")'
	R -e 'knitr::knit("README.Rmd", output = "README.md")'

clean:
	rm -rf build
