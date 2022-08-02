.PHONY: build vignettes docs

version := $(shell less DESCRIPTION | grep 'Version' | sed 's/Version: \(.*\)$$/\1/')
pkg := build/admixr_$(version).tar.gz

docs:
	rm -rf docs/reference
	R -e 'devtools::document()'
	R -e 'pkgdown::build_reference()'
	R -e 'pkgdown::build_reference_index()'
	R -e 'pkgdown::build_news()'

website:
	rm -rf docs/
	R -e 'devtools::install(upgrade = "never")'
	R -e 'knitr::knit("README.Rmd", output = "README.md")'
	R -e 'devtools::document()'
	R -e 'pkgdown::build_site()'

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
	R -e 'rhub::check_for_cran(platforms = c("macos-highsierra-release", "macos-highsierra-release-cran", "macos-m1-bigsur-release"))'

$(pkg): README.md
	R -e 'devtools::document()'
	mkdir -p build; cd build; R CMD build --log ../../slendr

README.md: README.Rmd $(logo)
	R -e 'devtools::install(upgrade = "never")'
	R -e 'knitr::knit("README.Rmd", output = "README.md")'

example_data:
	Rscript generate_examples.R

clean:
	rm -rf build
