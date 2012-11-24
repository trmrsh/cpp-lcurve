## Process this file with automake to generate Makefile.in
##

EXTRA_DIST = docgen.in document.pl Doxyfile.in style.css guide.html example_model_file example_data_file

install-data-local:
	$(INSTALL) -d html
	sed -e 's%ZzZzZ%$(VERSION)%' docgen.in > docgen.pl
	perl docgen.pl  html
	cp style.css    html/.
	cp guide.html   html/.
	cp example_model_file  html/.
	cp example_data_file  html/.
	sed -e 's%ZzZzZ%html/api%' Doxyfile.in > Doxyfile
	doxygen Doxyfile || :
	chmod 0755 html
	@for file in `find html | sed -e 's% %\\\%'`; do \
	nfile=`echo $$file | sed -e 's%\\\%\\ %'`; \
	if test -d "$$nfile"; then \
	  chmod 0755 "$$nfile"; \
	else \
	  chmod 0644 "$$nfile"; \
	fi; \
	done
	-rm -r     $(prefix)/html/$(PACKAGE)
	-mkdir -p  $(prefix)/html
	cp -r html $(prefix)/html/$(PACKAGE)

uninstall-local:
	rm -r $(prefix)/html

clean-local:
	-rm -r html
