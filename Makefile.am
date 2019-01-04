## Process this file with automake to generate Makefile.in

ACLOCAL_AMFLAGS = -I m4

SUBDIRS = src doc include

export: dist
	scp $(PACKAGE)-$(VERSION).tar.gz $(WEB_SERVE):$(WEB_PATH)/software/.
	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software ; \rm -f $(PACKAGE).tar.gz ; ln -s $(PACKAGE)-$(VERSION).tar.gz $(PACKAGE).tar.gz; \rm -fr $(PACKAGE); mkdir $(PACKAGE)"
	scp -r doc/html $(WEB_SERVE):$(WEB_PATH)/software/$(PACKAGE)/.
#	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- lcurve version --><td>.*</td>%<\!-- lcurve version --><td>$(VERSION)</td>%' index.html > junk.html; mv junk.html index.html"
#	ssh $(WEB_SERVE) "cd $(WEB_PATH)/software; sed -e 's%<\!-- lcurve date --><td>.*</td>%<\!-- lcurve date --><td>$(DATE)</td>%' index.html > junk.html; mv junk.html index.html"

ALIASES = Lcurve

.PHONY : $(ALIASES)

DATE    = $(shell date)

progdir = @bindir@/@PACKAGE@

prog_DATA = $(ALIASES)

$(ALIASES):
	@echo '#' > $(ALIASES)
	@echo '# This file defines aliases for the $(PACKAGE) v$(VERSION) software.' >> $(ALIASES)
	@echo '# Source it to define them.' >> $(ALIASES)
	@echo '#' >> $(ALIASES)
	@echo '# this to allow both bash and csh to work' >> $(ALIASES)
	@echo 'test "$$?BASH_VERSION" = "0" || eval '\''alias() { command alias "$$1=$$2"; }'\' >> $(ALIASES)
	@echo '#' >> $(ALIASES)
	@echo 'alias lprofile  $(progdir)/lprofile'       >> $(ALIASES)
	@echo 'alias levmarq   $(progdir)/levmarq'        >> $(ALIASES)
	@echo 'alias lroche    $(progdir)/lroche'         >> $(ALIASES)
	@echo 'alias lroches   $(progdir)/lroches'         >> $(ALIASES)
	@echo 'alias picture   $(progdir)/picture'        >> $(ALIASES)
	@echo 'alias simplex   $(progdir)/simplex'        >> $(ALIASES)
	@echo 'alias visualise $(progdir)/visualise'      >> $(ALIASES)
	@echo ' ' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Welcome to $(PACKAGE) version $(VERSION), built $(DATE)"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Commands available are: lprofile, levmarq, lroche, lroches, picture, simplex and visualise"' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
	@echo 'echo "Point a browser at ${prefix}/html/$(PACKAGE)/index.html for help."' >> $(ALIASES)
	@echo 'echo " "' >> $(ALIASES)
