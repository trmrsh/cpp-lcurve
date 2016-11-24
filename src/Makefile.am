## Process this file with automake to generate Makefile.in
##
## This is the file that must be edited if you are changing anything in the source directory

## Programs

progdir = @bindir@/@PACKAGE@

prog_PROGRAMS   = lprofile visualise picture simplex lroche levmarq rotprof lroches

lprofile_SOURCES = lprofile.cc numface.cc set_star_grid.cc set_star_emission.cc light_curve_comp.cc comp_light.cc trm_lcurve.cc \
set_star_continuum.cc re_scale.cc set_disc_grid.cc set_bright_spot_grid.cc set_disc_continuum.cc star_eclipse.cc

levmarq_SOURCES = levmarq.cc trm_lcurve.cc set_star_grid.cc comp_light.cc set_star_continuum.cc numface.cc light_curve_comp.cc re_scale.cc \
set_disc_grid.cc set_disc_continuum.cc set_bright_spot_grid.cc star_eclipse.cc
levmarq_CPPFLAGS = -fopenmp
levmarq_LDFLAGS  = -fopenmp

lroche_SOURCES = lroche.cc numface.cc set_star_grid.cc set_star_continuum.cc comp_light.cc trm_lcurve.cc light_curve_comp.cc re_scale.cc \
set_disc_grid.cc set_disc_continuum.cc set_bright_spot_grid.cc star_eclipse.cc
lroche_CPPFLAGS = -fopenmp
lroche_LDFLAGS  = -fopenmp

lroches_SOURCES = lroche.cc numface.cc set_star_grid.cc set_star_continuum.cc comp_light.cc trm_lcurve.cc light_curve_comp.cc re_scale.cc \
set_disc_grid.cc set_disc_continuum.cc set_bright_spot_grid.cc star_eclipse.cc

picture_SOURCES = picture.cc pos_disc.cc disc_eclipse.cc

rotprof_SOURCES = rotprof.cc trm_lcurve.cc set_star_grid.cc light_curve_comp.cc re_scale.cc numface.cc set_star_continuum.cc \
comp_light.cc set_disc_grid.cc star_eclipse.cc set_disc_continuum.cc set_bright_spot_grid.cc

simplex_SOURCES = simplex.cc trm_lcurve.cc set_star_grid.cc comp_light.cc set_star_continuum.cc numface.cc light_curve_comp.cc re_scale.cc \
set_disc_grid.cc set_disc_continuum.cc set_bright_spot_grid.cc star_eclipse.cc
simplex_CPPFLAGS = -fopenmp
simplex_LDFLAGS  = -fopenmp

visualise_SOURCES = visualise.cc numface.cc set_star_grid.cc set_disc_grid.cc trm_lcurve.cc light_curve_comp.cc re_scale.cc set_star_continuum.cc \
comp_light.cc set_disc_continuum.cc set_bright_spot_grid.cc star_eclipse.cc

INCLUDES = -I../include

DATE = $(shell date)

install-data-hook: 
	echo "This is $(PACKAGE)-$(VERSION), built on $(DATE)" > $(progdir)/VERSION



