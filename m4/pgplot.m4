# TRM_LIB_PGPLOT
#---------------
# macro to try to get linking to PGPLOT with gcc right the problem being that
# it needs to know the correct libraries to cope with the underlying fortran. 
# Basically works through a few possible library combinations. It assumes
# that both the header file and the library are in standard locations and 
# will fail if they are not. At its end it has defineS HAVE_PGPLOT if it has 
# worked, and PGPLOT_LIBS are the libraries to use. It simply bombs out at
# the moment if it does not work.

AC_DEFUN([TRM_LIB_PGPLOT],
[

  AC_CHECK_HEADERS([cpgplot.h], [], AC_MSG_ERROR([PGPLOT header not found; please fix.]))

  usual_libs="-lcpgplot -lpgplot -L/usr/X11R6/lib -lX11 -lpng -lz"

  save_libs="$LIBS"

  AC_MSG_CHECKING([pgplot linking assuming simplest set of libraries])

  LIBS="$usual_libs"

  AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "cpgplot.h"]], [cpgopen("/null");])],
    [AC_MSG_RESULT([yes]) 
     PGPLOT_LIBS="$LIBS" 
     pgplot_found=yes], 
    [AC_MSG_RESULT([no]) 
      pgplot_found=no])

  if test $pgplot_found = no; then
    AC_MSG_CHECKING([pgplot linking assuming g77 or g95])

    LIBS="$usual_libs -lg2c"

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "cpgplot.h"]], [cpgopen("/null");])],
      [AC_MSG_RESULT([yes]) 
       PGPLOT_LIBS="$LIBS" 
       pgplot_found=yes], 
      [AC_MSG_RESULT([no]) 
       pgplot_found=no])
  fi

  if test $pgplot_found = no; then
    AC_MSG_CHECKING([pgplot linking assuming gfortran])

    LIBS="$usual_libs -lgfortran"

    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "cpgplot.h"]], [cpgopen("/null");])],
      [AC_MSG_RESULT([yes]) 
       PGPLOT_LIBS="$LIBS" 
       pgplot_found=yes], 
      [AC_MSG_RESULT([no]) 
       pgplot_found=no])
  fi

  if test $pgplot_found = yes; then
     HAVE_PGPLOT=1
     AC_SUBST([PGPLOT_LIBS])
     LIBS="$LIBS $save_libs"
  else
     AC_MSG_ERROR([could not find or work out libraries for PGPLOT])
     LIBS="$save_libs"
  fi

  unset save_libs
  unset usual_libs
  unset pgplot_found

])
