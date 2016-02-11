TESTS = check_ptv check_epi 
check_PROGRAMS = check_ptv check_epi 

check_ptv_SOURCES = check_ptv.c $(top_builddir)/src_c/ptv.h
check_ptv_CFLAGS = @CHECK_CFLAGS@ 
check_ptv_LDADD = $(top_builddir)/src_c/libptv.la @CHECK_LIBS@ -lm -loptv

check_epi_SOURCES = check_epi.c $(top_builddir)/src_c/calibration.h
check_epi_CFLAGS = @CHECK_CFLAGS@
check_epi_LDADD = $(top_builddir)/src_c/libptv.la @CHECK_LIBS@ -lm -loptv

