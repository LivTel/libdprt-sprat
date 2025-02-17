# $Header: /home/cjm/cvs/libdprt-sprat/c/Makefile,v 1.2 2025-02-06 16:36:24 cjm Exp $

include	../../../Makefile.common
include	../../Makefile.common
include	../Makefile.common

LIBNAME			= $(LIBDPRT_HOME)_$(SPRAT_HOME)
INCDIR 			= $(LIBDPRT_SPRAT_SRC_HOME)/include
CCDIMAGERINCDIR		= $(LIBDPRT_SRC_HOME)/ccd_dprt
LTFILENAMESINCDIR 	= $(LT_SRC_HOME)/lt_filenames
OBJECTINCDIR		= $(LIBDPRT_SRC_HOME)/object/include
DOCSDIR 		= $(LIBDPRT_SPRAT_DOC_HOME)/cdocs
BINDIR			= $(LIBDPRT_SPRAT_BIN_HOME)/c/$(HOSTTYPE)
CFLAGS 			= -g $(CCHECKFLAG) -I$(INCDIR) -I$(CFITSIOINCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR) \
			-I$(JNIGENERALINCDIR) -I$(CCDIMAGERINCDIR) -I$(LTFILENAMESINCDIR) -I$(OBJECTINCDIR) \
			-L$(LT_LIB_HOME)
LINTFLAGS 		= -I$(INCDIR) -I$(JNIINCDIR) -I$(JNIMDINCDIR)
DOCFLAGS 		= -static
SRCS 			= dprt.c ngat_dprt_sprat_DpRtLibrary.c
HEADERS			= $(SRCS:%.c=%.h)
OBJS			= $(SRCS:%.c=$(BINDIR)/%.o)
DOCS 			= $(SRCS:%.c=$(DOCSDIR)/%.html)
LIBS			= -lcfitsio -ldprt_object -ldprt_libfits -llt_filenames -ldprt_jni_general -lsprat_ccd_dprt 

top: shared docs

shared: $(INCDIR)/ngat_dprt_sprat_DpRtLibrary.h $(LT_LIB_HOME)/$(LIBNAME).so

$(INCDIR)/ngat_dprt_sprat_DpRtLibrary.h: $(DPRT_BIN_HOME)/java/ngat/dprt/sprat/DpRtLibrary.class
	$(JAVAH) -jni -d $(INCDIR) -classpath $(DPRT_BIN_HOME)/java:$(CLASSPATH) ngat.dprt.sprat.DpRtLibrary
	$(TOUCH) $@

$(LT_LIB_HOME)/$(LIBNAME).so: $(OBJS)
	$(CC) $(CCSHAREDFLAG) $(CFLAGS) $(OBJS) -o $@ $(TIMELIB) $(LIBS) $(SOCKETLIB)

$(BINDIR)/%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

static: $(LT_LIB_HOME)/$(LIBNAME).a

$(LT_LIB_HOME)/$(LIBNAME).a: $(OBJS)
	ar rcv $@ $?

docs: $(DOCS)

$(DOCS) : $(SRCS)
	-$(CDOC) -d $(DOCSDIR) -h $(INCDIR) $(DOCFLAGS) $(SRCS)

# dont checkout ngat_dprt_o_DpRtLibrary.h - it is a machine built header
checkout:
	$(CO) $(CO_OPTIONS) $(SRCS)
	cd $(INCDIR); $(CO) $(CO_OPTIONS) dprt.h;

# dont checkin ngat_dprt_o_DpRtLibrary.h - it is a machine built header
checkin:
	-$(CI) $(CI_OPTIONS) $(SRCS)
	-(cd $(INCDIR); $(CI) $(CI_OPTIONS) dprt.h;)

staticdepend:
	makedepend $(MAKEDEPENDFLAGS) -p$(BINDIR)/ -- $(CFLAGS)  -- $(SRCS)

depend:
	makedepend $(MAKEDEPENDFLAGS) -p$(BINDIR)/ -- $(CFLAGS)  -- $(SRCS)

lint:
	$(LINT)	$(LINTFLAGS) $(SRCS)

clean:
	-$(RM) $(RM_OPTIONS) $(OBJS) $(LT_LIB_HOME)/$(LIBNAME).so $(LT_LIB_HOME)/$(LIBNAME).a 
	-$(RM) $(RM_OPTIONS) $(INCDIR)/ngat_dprt_sprat_DpRtLibrary.h $(TIDY_OPTIONS)

tidy:
	-$(RM)  $(RM_OPTIONS) $(TIDY_OPTIONS)
	-cd $(INCDIR); $(RM) $(RM_OPTIONS) $(TIDY_OPTIONS);

backup: tidy
	-$(RM) $(RM_OPTIONS) $(OBJS)

#
# $Log: not supported by cvs2svn $
# Revision 1.1  2014/09/03 14:07:41  cjm
# Initial revision
#
#
