PKGNAME   := YAM2
VERSION   := 0.2.2
ARCHIVE   := $(PKGNAME)-$(VERSION)
SRCDIR    := src
LIBDIR    := lib
CXXFLAGS  := -g -O2 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LIBS      := -lm
AR        := ar crs
MKDIR     := mkdir -p
CP        := cp -r
RM        := rm -f
UNAME     := $(shell uname -s)

LIB    	  := $(LIBDIR)/lib$(PKGNAME).a
LIBSRC 	  := $(wildcard $(SRCDIR)/*.cc)
LIBOBJ 	  := $(LIBSRC:.cc=.o)
EXE    	  := examples/m2 examples/m2cons
ifeq ($(UNAME), Darwin)
SHAREDLIB := $(LIBDIR)/lib$(PKGNAME).dylib
else
SHAREDLIB := $(LIBDIR)/lib$(PKGNAME).so
endif

DESTDIR   ?= /usr/local
HEADERS   := $(wildcard $(SRCDIR)/*.h)

# NLopt (https://nlopt.readthedocs.io/)
NLOPT ?= /usr
LIBS  += -L$(NLOPT)/lib -lnlopt -Wl,-rpath $(NLOPT)/lib

ROOT := $(shell command -v root-config 2> /dev/null)
ifdef ROOT
	CXXFLAGS += -I$(shell root-config --incdir) -DHAS_ROOT
endif

.PHONY: all lib install dist clean

all: $(LIB)

$(LIB): CXXFLAGS += -I$(NLOPT)/include -fPIC
$(LIB): $(LIBOBJ)
	$(MKDIR) $(LIBDIR)
	$(AR) $@ $^
	ranlib $@

ifeq ($(UNAME), Darwin)
lib: LDFLAGS += -dynamiclib -undefined dynamic_lookup
else
lib: LDFLAGS += -shared
endif
lib: CXXFLAGS += -fPIC
lib: $(LIBOBJ)
	$(MKDIR) $(LIBDIR)
	$(CXX) $(LDFLAGS) -o $(SHAREDLIB) $^

examples/%: examples/%.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $< $(LIB) $(LIBS)

install: $(LIB) $(HEADERS)
	install -d $(DESTDIR)/lib $(DESTDIR)/include/$(PKGNAME)
	install -m644 $(LIB) $(DESTDIR)/lib
ifeq ($(UNAME), Darwin)
	install -m644 $(HEADERS) $(DESTDIR)/include/$(PKGNAME)
else
	install -D -m644 $(HEADERS) $(DESTDIR)/include/$(PKGNAME)
endif

dist:
	@$(MKDIR) $(ARCHIVE)
	@$(CP) LICENSE Makefile README.md $(ARCHIVE)
	@$(MKDIR) $(ARCHIVE)/{examples,src}
	@$(CP) examples/*.cc $(ARCHIVE)/examples
	@$(CP) src/*.cc src/*.h $(ARCHIVE)/src
	@tar -czf $(ARCHIVE).tar.gz $(ARCHIVE)
	@$(RM) -r $(ARCHIVE)
	@echo dist tarball created: $(ARCHIVE).tar.gz

clean::
	$(RM) $(EXE) $(LIBOBJ) $(LIB) $(SHAREDLIB)
	$(RM) -r $(LIBDIR)
