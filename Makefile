PKGNAME  := YAM2
VERSION  := 0.2.0.1
SRCDIR 	 := src
LIBDIR 	 := lib
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LIBS     := -lm
AR       := ar crs
MKDIR    := mkdir -p
RM       := rm -f
UNAME    := $(shell uname -s)

LIB    := $(LIBDIR)/lib$(PKGNAME).a
LIBSRC := $(wildcard $(SRCDIR)/*.cc)
LIBOBJ := $(LIBSRC:.cc=.o)
EXE    := examples/m2 \
	  examples/m2cons
ifeq ($(UNAME), Darwin)
SHAREDLIB  := $(LIBDIR)/lib$(PKGNAME).dylib
else
SHAREDLIB  := $(LIBDIR)/lib$(PKGNAME).so
endif

DESTDIR ?= /usr/local
HEADERS := $(wildcard $(SRCDIR)/*.h)

# NLopt (https://nlopt.readthedocs.io/
NLOPT    ?= /usr
CXXFLAGS += -I$(NLOPT)/include
LIBS     += -L$(NLOPT)/lib -lnlopt -Wl,-rpath $(NLOPT)/lib

.PHONY: all lib install clean

all: $(LIB)

$(LIB): CXXFLAGS += -fPIC
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

install: $(LIB) lib $(HEADERS)
	install -d $(DESTDIR)/lib $(DESTDIR)/include/$(PKGNAME)
	install -m755 $(SHAREDLIB) $(DESTDIR)/lib
	install -m644 $(LIB) $(DESTDIR)/lib
ifeq ($(UNAME), Darwin)
	install -m644 $(HEADERS) $(DESTDIR)/include/$(PKGNAME)
else
	install -D -m644 $(HEADERS) $(DESTDIR)/include/$(PKGNAME)
endif

clean::
	$(RM) $(EXE) $(LIBOBJ) $(LIB) $(SHAREDLIB)
	$(RM) -r $(LIBDIR)
