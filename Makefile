PKGNAME  := YAM2
VERSION  := 0.0.0.0
SRCDIR 	 := src
LIBDIR 	 := lib
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LDFLAGS  := -O2
LIBS     := -lm
AR       := ar crs
MKDIR    := mkdir -p
RM       := rm -f

LIB    := $(LIBDIR)/lib$(PKGNAME).a
LIBSRC := $(wildcard $(SRCDIR)/*.cc)
LIBOBJ := $(LIBSRC:.cc=.o)
EXE    := examples/m2

# NLopt (https://nlopt.readthedocs.io/
NLOPT    ?= /usr
CXXFLAGS += -I$(NLOPT)/include
LIBS     += -L$(NLOPT)/lib -lnlopt -Wl,-rpath $(NLOPT)/lib

.PHONY: all install clean

all: $(LIB)

$(LIB): CXXFLAGS += -fPIC
$(LIB): $(LIBOBJ)
	$(MKDIR) $(LIBDIR)
	$(AR) $@ $^
	ranlib $@

examples/%: examples/%.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $< -L$(LIBDIR) -l$(PKGNAME) $(LIBS)

clean::
	$(RM) $(EXE) $(LIBOBJ) $(LIB)
	$(RM) -r $(LIBDIR)
