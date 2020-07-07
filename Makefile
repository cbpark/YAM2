PKGNAME  := YAM2
SRCDIR 	 := src
LIBDIR 	 := lib
CXXFLAGS := -g -O2 -Wall -Wextra -std=c++14 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LDFLAGS  := -O2
AR       := ar crs
MKDIR    := mkdir -p
RM       := rm -f

LIB    := $(LIBDIR)/lib$(PKGNAME).a
LIBSRC := $(wildcard $(SRCDIR)/*.cc)
LIBOBJ := $(LIBSRC:.cc=.o)
EXE    := analysis/yam2.exe

# NLopt (https://nlopt.readthedocs.io/
NLOPT    ?= /usr
CXXFLAGS += -I$(NLOPT)/include
LDFLAGS  += -L$(NLOPT)/lib -lnlopt

.PHONY: all install clean

all: $(EXE)

$(LIB): CXXFLAGS += -fPIC
$(LIB): $(LIBOBJ)
	$(MKDIR) $(LIBDIR)
	$(AR) $@ $^
	ranlib $@

analysis/%.exe: analysis/%.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $< -L$(LIBDIR) -l$(PKGNAME)

clean::
	$(RM) $(EXE) $(LIBOBJ) $(LIB)
	$(RM) -r $(LIBDIR)
