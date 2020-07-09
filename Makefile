PKGNAME  := YAM2
SRCDIR 	 := src
LIBDIR 	 := lib
CXXFLAGS := -g -O0 -Wall -Wextra -std=c++17 -pedantic -I$(SRCDIR) $(CXXFLAGS)
LDFLAGS  := -O0
LIBS     := -lm
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
LIBS     += -L$(NLOPT)/lib -lnlopt

.PHONY: all install clean

all: $(EXE)

$(LIB): CXXFLAGS += -fPIC
$(LIB): $(LIBOBJ)
	$(MKDIR) $(LIBDIR)
	$(AR) $@ $^
	ranlib $@

analysis/%.exe: analysis/%.o $(LIB)
	$(CXX) $(LDFLAGS) -o $@ $< -L$(LIBDIR) -l$(PKGNAME) $(LIBS)

clean::
	$(RM) $(EXE) $(LIBOBJ) $(LIB)
	$(RM) -r $(LIBDIR)
