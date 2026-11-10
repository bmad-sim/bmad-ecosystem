# Build the standalone test program for the GTK plot window (../qp_gtk_window.c)
# with the Bmad gcc toolchain.
#   make                         (uses gcc from PATH)
#   make CC=/opt/local/bin/gcc   (MacPorts gcc)
# Needs GTK 3 (MacPorts: gtk3 +quartz), Cairo, and readline.

CC ?= gcc
PKG_CONFIG ?= pkg-config

CFLAGS  = -O2 -Wall -std=gnu99 -DCESR_GTK -I.. $(shell $(PKG_CONFIG) --cflags gtk+-3.0)
LDLIBS  = $(shell $(PKG_CONFIG) --libs gtk+-3.0) -lreadline -lm

qp_gtk_prototype: qp_gtk_prototype.c ../qp_gtk_window.c ../qp_gtk_window.h
	$(CC) $(CFLAGS) -o $@ qp_gtk_prototype.c ../qp_gtk_window.c $(LDFLAGS) $(LDLIBS)

clean:
	rm -f qp_gtk_prototype

.PHONY: clean
