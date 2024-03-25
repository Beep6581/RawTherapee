#!/bin/sh

if [ -d CVS ] ; then
  echo 'Use mkdist script in cvs export-ed dirs'
else
  if [ -f Makefile.devel ] ; then
    autoreconf --install
    mkdir dcraw
    cd dcraw
    wget https://www.dechifro.org/dcraw/dcraw.c
    cd ..
    rm -f clist2c.pl clist2html.pl
    rm -f Makefile.devel 
    rm mkdist.sh export-dist.sh
  else
   echo 'Wrong directory or mkdist.sh already executed'
  fi
fi
