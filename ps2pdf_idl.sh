#!/bin/sh
#
# Wrapper of PS2PDF. Before directly calling ps2pdf, ps file
# will be modify so that pdf results will be proper direction
# to view.
# Warning:
#   this script will modify orignal postscritp file.
#
# Usage:
# > ps2pdf infile.ps outfile.pdf
#   infile.ps   : Postscript input file name
#   outfile.pdf : PDF output file name
#
# hiro@Caltech
# 05/09/10 alpha version

tmpfile1=".ps2pdf_idl.ps"

bb_rng=`grep "%%BoundingBox" $1`
bb_rng1=`echo $bb_rng | awk '{print $2}'`
bb_rng2=`echo $bb_rng | awk '{print $3}'`
bb_rng3=`echo $bb_rng | awk '{print $4}'`
bb_rng4=`echo $bb_rng | awk '{print $5}'`

sed -e "s/270 rotate/90 rotate/" $1 > $tmpfile1
sed -e "s/$bb_rng1 $bb_rng4 translate/$bb_rng3 $bb_rng2 translate/" $tmpfile1 > $1

ps2pdf $1 $2

exit
