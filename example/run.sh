#!/bin/sh

../../rosetta/main/source/bin/AbinitioRelax.default.linuxgccrelease @flags

rm score.fsc
rm S_00000001.pdb
rm default.out
