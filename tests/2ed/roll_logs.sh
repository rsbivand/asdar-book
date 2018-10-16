#! /bin/bash
export RUN_IN="/home/rsb/proj/other/book/tests2ed"

cd $RUN_IN

CHAPTERS="hello cm vis die cm2 std sppa geos lat dismap"

for ch in $CHAPTERS; do
  cd ${ch}
  cp ${ch}.log ${ch}.log.save
  echo "rolled ${ch}.log"
  cd ..
done
