#! /bin/bash
export ASDAR="https://github.com/rsbivand/asdar-book/raw/master/dump_from_bitbucket_repo/"
#http://www.asdar-book.org/"
export ASDAR_DOWNLOAD="${ASDAR}datasets"
export ASDAR_BUNDLES="${ASDAR}bundles"
export R_LIBS="/home/rsb/lib/r_libs"
export R_CMD="/home/rsb/topics/R/R351-share/bin/R"
export RUN_IN="/home/rsb/proj/other/book/tests"

cd $RUN_IN

CHAPTERS="hello cm vis die cm2 csdacm sppa geos lat1 lat2 dismap"

for ch in $CHAPTERS; do
  if test -d ${ch}
    then cd ${ch}
  else
    mkdir ${ch}
    cd ${ch}
  fi
  wget -N "${ASDAR_BUNDLES}/${ch}_bundle.zip"
  unzip -o "${ch}_bundle.zip"
  echo "source(\"${ch}_mod.R\", echo=TRUE)" > scrpt
  $R_CMD --vanilla < scrpt 2>&1 > ${ch}.log
  if test $? -ne 0
  then 
    echo "ASDAR $ch test failure" | mail -s "$ch test failure" Roger.Bivand@nhh.no
    echo "ASDAR $ch test failure"
  fi
  cd ..
done

