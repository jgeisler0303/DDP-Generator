#!/bin/bash
CUR_DIR=$(pwd)
MY_PATH="`dirname \"$0\"`"

cd $MY_PATH
maxima --batch-string="problem_file: \"$CUR_DIR/$1\"; batchload(\"make_iLQG_MMex.mac\");"

cd $CUR_DIR