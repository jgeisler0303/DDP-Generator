#!/bin/bash
CUR_DIR=$(pwd)
MY_PATH="`dirname \"$0\"`"

cd $MY_PATH
maxima -q --batch-string="problem_file: \"$CUR_DIR/$1\"$ batchload(\"make_iLQG.mac\")$"

cd $CUR_DIR