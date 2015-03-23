#! /bin/bash
DIR=$( cd "$( dirname "$0" )" && pwd )
mkdir -p $DIR/../docs/html
rst2html $DIR/../docs/src/readme.txt $DIR/../docs/html/readme.html
rst2html $DIR/../docs/src/tutorial.txt $DIR/../docs/html/tutorial.html
