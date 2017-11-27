#!/usr/bin/env bash

if [ "$1" == "major" ] || [ "$1" == "minor" ] || [ "$1" == "patch" ] || [ "$1" == "p" ] ;
then

sphinx-build -b html docs/sphinx_files/ docs/ -a
git add .
git commit -m "deploy $1"

if [ "$1" ==  "major" ]
then
bumpversion major
elif [ "$1" ==  "minor" ]
then
bumpversion minor
elif [ "$1" == "patch" ] || [ "$1" == "p" ];
then
bumpversion patch
fi

git push origin master && git push origin master --tags
pip install tribology --upgrade --no-cache-dir
echo "deployed successfully as $1"

else
echo "parameter '$1' not defined"
fi
