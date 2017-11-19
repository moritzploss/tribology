sphinx-build -b html docs/sphinx_files/ docs/ -a
git add .
git commit -m "deploy $1"

if [ "$2" == "major" ] || [ "$2" == "minor" ] || [ "$2" == "patch" ] || [ "$2" == "p" ] ;
then

if [ "$2" ==  "major" ]
then
    echo $2
    #bumpversion major
elif [ "$2" ==  "minor" ]
then
    echo $2
    #bumpversion minor
elif [ "$2" == "patch" ] || [ "$2" == "p" ];
then
    echo $2
    #bumpversion patch
fi
git push origin master && git push origin master --tags
echo "deployed successfully"
else
echo "parameter '$2' not defined"
fi
