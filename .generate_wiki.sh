#!/bin/sh

####################################################################################################
# The original script is signed by Jeroen de Bruijn. A later version is signed by Francesco Romano.
# This script is a further adaption of the script of Francesco Romano, which can be found under
# https://gist.github.com/francesco-romano/351a6ae457860c14ee7e907f2b0fc1a5 (date: 20. Dec. 2020)
#
# Authors: Andreas Rupp and Guido Kanschat, Heidelberg University, 2020
###################################################################################################

# Setup this script and get the current gh-pages branch.
echo 'Setting up the script...'

# Exit with nonzero exit code if anything fails.
set -e

# Set global variables.
# GH_REPO_ORG=`echo $TRAVIS_REPO_SLUG | cut -d "/" -f 1`
# GH_REPO_NAME=`echo $TRAVIS_REPO_SLUG | cut -d "/" -f 2`
GH_REPO_ORG=AndreasRupp
GH_REPO_NAME=HyperHDG.wiki

# Retrieve master branch of the repositoy containing the wiki.
git clone https://AndreasRuppTravis:$REPO_TOKEN@github.com/$GH_REPO_ORG/$GH_REPO_NAME.git code_docs
cd code_docs

# Set the push default to simple i.e. push only the current branch.
git config --global push.default simple

# Go back to first commit.
git reset --hard `git rev-list --max-parents=0 --abbrev-commit HEAD`

# Copy wiki into current branch and remove auxiliary files.
cp  ../* .
rm -f .travis.yml .generate_wiki.sh .gitignore

# Start uploading the wiki pages.

echo 'Uploading documentation to the wiki branch...'

# Add everything in this directory (the Doxygen code documentation) to the gh-pages branch.
git add --all

# Commit the added files with a title and description containing the Travis CI build number and
# the GitHub commit reference that issued this build.
git commit -m "Deploy code docs to GitHub Pages Travis build: ${TRAVIS_BUILD_NUMBER}" \
  -m "Commit: ${TRAVIS_COMMIT}"

# Force push to the remote GitHub pages branch.
git push --force https://AndreasRuppTravis:$REPO_TOKEN@github.com/$GH_REPO_ORG/$GH_REPO_NAME.git

