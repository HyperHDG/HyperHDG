#!/bin/sh

####################################################################################################
# The original script is signed by Jeroen de Bruijn. A later version is signed by Francesco Romano.
# This script is a further adaption of the script of Francesco Romano, which can be found under
# https://gist.github.com/francesco-romano/351a6ae457860c14ee7e907f2b0fc1a5 (date: 20. Dec. 2020)
#
# Authors: Andreas Rupp and Guido Kanschat, Heidelberg University, 2020
####################################################################################################

# Setup this script and get the current gh-pages branch.
echo 'Setting up the script...'

# Exit with nonzero exit code if anything fails.
set -e

# Set global variables.
GH_REPO_ORG=AndreasRupp
GH_REPO_NAME=HyperHDG.wiki

# Pretend to be user Andreas Rupp CI.
git config --global user.name "Andreas Rupp CI"
git config --global user.email "HyperHDG@rupp.ink"

# Retrieve master branch of the repositoy containing the wiki.
git clone https://AndreasRuppCI:$REPO_TOKEN@github.com/$GH_REPO_ORG/$GH_REPO_NAME.git code_docs
cd code_docs

# Set the push default to simple i.e. push only the current branch.
git config --global push.default simple

# Go back to first commit.
git reset --hard `git rev-list --max-parents=0 --abbrev-commit HEAD`

# Copy wiki (consisiting of .md files) into current branch and remove auxiliary files.
cp  ../*.md .
rm -f README.md

# Start uploading the wiki pages.
echo 'Uploading documentation to the wiki branch...'

# Add everything in this directory (the Doxygen code documentation) to the gh-pages branch.
git add --all

# Commit the added files with a title and description containing the GitHub actions build number and
# the GitHub commit reference that issued this build.
git commit -m "Deploy Wiki to GitHub gh-wiki build: ${GITHUB_RUN_NUMBER}" \
  -m "Commit: $(git rev-parse --short "$GITHUB_SHA")"

# Force push to the remote GitHub pages branch.
git push --force https://AndreasRuppCI:$REPO_TOKEN@github.com/$GH_REPO_ORG/$GH_REPO_NAME.git
