#!/bin/bash
# based on https://stackoverflow.com/a/49516361/7728169

set -e

if [[ "$PY" != "3.7" ]]; then
  # we only deploy with the python 3.7 build
  exit 0
fi

if [[ "$TRAVIS_BRANCH" == "develop" && "$TRAVIS_PULL_REQUEST" == "false" ]]; then
  export DOCS_VERSION="latest"
elif [[ "$TRAVIS_BRANCH" == "master" && "$TRAVIS_PULL_REQUEST" == "false" ]]; then
  export DOCS_VERSION="stable"
elif [ -n "$TRAVIS_TAG" ]; then
  # this is a tag build
  export DOCS_VERSION="$TRAVIS_TAG"
else
  # nothing to deploy
  exit 0
fi

echo "Docs version: $DOCS_VERSION"

PUBLICATION_BRANCH=gh-pages
# Checkout the branch
pwd
ls
ls docs
ls docs/_build
ls docs/_build/html

REPO_PATH=$PWD
pushd $HOME || exit 1
git clone --branch=$PUBLICATION_BRANCH https://${GITHUB_TOKEN}@github.com/$TRAVIS_REPO_SLUG publish
cd publish || exit 1

# Update pages
if [[ -d "$DOCS_VERSION" ]]; then
  git rm -rf "$DOCS_VERSION"
fi
mkdir "$DOCS_VERSION"
cp -r "$REPO_PATH/docs/_build/html/*" "$DOCS_VERSION"
# Commit and push latest version
git add .
git config user.name  "Travis"
git config user.email "travis@travis-ci.org"
git commit -m "Updated distribution."
git push -fq origin $PUBLICATION_BRANCH
popd || exit 1