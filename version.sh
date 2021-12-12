#!/bin/sh -e

if [ -e .git ]; then
    version=$(git describe --tags | cut -d - -f 1-2 | tr - .)
elif [ -n "$PORTVERSION" ]; then
    version=$PORTVERSION
elif [ -n "$PKGVERSION" ]; then
    version=$PKGVERSION
else
    version="Unknown-version"
fi
echo $version

