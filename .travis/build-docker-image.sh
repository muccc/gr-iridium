#!/bin/bash
set -ev

# Increment version number to flush cache if required
TAG="${1}"
echo "Using tag ${TAG}"

# Create cache folder when running outside Travis
if [ -d "docker" ]; then
    echo "Found existing docker cache folder"
else
    mkdir docker
fi

# Load any image from cache if available
CACHE_FROM=""
if [ -f "docker/${TAG}.tar.gz" ]; then
    gzip -dc "docker/${TAG}.tar.gz" | docker load
    CACHE_FROM="${TAG}"
fi

# build the image
docker build --cache-from="${CACHE_FROM}" --tag="${TAG}" --file=".travis/builder.dockerfile" .travis

# Save Docker image to cache
if [ -f "docker/${TAG}.tar.gz" ]; then
    echo "Leaving cache unchanged"
else
    docker save "${TAG}" | gzip > docker/${TAG}.tar.gz
fi

# Clean cache of any previous versions
ls -ahls docker
find docker/ -type f ! -name "${TAG}.tar.gz" -delete
ls -ahls docker