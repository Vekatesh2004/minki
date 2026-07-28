#!/usr/bin/env bash
set -euo pipefail

REPO="venky04/minki"
TAG="${1:-latest}"

if ! command -v docker >/dev/null 2>&1; then
  echo "Docker is not installed or not on PATH." >&2
  exit 1
fi

if ! docker info >/dev/null 2>&1; then
  echo "Docker daemon is not running. Start Docker Desktop or docker service first." >&2
  exit 1
fi

echo "Building backend image..."
docker build -f docker/Dockerfile.backend -t "$REPO:backend-$TAG" .

echo "Building frontend image..."
docker build -f docker/Dockerfile.frontend -t "$REPO:frontend-$TAG" ./frontend

echo "Logging in to Docker Hub..."
docker login

echo "Pushing backend image..."
docker push "$REPO:backend-$TAG"

echo "Pushing frontend image..."
docker push "$REPO:frontend-$TAG"

echo "Done. Images published to $REPO"
