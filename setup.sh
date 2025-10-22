#!/usr/bin/env bash
# setup.sh - Initial setup for pipeline

set -euo pipefail

echo "==> Setting up Mitochondrial Variant Pipeline"
echo

# Create directories
echo "Creating directory structure..."
mkdir -p config scripts tools logs
echo "  ✓ Directories created"

# Download Picard
PICARD_VERSION="3.1.1"
PICARD_URL="https://github.com/broadinstitute/picard/releases/download/${PICARD_VERSION}/picard.jar"

if [[ ! -f tools/picard.jar ]]; then
  echo
  echo "Downloading Picard ${PICARD_VERSION}..."
  if command -v wget >/dev/null 2>&1; then
    wget -q -O tools/picard.jar "$PICARD_URL"
  elif command -v curl >/dev/null 2>&1; then
    curl -sL -o tools/picard.jar "$PICARD_URL"
  else
    echo "ERROR: Neither wget nor curl found. Please install one."
    exit 1
  fi
  echo "  ✓ Picard downloaded"
else
  echo "  ✓ Picard already exists"
fi

# Make scripts executable
echo
echo "Making scripts executable..."
chmod +x run_mito_pipeline.sh scripts/*.sh 2>/dev/null || true
echo "  ✓ Scripts are executable"

# Check dependencies
echo
echo "Checking dependencies..."
bash check_dependencies.sh || true

echo
echo "=========================================="
echo "Setup complete!"
echo "=========================================="
echo
echo "Next steps:"
echo "  1. Place species data in directories (Species_name/)"
echo "  2. Run: ./run_mito_pipeline.sh --list"
echo "  3. Run: ./run_mito_pipeline.sh"
echo
