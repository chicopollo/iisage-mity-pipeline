#!/usr/bin/env bash
# check_dependencies.sh - Verify all dependencies

set -euo pipefail

echo "==> Checking Pipeline Dependencies"
echo

ALL_GOOD=1

check_cmd() {
  local cmd=$1
  local msg=${2:-}
  
  if command -v "$cmd" >/dev/null 2>&1; then
    echo "✅ $cmd: $(command -v "$cmd")"
    return 0
  else
    echo "❌ $cmd: NOT FOUND ${msg}"
    ALL_GOOD=0
    return 1
  fi
}

check_docker() {
  if command -v docker >/dev/null 2>&1; then
    if docker info >/dev/null 2>&1; then
      echo "✅ Docker: Running"
      return 0
    elif sudo docker info >/dev/null 2>&1; then
      echo "⚠️  Docker: Running (requires sudo)"
      return 0
    else
      echo "❌ Docker: Cannot connect to daemon"
      ALL_GOOD=0
      return 1
    fi
  else
    echo "❌ Docker: NOT FOUND"
    ALL_GOOD=0
    return 1
  fi
}

echo "Core tools:"
check_cmd bwa "(install: apt-get install bwa)"
check_cmd samtools "(install: apt-get install samtools)"
check_cmd bcftools "(install: apt-get install bcftools)"
check_cmd bgzip "(install: apt-get install tabix)"
check_cmd tabix "(install: apt-get install tabix)"
check_cmd seqtk "(install: apt-get install seqtk)"
check_cmd java "(install: apt-get install openjdk-11-jre)"
check_cmd ssconvert "(install: apt-get install gnumeric)"

echo
echo "Docker:"
check_docker

echo
echo "Pipeline files:"
if [[ -f tools/picard.jar ]]; then
  echo "✅ Picard: tools/picard.jar"
else
  echo "❌ Picard: NOT FOUND (run: bash setup.sh)"
  ALL_GOOD=0
fi

if [[ -f config/pipeline.conf ]]; then
  echo "✅ Config: config/pipeline.conf"
else
  echo "❌ Config: NOT FOUND"
  ALL_GOOD=0
fi

echo
echo "=========================================="
if [[ $ALL_GOOD -eq 1 ]]; then
  echo "✅ All dependencies satisfied!"
  echo "=========================================="
  exit 0
else
  echo "❌ Some dependencies missing"
  echo "=========================================="
  echo
  echo "Installation guide:"
  echo "  Ubuntu/Debian:"
  echo "    sudo apt-get install bwa samtools bcftools tabix seqtk gnumeric docker.io openjdk-11-jre"
  echo
  echo "  After installing Docker:"
  echo "    sudo usermod -aG docker \$USER"
  echo "    # Log out and back in"
  echo
  exit 1
fi
