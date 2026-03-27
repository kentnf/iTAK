#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")/.."

python_files=()
while IFS= read -r file; do
  python_files+=("$file")
done < <(rg --files -g '*.py' -g '!database/**' -g '!bin/**' -g '!__pycache__/**')

python -m py_compile "${python_files[@]}"
python -m unittest discover -s tests
