#!/usr/bin/env bash
set -euo pipefail

if [[ $# -gt 2 ]]; then
  echo "Usage: $0 [db-tag] [output-dir]" >&2
  exit 1
fi

db_tag="${1:-db-v1}"
output_dir="${2:-dist}"

repo_root="$(cd "$(dirname "$0")/.." && pwd)"
database_dir="${repo_root}/database"
root_dir_name="iTAK-${db_tag}"
asset_name="itak-${db_tag}.tar.gz"
sha256_name="${asset_name}.sha256"

required_files=(
  "TF_Rule.txt"
  "GA_table.txt"
  "PK_class_desc.txt"
  "Tfam_domain.hmm"
  "TF_selfbuild.hmm"
  "PlantsPHMM3_89.hmm"
  "Plant_Pkinase_fam.hmm"
  "Pkinase_sub_WNK1.hmm"
  "Pkinase_sub_MAK.hmm"
)

for file_name in "${required_files[@]}"; do
  if [[ ! -f "${database_dir}/${file_name}" ]]; then
    echo "Missing required database file: ${database_dir}/${file_name}" >&2
    exit 1
  fi
done

mkdir -p "${repo_root}/${output_dir}"

work_dir="$(mktemp -d)"
trap 'rm -rf "${work_dir}"' EXIT

stage_root="${work_dir}/${root_dir_name}"
mkdir -p "${stage_root}"
cp -R "${database_dir}" "${stage_root}/database"

archive_path="${repo_root}/${output_dir}/${asset_name}"
sha256_path="${repo_root}/${output_dir}/${sha256_name}"

rm -f "${archive_path}" "${sha256_path}"

tar -czf "${archive_path}" -C "${work_dir}" "${root_dir_name}"
shasum -a 256 "${archive_path}" | awk '{print $1 "  '"${asset_name}"'"}' > "${sha256_path}"

echo "Created release assets:"
echo "  ${archive_path}"
echo "  ${sha256_path}"
echo
echo "Suggested next steps:"
echo "  1. Create or update GitHub release tag: ${db_tag}"
echo "  2. Upload both files as release assets"
echo "  3. Verify: python -m itak db download --url <asset-url> --sha256-url <sha256-url>"
