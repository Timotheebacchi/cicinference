#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
release_dir="${script_dir}/release"
package_dir="${script_dir}/package"

required_files=(
  "${script_dir}/ado/cicinference.ado"
  "${script_dir}/ado/cicinference.sthlp"
  "${script_dir}/mata/cicinference_quantiles.mata"
  "${script_dir}/mata/cicinference_density.mata"
  "${script_dir}/mata/cicinference_variance.mata"
  "${script_dir}/mata/cicinference_main.mata"
  "${script_dir}/mata/cicinference_bootstrap.mata"
  "${package_dir}/cicinference.pkg"
  "${package_dir}/stata.toc"
)

for file in "${required_files[@]}"; do
  if [[ ! -f "${file}" ]]; then
    echo "Missing required runtime file: ${file}" >&2
    exit 1
  fi
done

rm -rf "${release_dir}"
mkdir -p "${release_dir}"

cp "${script_dir}/ado/cicinference.ado" "${release_dir}/"
cp "${script_dir}/ado/cicinference.sthlp" "${release_dir}/"
cp "${script_dir}/mata/cicinference_quantiles.mata" "${release_dir}/"
cp "${script_dir}/mata/cicinference_density.mata" "${release_dir}/"
cp "${script_dir}/mata/cicinference_variance.mata" "${release_dir}/"
cp "${script_dir}/mata/cicinference_main.mata" "${release_dir}/"
cp "${script_dir}/mata/cicinference_bootstrap.mata" "${release_dir}/"
cp "${package_dir}/cicinference.pkg" "${release_dir}/"
cp "${package_dir}/stata.toc" "${release_dir}/"

(
  cd "${release_dir}"
  zip -q cicinference-stata-release.zip \
    cicinference.ado \
    cicinference.sthlp \
    cicinference_quantiles.mata \
    cicinference_density.mata \
    cicinference_variance.mata \
    cicinference_main.mata \
    cicinference_bootstrap.mata \
    cicinference.pkg \
    stata.toc
)

echo "Distributed files:"
find "${release_dir}" -maxdepth 1 -type f -print | sort
