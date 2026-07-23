#!/usr/bin/env bash
# Verify the maintainer staging script fails closed without the RACMO grid mask.

set -euo pipefail

script="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/stage_fast_nc_fixtures.sh"

# Syntax plus all three source-guard lines form the static maintainer contract.
bash -n "$script"
grep -Fq 'if [[ ! -f "$racmo_topo" ]]; then' "$script"
grep -Fq 'Missing RACMO FGRN11 topography/mask source file:' "$script"
grep -Fq 'cp -p "$racmo_topo"' "$script"
