#!/usr/bin/env bash
# Move rendered Quarto reports (.html) from the analysis folders into docs/,
# keeping the folder structure, and rebuild docs/index.html for GitHub Pages.
#
# Usage (from anywhere):  bash move_reports_to_docs.sh
# Re-run after every render; existing reports in docs/ are replaced by the new ones.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

SOURCE_DIRS=(Organoid_analysis Human_analysis compare_organoid_human)
DOCS_DIR="docs"

mkdir -p "$DOCS_DIR"
touch "$DOCS_DIR/.nojekyll"   # serve files as-is on GitHub Pages (no Jekyll processing)

n_moved=0
for src in "${SOURCE_DIRS[@]}"; do
  [[ -d "$src" ]] || continue
  while IFS= read -r -d '' html; do
    dest="$DOCS_DIR/$html"
    mkdir -p "$(dirname "$dest")"
    mv -f "$html" "$dest"
    # non-self-contained renders put their assets in <name>_files/
    assets="${html%.html}_files"
    if [[ -d "$assets" ]]; then
      rm -rf "${dest%.html}_files"
      mv "$assets" "$(dirname "$dest")/"
    fi
    echo "moved: $html -> $dest"
    n_moved=$((n_moved + 1))
  done < <(find "$src" -path "$src/logs" -prune -o -type f -name '*.html' -print0)
done
echo "$n_moved report(s) moved."

# Rebuild docs/index.html listing every report
{
  echo '<!doctype html><html><head><meta charset="utf-8">'
  echo '<meta name="viewport" content="width=device-width, initial-scale=1">'
  echo '<title>Analysis reports</title>'
  echo '<style>body{font-family:sans-serif;max-width:800px;margin:2em auto;padding:0 1em;line-height:1.6}h2{margin-top:1.5em}</style>'
  echo '</head><body><h1>Analysis reports</h1>'
  echo '<p>Rendered Quarto reports for <a href="https://github.com/">aging_human_organoid_comparisson</a>.</p>'
  for src in "${SOURCE_DIRS[@]}"; do
    [[ -d "$DOCS_DIR/$src" ]] || continue
    echo "<h2>$src</h2><ul>"
    find "$DOCS_DIR/$src" -type f -name '*.html' | sort | while IFS= read -r f; do
      rel="${f#"$DOCS_DIR"/}"
      echo "<li><a href=\"$rel\">$(basename "$f" .html)</a></li>"
    done
    echo "</ul>"
  done
  echo '</body></html>'
} > "$DOCS_DIR/index.html"
echo "updated: $DOCS_DIR/index.html"
