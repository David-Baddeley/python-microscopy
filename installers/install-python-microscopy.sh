#!/usr/bin/env bash
# Linux/Mac installer for PYME. Also serves as the Mac installer pending a native .app.
# Usage: install-python-microscopy.sh [--dest <path>]   (default: ~/PYME)
#
# Context A (CI/dev): uv is expected to be in PATH already.
# Context B (end-user): uv is bootstrapped via astral.sh if not found.
set -euo pipefail

# --- Configuration ---
TARGET_PYTHON=3.13
PACKAGE_NAME=python-microscopy
ENTRY_POINTS=(PYMEAcquire PYMEImage PYMEVis PYMEClusterOfOne)
DEFAULT_DEST="$HOME/PYME"

DEST="$DEFAULT_DEST"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dest)   DEST="$2"; shift 2 ;;
        --dest=*) DEST="${1#--dest=}"; shift ;;
        -h|--help)
            echo "Usage: $0 [--dest <path>]"
            echo "  Installs PYME into a self-contained directory using uv."
            echo "  Default: $DEFAULT_DEST"
            exit 0 ;;
        *) echo "Unknown argument: $1 (run '$0 --help' for usage)" >&2; exit 1 ;;
    esac
done

DEST="${DEST/#\~/$HOME}"
echo "==> Installing PYME to: $DEST"
mkdir -p "$DEST"

# --- Ensure uv is available ---
if ! command -v uv &>/dev/null; then
    echo "==> uv not found; downloading via astral.sh..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    # Add the two common install locations; one will contain the binary.
    export PATH="$HOME/.local/bin:${CARGO_HOME:-$HOME/.cargo}/bin:$PATH"
    if ! command -v uv &>/dev/null; then
        echo "ERROR: uv installation failed or installed to an unexpected location." >&2
        exit 1
    fi
fi

# --- Standalone Python installation via py-app-standalone ---
echo "==> Creating standalone Python installation..."
uv tool run py-app-standalone --python-version "$TARGET_PYTHON" --target "$DEST" --force "$PACKAGE_NAME"

# --- Normalize cpython-* to a stable directory name ---
PYTHON_DIRS=( "$DEST"/cpython-* )
[[ -d "${PYTHON_DIRS[0]}" ]] || { echo "ERROR: standalone Python directory not found" >&2; exit 1; }
mv "${PYTHON_DIRS[0]}" "$DEST/python"
# Remove the bare-venv temp directory that py-app-standalone leaves behind.
rm -rf "$DEST/bare-venv"

# --- Top-level entry point symlinks ---
for ep in "${ENTRY_POINTS[@]}"; do
    ln -sf "$DEST/python/bin/$ep" "$DEST/$ep"
done

# --- Shell helper (adds python/bin to PATH) ---
cat > "$DEST/pyme-shell" <<'SHELL_SCRIPT'
#!/usr/bin/env bash
PYME_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PATH="$PYME_DIR/python/bin:$PATH"
exec "${SHELL:-bash}"
SHELL_SCRIPT
chmod +x "$DEST/pyme-shell"

echo ""
echo "==> Done."
echo "    Entry points: ${ENTRY_POINTS[*]}"
echo "    Add to PATH:       export PATH=\"$DEST:\$PATH\""
echo "    Activated shell:   $DEST/pyme-shell"
