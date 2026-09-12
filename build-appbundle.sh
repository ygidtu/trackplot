#!/bin/bash
# trackplot AppBundle creation script
# This script creates a complete AppDir and packages it as AppBundle

set -e

APPDIR_NAME="trackplot.AppDir"
APPDIR="$(pwd)/$APPDIR_NAME"

# Extract version from pyproject.toml
PYPROJECT="$(dirname "$(readlink -f "$0")")/pyproject.toml"
if [ -f "$PYPROJECT" ]; then
    VERSION="$(grep -oP '(?<=^version = ")[^"]+' "$PYPROJECT")"
else
    echo "Error: pyproject.toml not found"
    exit 1
fi

MAINTAINER="ygidtu"
APPBUNDLE_ID="trackplot-${VERSION}-${MAINTAINER}"
OUTPUT_FILE="trackplot-${VERSION}.dwfs.AppBundle"

echo "=== Creating AppDir structure ==="
mkdir -p "$APPDIR/usr/bin"
mkdir -p "$APPDIR/usr/lib/python3/site-packages"
mkdir -p "$APPDIR/usr/share/applications"
mkdir -p "$APPDIR/usr/share/icons/hicolor/256x256/apps"
mkdir -p "$APPDIR/usr/share/trackplot"

echo "=== Installing trackplot to AppDir ==="
# Create a virtual environment in AppDir
python3 -m venv "$APPDIR/.venv"
source "$APPDIR/.venv/bin/activate"

# Install trackplot and dependencies
pip install --upgrade pip
pip install trackplot

# Deactivate venv
deactivate

echo "=== Creating AppRun script ==="
cat > "$APPDIR/AppRun" << 'EOF'
#!/bin/bash
APPDIR="$(readlink -f "${BASH_SOURCE[0]}")"
APPDIR="$(dirname "$APPDIR")"

# Activate the bundled Python environment
source "$APPDIR/.venv/bin/activate"

# Execute trackplot
exec trackplot "$@"
EOF
chmod +x "$APPDIR/AppRun"

echo "=== Creating desktop file ==="
cat > "$APPDIR/usr/share/applications/trackplot.desktop" << EOF
[Desktop Entry]
Name=trackplot
Comment=NGS data visualization toolkit - Track genome sequencing data
Exec=trackplot
Icon=trackplot
Type=Application
Categories=Science;Biology;Bioinformatics;
Terminal=true
EOF

echo "=== Cleaning up venv ==="
# Remove development files to reduce size
rm -rf "$APPDIR/.venv/lib/python*/site-packages/pip"
rm -rf "$APPDIR/.venv/lib/python*/site-packages/setuptools"
rm -rf "$APPDIR/.venv/lib/python*/site-packages/easy_install.py"
rm -rf "$APPDIR/.venv/include"
rm -rf "$APPDIR/.venv/share"

echo "=== AppDir created! ==="
echo "Contents:"
ls -la "$APPDIR/"
echo ""
echo "Size: $(du -sh "$APPDIR" | cut -f1)"

echo ""
echo "=== Next steps ==="
echo "1. Install PELF tool: https://github.com/nakedib/pelf"
echo "2. Package as AppBundle:"
echo "   ./pelf --add-appdir \"$APPDIR\" --appbundle-id \"$APPBUNDLE_ID\" --output-to \"$OUTPUT_FILE\""
echo ""
echo "Or use pelfCreator for automated build:"
echo "   pelfCreator --maintainer \"$MAINTAINER\" --name \"org.trackplot.Trackplot\" --entrypoint \"trackplot.desktop\""
echo ""
echo "=== Output ==="
echo "  Version:     $VERSION"
echo "  AppBundle:   $OUTPUT_FILE"
echo "  AppBundleID: $APPBUNDLE_ID"
