#!/bin/bash
# Publish a new version of `PID_pendulum' to PyPI. Requires a valid
# username and password as well as a new version as an input parameter.
DOC_CONFIG='PID_pendulum/conf.py'
PKG_CONFIG='setup.py'
TWINE=$(pip3 show twine)
ANSI_GREEN='\x1b[0;32m'
ANSI_RED='\x1b[0;31m'
ANSI_RESET='\x1b[0m'
DELIM='*--------------------*'

if [[ -z "$1" ]]; then
    echo -e "${ANSI_RED}Error: Provide a new version number as a parameter.${ANSI_RESET}"
    exit 1
elif [[ "$1" != v*.*.* ]];then
    echo -e "${ANSI_RED}Error: Version number must be in format 'vX.Y.Z'.${ANSI_RESET}"
    exit 1
fi
VERSION="$1"

if [[ "$TWINE" = '' ]]; then
    echo -e "${ANSI_GREEN}Installing dependency 'twine' with pip.\n${DELIM}${ANSI_RESET}"
    sudo pip3 install twine
fi

# make the python module after replace old with new version information
echo -e "${ANSI_GREEN}Creating distribution.\n${DELIM}${ANSI_RESET}"
sed -i "s/version = \"[a-zA-Z0-9.'\" ]*/version = \"$VERSION\"/" "$PKG_CONFIG"
make clean
make all

# create documentation
echo -e "${ANSI_GREEN}Building documentation.\n${DELIM}${ANSI_RESET}"
# replace exact version (if version > 0.x.y, major version needs to
# be modified separately though)
sed -i "s/release = u'[a-zA-Z0-9.'\" ]*/release = u'$VERSION'/" "$DOC_CONFIG"
make docs

# publish the module to PyPI
echo -e "${ANSI_GREEN}Publishing to PyPI.\n${DELIM}${ANSI_RESET}"
python3 -m twine check dist/*
python3 -m twine upload dist/*

# clean up build artifacts
echo -e "${ANSI_GREEN}Cleaning up.\n${DELIM}${ANSI_RESET}"
make clean
