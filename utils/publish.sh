#!/bin/bash
# Publish a new version of `PID_pendulum' to PyPI. Requires a valid
# username and password.
TWINE=$(pip3 show twine)
ANSI_GREEN='\x1b[0;32m'
ANSI_RESET='\x1b[0m'
DELIM='*--------------------*'

if [ "$TWINE" = '' ]; then
    echo -e "${ANSI_GREEN}Installing dependency 'twine' with pip.\n${DELIM}${ANSI_RESET}"
    sudo pip3 install twine
fi

echo -e "${ANSI_GREEN}Creating distribution.\n${DELIM}${ANSI_RESET}"
make all
echo -e "${ANSI_GREEN}Publishing to PyPI.\n${DELIM}${ANSI_RESET}"
python3 -m twine upload dist/*
