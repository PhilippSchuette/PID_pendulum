#!/bin/bash
# Deploy to Heroku.
JS_SRC='static/scripts/main.js'
REQ_TXT='requirements.txt'
URL_DEPLOY='"https://pid-pendulum-demo.herokuapp.com"'
URL_LOCAL='"http://localhost:5000"'
ANSI_GREEN='\x1b[0;32m'
ANSI_RED='\x1b[0;31m'
ANSI_RESET='\x1b[0m'
DELIM='*--------------------*'

if [[ -z "$1" ]]; then
    echo -e "${ANSI_RED}Error: Provide a version of PID_pendulum as a parameter.\n${DELIM}${ANSI_RESET}"
    exit 1
elif [[ "$1" != v*.*.* ]];then
    echo -e "${ANSI_RED}Error: Version number must be in format 'vX.Y.Z'.\n${DELIM}${ANSI_RESET}"
    exit 1
fi
PID_VERSION="$1"

if [ -z "$2" ]; then
    echo -e "${ANSI_RED}Error: Provide the path to the app as a parameter.\n${DELIM}${ANSI_RESET}"
    exit 1
fi
APP_PATH="$2"

# change in-code constants to deployment settings
echo -e "${ANSI_GREEN}Setting up for deployment.\n${DELIM}${ANSI_RESET}"
cd "$APP_PATH" || exit 1
sed -i "s/const baseURL = $URL_LOCAL/const baseURL = $URL_DEPLOY/" "$JS_SRC"
sed -i "s/PID-pendulum==[a-zA-Z0-9.'\" ]*/PID-pendulum==$PID_VERSION/" "$REQ_TXT"

# revert everything back to local development settings (except PID_pendulum version)
sed -i "s/const baseURL = $URL_DEPLOY/const baseURL = $URL_LOCAL/" "$JS_SRC"
