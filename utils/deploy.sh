#!/bin/bash
# Deploy to Heroku. The commit message is generic, though.
JS_SRC='static/scripts/main.js' # relative to app root
REQ_TXT='requirements.txt' # relative to app root
ANSI_GREEN='\x1b[0;32m'
ANSI_RED='\x1b[0;31m'
ANSI_RESET='\x1b[0m'
DELIM='*--------------------*'

if [[ -z "$1" ]]; then
    echo -e "${ANSI_RED}Error: Provide a version of PID_pendulum as a first parameter.${ANSI_RESET}"
    exit 1
elif [[ "$1" != v*.*.* ]];then
    echo -e "${ANSI_RED}Error: Version number must be in format 'vX.Y.Z'.${ANSI_RESET}"
    exit 1
fi
PID_VERSION="$1"

if [ -z "$2" ]; then
    echo -e "${ANSI_RED}Error: Provide the path to the app as a second parameter.${ANSI_RESET}"
    exit 1
fi
APP_PATH="$2"

# change in-code constants to deployment settings
echo -e "${ANSI_GREEN}Setting up for deployment.\n${DELIM}${ANSI_RESET}"
cd "$APP_PATH" || exit 1
sed -i "s%http://localhost:5000%https://pid-pendulum-demo\.herokuapp\.com%" "$JS_SRC"
sed -i "s/PID-pendulum==[a-zA-Z0-9.'\" ]*/PID-pendulum==$PID_VERSION/" "$REQ_TXT"

# deploy
echo -e "${ANSI_GREEN}Deploying app to Heroku.\n${DELIM}${ANSI_RESET}"
git add --all
git commit -m "update to PID-pendulum $PID_VERSION"
git push heroku master

# revert everything back to local development settings (except PID_pendulum version)
echo -e "${ANSI_GREEN}Cleaning up.\n${DELIM}${ANSI_RESET}"
sed -i "s%https://pid-pendulum-demo\.herokuapp\.com%http://localhost:5000%" "$JS_SRC"
