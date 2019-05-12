#!/bin/sh
# Make a test request to the PID pendulum controller API.
echo 'Response:'
if [ "$1" = 'err' ]; then
    curl 'http://localhost:5000/api/v1/?invalid_param'
else
    curl 'http://localhost:5000/api/v1/?alpha=0&beta=0&mu=0&phi0=0&phi0_dot=0&max_control=0&frequency=0&deadband=0&set_point=0&precision=0&key=0&N=10'
fi
echo
