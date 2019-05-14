#!/bin/sh
# Make a test request to the PID pendulum controller API.
echo 'Response:'
if [ "$1" = 'err' ]; then
    curl 'http://localhost:5000/api/v1/?invalid_param'
else
    curl 'http://localhost:5000/api/v1/?alpha=4.0&beta=1.5&mu=0.8&phi0=0.25&phi0_dot=0.25&max_control=3.0&frequency=10&deadband=0.01&set_point=0&precision=5&key=0&N=1000'
fi
echo
