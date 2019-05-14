#!/bin/python
#
# Author: Daniel Schuette
# License: GPL-3.0
# Date: 12/05/2019
import numpy as np
from flask import Flask, jsonify, request

from PID_control import Pendulum

app = Flask(__name__, static_url_path="/static")


# this is the API route that retrieves the pendulum data, required URL
# parameters are described below
@app.route("/api/v1/", methods=["GET"])
def api():
    """
    The following parameters are extracted from the URL:
        - alpha, beta, mu
        - phi0, phi0_dot
        - max_control, deadband, set_point
        - t_start=0.0, t_end=30.0
        - N=(t_end-t_start)*100
        - nonlinear=True
        - L=0.1
        - precision=5, frequency=10
        - API key
    """
    required = [
        "alpha", "beta", "mu", "phi0", "phi0_dot", "max_control",
        "deadband", "set_point", "key"
    ]  # list of required params
    optional = {
        "t_start": 0, "t_end": 30, "N": 3000, "nonlinear": True, "L": 0.1,
        "frequency": 10, "precision": 5
    }  # dict of optional params and default values

    # get params from URL and save them to a dict
    # only params with defaults don't result in an error if not provided
    params = {}

    # TODO: check API key before performing any actions
    # TODO: validate user request (some values must be within certain limits)
    # TODO: use different methods (`solve' and `solve_from_angle' must be
    #       parameterized)

    for param in required:
        if param in request.args:
            params[param] = request.args[param]
        else:
            return jsonify("Error: Parameter {} is required.".format(param))

    for param, default in optional.items():
        if param in request.args:
            params[param] = request.args[param]
        else:
            params[param] = default

    # instantiate PIDControl object and return its results
    pendulum = Pendulum(
        float(params["t_start"]), float(params["t_end"]),
        int(params["N"]), np.sin, 0.1
    )
    pendulum.solve(
        float(params["phi0"]), float(params["phi0_dot"]),
        float(params["alpha"]), float(params["beta"]), float(params["mu"]),
        float(params["max_control"]), int(params["frequency"]),
        float(params["deadband"]), float(params["set_point"]),
        int(params["precision"])
    )

    # collect angles and support values into a dictionary and return the JSON
    angles = {
        "angles": pendulum.get_func_values(),
        "support_values": pendulum.get_support_values()
    }
    return jsonify(angles)


# the following routes serve static web pages for the actual web app
@app.route("/")
def home():
    return app.send_static_file("index.html")


@app.route("/about")
def about():
    return app.send_static_file("about.html")


@app.route("/bibliography")
def bibliography():
    return app.send_static_file("bibliography.html")


@app.route("/license")
def license():
    return app.send_static_file("license.html")


@app.route("/contribute")
def contribute():
    return app.send_static_file("contribute.html")


if __name__ == "__main__":
    app.run(debug=True)
