#!/bin/python
#
# Author: Daniel Schuette
# License: GPL-3.0
# Date: 12/05/2019
import numpy as np
from flask import Flask, jsonify, render_template, request
from PID_pendulum.PID_control import Pendulum

app = Flask(__name__, static_url_path="/static")


# this is the API route that retrieves the pendulum data, required URL
# parameters are described below
@app.route("/api/v1/", methods=["GET"])
def api():
    """
    An API wrapper around the `PID_pendulum' module. Optional and required
    parameters can be seen below. Refer to in-line documentation for more
    implementation details. Unreasonable values do *not* result in errors
    but unreasonable output.
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
    # TODO: check API key before performing any actions
    # TODO: use different methods (`solve' and `solve_from_angle' must be
    #       parameterized)
    angles = {
        "angles": pendulum.get_func_values(),
        "support_values": pendulum.get_support_values(),
        "xy_coordinates": pendulum.get_xy_coordinates()
    }
    return jsonify(angles)


# the following routes serve static web pages for the actual web app
@app.route("/")
@app.route("/home")
def home():
    return render_template("index.html", navbar="home")


@app.route("/about")
def about():
    return render_template("about.html", navbar="about")


@app.route("/bibliography")
def bibliography():
    return render_template("bibliography.html", navbar="bibliography")


@app.route("/license")
def license():
    return render_template("license.html", navbar="dropdown")


@app.route("/contribute")
def contribute():
    return render_template("contribute.html", navbar="dropdown")


if __name__ == "__main__":
    app.run(debug=True)
