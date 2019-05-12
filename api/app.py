#!/bin/python
#
# Author: Daniel Schuette
# License: GPL-3.0
# Date: 12/05/2019
from flask import Flask, jsonify, request

app = Flask(__name__, static_url_path="/static")


# this is the API route that retrieves the pendulum data, required URL
# parameters are described below
@app.route("/api/v1/", methods=["GET"])
def api():
    """
    The following parameters are extracted from the URL:
        - alpha, beta, mu
        - t_start=0.0, t_end=30.0
        - N=(t_start - t_end)*100
        - nonlinear=True
        - phi0, phi0_dot
        - max_control, frequency, deadband, set_point, precision
        - API key
    TODO: Maybe make this a singleton? Get the API key from the environment.
    """
    params = []
    params.append("something")
    if "id" in request.args:
        return jsonify("You provided {}.".format(request.args["id"]))


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
