import flask
from flask import jsonify, request

app = flask.Flask(__name__)
app.config["DEBUG"] = True

# TODO: call the methods here


@app.route("/", methods=["GET"])
def home():
    return """
<h1>Hello Python!</h1>
"""


@app.route("/api/v1/", methods=["GET"])
def api():
    """
    The following parameters are extracted from the URL:
        - alpha
        - beta
        - mu
        - t_start=0.0
        - t_end=30.0
        - N=(t_start - t_end)*100
        - nonlinear=True
        - phi0
        - phi0_dot
        - max_control
        - frequency
        - deadband
        - set_point
        - precision
    """
    params = []
    params.append("something")
    if "id" in request.args:
        return jsonify("You provided {}.".format(request.args["id"]))


if __name__ == "__main__":
    app.run()
