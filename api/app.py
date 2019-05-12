from flask import Flask, jsonify, request

app = Flask(__name__, static_url_path="/static")


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
    Maybe make this a singleton?
    """
    params = []
    params.append("something")
    if "id" in request.args:
        return jsonify("You provided {}.".format(request.args["id"]))


@app.route("/")
def home():
    return app.send_static_file("index.html")


if __name__ == "__main__":
    app.run(debug=True)
