from flask import Flask, jsonify

app = Flask(__name__, static_folder="public", static_url_path="/public")


@app.route("/heartbeat")
def heartbeat():
    return jsonify({"status": "healthy"})


@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def catch_all(path):
    if path == "":
        path = "index.html"
    return app.send_static_file(path)


def main(debug=False):
    # db.create_all()
    app.run(debug=debug)


if __name__ == "__main__":
    main(debug=True)
