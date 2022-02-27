import logging

from flask import Flask, _app_ctx_stack, jsonify
from sqlalchemy.orm import scoped_session

from regiosqm_api import database, models

logging.basicConfig(level=logging.DEBUG)

_logger = logging.getLogger(__name__)

app = Flask(__name__, static_folder="public", static_url_path="/public")
app.session = scoped_session(database.SessionLocal, _app_ctx_stack.__ident_func__)


@app.route("/heartbeat")
def heartbeat():
    return jsonify({"status": "healthy"})


@app.route("/submit")
def work():
    return jsonify({"status": "healthy"})


@app.route("/records")
def get_records():
    records = app.session.query(models.Prediction).all()
    _logger.info(len(records))
    for record in records:
        _logger.info(record)
    return jsonify()


@app.route("/", defaults={"path": ""})
@app.route("/<path:path>")
def catch_all(path):
    if path == "":
        path = "index.html"
    return app.send_static_file(path)


def main(debug=False):
    # db.create_all()
    app.run(debug=debug)


@app.teardown_appcontext
def remove_session(*args, **kwargs):
    app.session.remove()


if __name__ == "__main__":
    main(debug=True)
