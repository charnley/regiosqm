import logging
from functools import wraps

import httpagentparser
from flask import Flask, _app_ctx_stack, jsonify, request
from flask.helpers import make_response
from sqlalchemy.orm import scoped_session

from regiosqm_api import database, models

logging.basicConfig(level=logging.DEBUG)

_logger = logging.getLogger(__name__)

app = Flask(__name__, static_folder="public", static_url_path="/public")
app.session = scoped_session(database.SessionLocal, scopefunc=_app_ctx_stack)


def browser_required(f):
    """Security by obscurity"""

    @wraps(f)
    def decorator(*args, **kwargs):
        browsers = [
            "camino",
            "chrome",
            "firefox",
            "galeon",
            "kmeleon",
            "konqueror",
            "links",
            "lynx",
            "msie",
            "msn",
            "netscape",
            "opera",
            "safari",
            "seamonkey",
            "webkit",
        ]

        headeragent = request.headers.get("User-Agent")
        agent = httpagentparser.detect(headeragent)

        _logger.debug(f"agent: {agent}")

        try:
            agent_browser = agent["browser"]["name"]
        except KeyError:
            agent_browser = None

        if agent_browser is None or agent_browser not in browsers:
            return make_response(
                jsonify({"message": "Invalid token! Use browser interface only"}), 401
            )

        return f(*args, **kwargs)

    return decorator


@app.route("/status")
@browser_required
def heartbeat():

    headeragent = request.headers.get("User-Agent")
    agent = httpagentparser.detect(headeragent)

    try:
        agent_browser = agent["browser"]["name"]
    except KeyError:
        agent_browser = None

    return jsonify(
        {
            "status": "healthy",
            "user-agent": str(headeragent),
            "browser": str(agent_browser),
        }
    )


@app.route("/submit")
@browser_required
def work():
    return jsonify({"status": "healthy"})


@app.route("/records")
@browser_required
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
