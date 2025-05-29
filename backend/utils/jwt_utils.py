from flask import current_app
import jwt

def decode_jwt(token):
    try:
        payload = jwt.decode(token, current_app.config["SECRET_KEY"], algorithms=["HS256"])
        return payload
    except Exception:
        return None