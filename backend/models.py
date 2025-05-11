from flask_login import UserMixin
from flask_bcrypt import Bcrypt
from bson import ObjectId

bcrypt = Bcrypt()

class User(UserMixin):
    def __init__(self, _id, username, password):
        self.id = str(_id)  # Flask-Login needs string type for self.id
        self.username = username
        self.password = password

    @staticmethod
    def get_user_by_username(mongo, username):
        user_data = mongo.db.users.find_one({'username': username})
        if user_data:
            return User(user_data['_id'], user_data['username'], user_data['password'])
        return None

    @staticmethod
    def get_user_by_id(mongo, user_id):
        user_data = mongo.db.users.find_one({'_id': ObjectId(user_id)})
        if user_data:
            return User(user_data['_id'], user_data['username'], user_data['password'])
        return None

    @staticmethod
    def create_user(mongo, username, password):
        if mongo.db.users.find_one({'username': username}):
            return None  # L'utilisateur existe déjà
        hashed_password = bcrypt.generate_password_hash(password).decode('utf-8')
        result = mongo.db.users.insert_one({
            'username': username,
            'password': hashed_password
        })
        return str(result.inserted_id)
