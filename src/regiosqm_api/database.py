from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from regiosqm_api.models import Base


def initialize_db(engine):
    Base.metadata.create_all(engine)


SQLALCHEMY_DATABASE_URL = "sqlite:///./database.db"
engine = create_engine(SQLALCHEMY_DATABASE_URL, connect_args={"check_same_thread": False})
initialize_db(engine)  # TODO Only if not-created
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
