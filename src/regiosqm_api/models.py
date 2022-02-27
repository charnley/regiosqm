import gzip
import io

import numpy as np
import sqlalchemy as sa
from sqlalchemy import Column, DateTime, Integer, LargeBinary, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.schema import MetaData

# Recommended naming convention used by Alembic, as various different database
# providers will autogenerate vastly different names making migrations more
# difficult. See: http://alembic.zzzcomputing.com/en/latest/naming.html
NAMING_CONVENTION = {
    "ix": "ix_%(column_0_label)s",
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(constraint_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s",
}

metadata = MetaData(naming_convention=NAMING_CONVENTION)
Base = declarative_base(metadata=metadata)


class CompressedString(sa.types.TypeDecorator):
    """Storage datatype for large blobs of text"""

    @staticmethod
    def compress(s):
        if type(s) == str:
            s = s.encode()
        b = gzip.compress(s)
        return b

    @staticmethod
    def decompress(b):
        s = gzip.decompress(b)
        return s

    impl = LargeBinary

    def process_bind_param(self, value, dialect):
        return self.compress(value)

    def process_result_value(self, value, dialect):
        return self.decompress(value)


class NumpyArray(sa.types.TypeDecorator):
    """Storage datatype for numpy arrays"""

    impl = LargeBinary

    @staticmethod
    def compress(s):
        if type(s) == str:
            s = s.encode()
        b = gzip.compress(s)
        return b

    @staticmethod
    def decompress(b):
        s = gzip.decompress(b)
        return s

    @staticmethod
    def save_array(arr):
        s = io.StringIO()
        np.savetxt(s, arr)
        return s.getvalue()

    @staticmethod
    def load_array(txt):
        s = io.StringIO(txt)
        arr = np.loadtxt(s)
        return arr

    def process_bind_param(self, value, dialect):
        value = self.save_array(value)
        return self.compress(value)

    def process_result_value(self, value, dialect):
        value = self.decompress(value)
        value = self.load_array(value)
        return value


class Prediction(Base):
    __tablename__ = "regiopredictions"
    id = Column(Integer, primary_key=True)

    # Basic descriptors
    hashkey = Column(String)
    created = Column(DateTime)
    name = Column(String)
    smiles = Column(String)
    sdf = Column(String)
    method = Column(String)
