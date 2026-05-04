import os

DB_CONFIG = {
    "host"     : os.getenv("DB_HOST", "localhost"),
    "database" : os.getenv("DB_NAME", "genomics_db"),
    "user"     : os.getenv("DB_USER", "genomics_user"),
    "password" : os.getenv("DB_PASSWORD", "genomics123"),
    "port"     : os.getenv("DB_PORT", "5432")
}