"""
    Connection to database and
        functions for executing queries

"""


import psycopg2 # for connecting to postgres
from psycopg2.extras import RealDictCursor # results returned as dictionaries
from database.credentials import DB_CONFIG


def get_connection():
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        return conn
    except psycopg2.OperationalError as e:
        print(f"connection failed: {e}")
        return None    

def execute_query(sql: str, params=None):
    conn = None
    cursor = None

    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute(sql, params)
        results = cursor.fetchall()
        conn.commit()

        return [dict(row) for row in results]
        
    except Exception as e:
        print(f"query faild {e}")
        return []
    
    finally:
        if cursor: cursor.close()
        if conn: conn.close()

def execute_insert(sql: str, params=None):
    conn = None
    cursor = None

    try:
        conn = psycopg2.connect(**DB_CONFIG)
        cursor = conn.cursor()
        cursor.execute(sql,params)
        conn.commit()
        try:
            return cursor.fetchone()[0]
        except Exception:
            return None
    
    except Exception as e:
        print(f"query failed: {e}")
        return []
    
    finally:
        if cursor: cursor.close()
        if conn: conn.close()