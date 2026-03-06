"""
database/connect.py
====================
Central database connection module.
Every other script in this project imports from here.

Usage:
    from database.connect import get_connection, execute_query
"""

import psycopg2
from psycopg2.extras import RealDictCursor

# ── Connection config ────────────────────────────────────────
DB_CONFIG = {
    "host"     : "localhost",
    "database" : "genomics_db",
    "user"     : "genomics_user",
    "password" : "genomics123",
    "port"     : "5432"
}


def get_connection():
    """
    Create and return a raw psycopg2 connection.
    Use this when you need full control over commits/rollbacks.

    Example:
        conn = get_connection()
        cursor = conn.cursor()
        cursor.execute("SELECT * FROM variants")
        conn.commit()
        conn.close()
    """
    try:
        conn = psycopg2.connect(**DB_CONFIG)
        return conn
    except psycopg2.OperationalError as e:
        print(f"❌ Database connection failed: {e}")
        print("   Is PostgreSQL running? Try: sudo systemctl start postgresql")
        raise


def execute_query(sql, params=None, fetch=True):
    """
    Run a SQL query and return results as a list of dictionaries.
    Each row is a dict like: {'gene_name': 'BRCA1', 'position': 41244429}

    Args:
        sql    : The SQL query string
        params : Optional tuple of parameters (use %s in query)
        fetch  : If True, return rows. If False, just execute (for INSERT/UPDATE)

    Example:
        rows = execute_query(
            "SELECT * FROM variants WHERE gene_name = %s",
            params=("BRCA1",)
        )
    """
    conn   = None
    cursor = None
    try:
        conn   = psycopg2.connect(**DB_CONFIG)
        # RealDictCursor returns rows as dicts instead of tuples
        cursor = conn.cursor(cursor_factory=RealDictCursor)
        cursor.execute(sql, params)

        if fetch:
            results = cursor.fetchall()
            # Convert from RealDict to plain list of dicts
            return [dict(row) for row in results]
        else:
            conn.commit()
            return cursor.rowcount   # number of rows affected

    except psycopg2.Error as e:
        if conn:
            conn.rollback()
        print(f"❌ Query failed: {e}")
        print(f"   SQL: {sql}")
        raise

    finally:
        # Always close cursor and connection even if error occurs
        if cursor:
            cursor.close()
        if conn:
            conn.close()


def execute_insert(sql, params=None):
    """
    Run an INSERT and return the new row's ID.
    Requires your INSERT to end with: RETURNING <id_column>

    Example:
        user_id = execute_insert(
            "INSERT INTO users (username, email) VALUES (%s, %s) RETURNING user_id",
            params=("younes", "younesadi18@gmail.com")
        )
    """
    conn   = None
    cursor = None
    try:
        conn   = psycopg2.connect(**DB_CONFIG)
        cursor = conn.cursor()
        cursor.execute(sql, params)
        new_id = cursor.fetchone()[0]
        conn.commit()
        return new_id

    except psycopg2.Error as e:
        if conn:
            conn.rollback()
        print(f"❌ Insert failed: {e}")
        raise

    finally:
        if cursor:
            cursor.close()
        if conn:
            conn.close()


def test_connection():
    """
    Quick sanity check — verifies connection and prints table row counts.
    Run this any time you want to confirm the database is healthy.
    """
    print("Testing database connection...\n")
    try:
        rows = execute_query("""
            SELECT
                table_name,
                (SELECT COUNT(*) FROM information_schema.columns
                 WHERE table_name = t.table_name
                 AND table_schema = 'public') as column_count
            FROM information_schema.tables t
            WHERE table_schema = 'public'
            ORDER BY table_name;
        """)

        print(f"✅ Connected to genomics_db\n")
        print(f"  {'Table':<30} {'Columns'}")
        print(f"  {'-'*40}")
        for row in rows:
            print(f"  {row['table_name']:<30} {row['column_count']}")

        # Print row counts for the two knowledge base tables
        print()
        for table in ['clinvar_annotations', 'pharmgkb_annotations', 'variants']:
            count = execute_query(f"SELECT COUNT(*) as count FROM {table}")
            print(f"  {table:<30} {count[0]['count']:,} rows")

    except Exception as e:
        print(f"❌ Connection test failed: {e}")


if __name__ == "__main__":
    test_connection()