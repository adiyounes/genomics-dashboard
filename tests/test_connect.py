import sys
from pathlib import Path
sys.path.append(str(Path(__file__).parent.parent))

from database.connect import get_connection, execute_query

def test_get_connection():
    conn = get_connection()
    assert conn is not None
    conn.close()

def test_execute_query_returns_list():
    assert isinstance(execute_query("SELECT 1"), list)