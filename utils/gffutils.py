import os.path
import sqlite3

import gffutils


def init_db(gff_path, db_path):
    if not os.path.exists(db_path):
        db = gffutils.create_db(gff_path, db_path, merge_strategy="create_unique")
    else:
        db = sqlite3.connect(db_path, check_same_thread=False)
        db = gffutils.FeatureDB(db)
    return db
