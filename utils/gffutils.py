import os.path

import gffutils


def init_db(gff_path, db_path):
    if not os.path.exists(db_path):
        db = gffutils.create_db(gff_path, db_path, merge_strategy="create_unique")
    else:
        db = gffutils.FeatureDB(db_path)
    return db
