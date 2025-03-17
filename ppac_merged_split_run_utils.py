import os.path
import pickle

from reannotation.pipelines import suspicious_orthologue_pipeline


def pickle_cache_suspicious_orthologue_pipeline(tool, *args, **kwargs):
    merged_path = os.path.join("data", "tmp", "ppac_{}_merged.pickle".format(tool))
    split_path = os.path.join("data", "tmp", "ppac_{}_split.pickle".format(tool))
    if os.path.isfile(merged_path) and os.path.isfile(split_path):
        with open(merged_path, "rb") as f:
            merged = pickle.load(f)
        with open(split_path, "rb") as f:
            split = pickle.load(f)
    else:
        merged, split = suspicious_orthologue_pipeline(*args, **kwargs)
        with open(merged_path, 'wb') as f:
            pickle.dump(merged, f, protocol=pickle.HIGHEST_PROTOCOL)
        with open(split_path, 'wb') as f:
            pickle.dump(split, f, protocol=pickle.HIGHEST_PROTOCOL)
    return merged, split
