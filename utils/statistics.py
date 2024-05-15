import statistics


def cov(array):
    """Coefficient of Variation"""
    m = statistics.mean(array)
    s = statistics.stdev(array)
    return s / m
