import os


# OPENFOLD_IGNORE_IMPORT=1 ignores all import to reduce I/O time
if str(os.getenv("OPENFOLD_IGNORE_IMPORT")) != "1":
    from . import model
    from . import utils
    from . import np
    from . import resources

    __all__ = ["model", "utils", "np", "data", "resources"]

else:
    __all__ = [""]
