import logging
import sys


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s - %(asctime)s - %(name)s: %(message)s",
        stream=sys.stdout,
    )


setup_logging()
