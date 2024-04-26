import logging.config
import logging.handlers
import os

import yaml


# configurate the logging utility
def setup_logging():
    cfg_file_path = (
        f"{os.path.dirname(os.path.realpath(__file__))}/../logging_config.yaml"
    )
    if os.path.isfile(cfg_file_path):
        with open(cfg_file_path, "rt") as f:
            config = yaml.safe_load(f.read())
        logging.config.dictConfig(config)

        for name in logging.root.manager.loggerDict:
            logging.getLogger("gplately").debug(f"logger: {name}")
            for h in logging.getLogger(name).handlers:
                logging.getLogger("gplately").debug(h)


def turn_on_debug_logging():
    debug_formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s [%(module)s:%(filename)s:%(lineno)s]"
    )
    gplately_logger = logging.getLogger("gplately")
    gplately_logger.setLevel(logging.DEBUG)
    for h in gplately_logger.handlers:
        h.setLevel(logging.DEBUG)
        h.setFormatter(debug_formatter)
    ptt_logger = logging.getLogger("ptt")
    ptt_logger.setLevel(logging.DEBUG)
    for h in ptt_logger.handlers:
        h.setLevel(logging.DEBUG)
        h.setFormatter(debug_formatter)
