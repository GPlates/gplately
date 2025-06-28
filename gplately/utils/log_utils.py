#
#    Copyright (C) 2024-2025 The University of Sydney, Australia
#
#    This program is free software; you can redistribute it and/or modify it under
#    the terms of the GNU General Public License, version 2, as published by
#    the Free Software Foundation.
#
#    This program is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
#    for more details.
#
#    You should have received a copy of the GNU General Public License along
#    with this program; if not, write to Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
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
    if get_debug_level() > 0:
        turn_on_debug_logging()


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

    gplately_logger.debug("The GPlately debug logging has been turned on.")


def get_debug_level():
    if "GPLATELY_DEBUG" in os.environ:
        if os.environ["GPLATELY_DEBUG"].lower() == "true":
            return 1
        try:
            return int(os.environ["GPLATELY_DEBUG"])
        except:
            return 0
    else:
        return 0
