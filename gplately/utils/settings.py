import logging
import os
from pathlib import Path

logger = logging.getLogger("gplately")

# Settings are resolved from three sources, in order of precedence:
# 1. environment variables
# 2. gplately-settings.toml in the current working directory, if it exists
# 3. settings.toml in ~/.gplately, if it exists
# If ~/.gplately/settings.toml does not exist, it is created with default values.

_DEFAULTS = {
    "disable_dev_warning": False,
    "model_repository": "~/.gplately/plate-models",
    "debug": False,
    "notebook_quick_run": False,
    "test_level": 0,
}

_ENV_VARS = {
    "disable_dev_warning": "GPLATELY_DISABLE_DEV_WARNING",
    "model_repository": "GPLATELY_MODEL_REPOSITORY",
    "debug": "GPLATELY_DEBUG",
    "notebook_quick_run": "GPLATELY_NOTEBOOK_QUICK_RUN",
    "test_level": "GPLATELY_TEST_LEVEL",
}

_GPLATELY_HOME = Path("~/.gplately").expanduser()
_HOME_SETTINGS_FILE = _GPLATELY_HOME / "settings.toml"
_CWD_SETTINGS_FILE = Path("gplately-settings.toml")


def _load_toml(path: Path) -> dict:
    if not path.is_file():
        return {}
    try:
        import tomllib  # Python 3.11+
    except ImportError:
        try:
            import tomli as tomllib  # backport for Python < 3.11
        except ImportError:
            logger.warning(
                f"Unable to read {path} -- reading a settings.toml file requires the "
                "'tomli' package on Python < 3.11. Install it with 'pip install tomli'."
            )
            return {}
    try:
        with open(path, "rb") as f:
            return tomllib.load(f)
    except Exception as exc:
        logger.warning(f"Failed to parse settings file {path}: {exc}")
        return {}


def _format_toml_value(value) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, (int, float)):
        return str(value)
    escaped = str(value).replace("\\", "\\\\").replace('"', '\\"')
    return f'"{escaped}"'


def _create_default_settings_file(path: Path):
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
        lines = [
            f"{key} = {_format_toml_value(value)}" for key, value in _DEFAULTS.items()
        ]
        path.write_text("\n".join(lines) + "\n")
    except OSError as exc:
        logger.warning(f"Unable to create default settings file at {path}: {exc}")


if not _HOME_SETTINGS_FILE.is_file():
    _create_default_settings_file(_HOME_SETTINGS_FILE)

_home_settings = _load_toml(_HOME_SETTINGS_FILE)
_cwd_settings = _load_toml(_CWD_SETTINGS_FILE)


def _get_setting(key: str):
    env_var = os.environ.get(_ENV_VARS[key])
    if env_var is not None:
        return env_var
    if key in _cwd_settings:
        return _cwd_settings[key]
    if key in _home_settings:
        return _home_settings[key]
    return _DEFAULTS[key]


def _get_bool(key: str) -> bool:
    value = _get_setting(key)
    if isinstance(value, str):
        return value.lower() == "true"
    return bool(value)


def _get_str(key: str) -> str:
    return str(_get_setting(key))


def _get_int(key: str) -> int:
    return int(_get_setting(key))


disable_dev_warning = _get_bool("disable_dev_warning")
model_repository = _get_str("model_repository")
debug = _get_bool("debug")
notebook_quick_run = _get_bool("notebook_quick_run")
test_level = _get_int("test_level")
