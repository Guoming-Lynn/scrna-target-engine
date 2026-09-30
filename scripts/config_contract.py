"""Required quicklook config keys. The batch-column default lives only here."""

REQUIRED_CONFIG = ("TARGET_GENE", "INPUT_PATH", "OUTPUT_ROOT")
DEFAULT_BATCH_COLUMN = "sample_id"


def missing_config(config: dict[str, object]) -> list[str]:
    missing: list[str] = []
    for key in REQUIRED_CONFIG:
        value = config.get(key)
        if not value or str(value).startswith("YOUR_"):
            missing.append(key)
    return missing


def batch_column(config: dict[str, object]) -> str:
    """Return the configured batch column, or the single shared default."""
    value = config.get("BATCH_COLUMN", DEFAULT_BATCH_COLUMN)
    if not isinstance(value, str) or not value.strip():
        return DEFAULT_BATCH_COLUMN
    return value
