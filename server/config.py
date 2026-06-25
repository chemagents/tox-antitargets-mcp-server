"""Configuration for the tox-antitargets MCP server (pydantic-settings)."""
from functools import lru_cache
from pathlib import Path

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict

PACKAGE_DIR = Path(__file__).resolve().parent
PROJECT_DIR = PACKAGE_DIR.parent
DEFAULT_ENV_FILE = PROJECT_DIR / ".env"
DEFAULT_DATA_FILE = PACKAGE_DIR / "data" / "antitargets_LD50_affinity.csv"
DEFAULT_ARTIFACTS_DIR = PROJECT_DIR / "artifacts"

DATASET_URL = (
    "https://raw.githubusercontent.com/chemagents/ld50-antitargets/"
    "main/antitargets_LD50_affinity.csv"
)


class Settings(BaseSettings):
    """Server settings. Values are read from environment / .env with TOX_ prefix."""

    model_config = SettingsConfigDict(
        env_prefix="TOX_",
        env_file=str(DEFAULT_ENV_FILE),
        env_file_encoding="utf-8",
        extra="ignore",
        case_sensitive=False,
    )

    # --- MCP transport (matches the other CoScientist MCP servers) ---
    mcp_host: str = Field(default="0.0.0.0")
    mcp_port: int = Field(default=7331)
    mcp_path: str = Field(default="/mcp")

    # --- dataset ---
    dataset_path: str = Field(default=str(DEFAULT_DATA_FILE))
    dataset_url: str = Field(default=DATASET_URL)

    # --- analysis parameters (defaults reproduce the paper) ---
    binder_threshold: float = Field(default=-7.0)        # kcal/mol; paper section 3.3
    tanimoto_threshold: float = Field(default=0.65)      # Butina; paper section 2.6
    morgan_nbits: int = Field(default=2048)

    # --- artifacts (figures): local dir or S3-compatible bucket ---
    artifacts_dir: str = Field(default=str(DEFAULT_ARTIFACTS_DIR))
    artifact_url_base: str = Field(default="")           # optional public prefix for local files
    s3_endpoint_url: str = Field(default="")
    s3_access_key: str = Field(default="")
    s3_secret_key: str = Field(default="")
    s3_bucket_name: str = Field(default="")
    s3_url_expiration: int = Field(default=3600)

    @property
    def use_s3(self) -> bool:
        return bool(self.s3_endpoint_url and self.s3_bucket_name and self.s3_access_key)


@lru_cache
def get_settings() -> Settings:
    return Settings()
