"""Configuration for the tox-antitargets MCP server (pydantic-settings)."""
from functools import lru_cache
from pathlib import Path
from typing import Literal

from pydantic import Field, model_validator
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
    s3_public_endpoint_url: str = Field(default="")
    s3_access_key: str = Field(default="")
    s3_secret_key: str = Field(default="")
    s3_bucket_name: str = Field(default="")
    s3_region: str = Field(default="us-east-1")
    s3_addressing_style: Literal["auto", "path", "virtual"] = Field(default="path")
    s3_ca_bundle: str = Field(default="")
    s3_allow_local_fallback: bool = Field(default=False)
    s3_url_expiration: int = Field(default=3600, ge=1, le=604800)

    @model_validator(mode="after")
    def validate_s3_configuration(self) -> "Settings":
        core = {
            "s3_endpoint_url": self.s3_endpoint_url.strip(),
            "s3_access_key": self.s3_access_key.strip(),
            "s3_secret_key": self.s3_secret_key.strip(),
            "s3_bucket_name": self.s3_bucket_name.strip(),
        }
        configured = [name for name, value in core.items() if value]
        if configured and len(configured) != len(core):
            missing = ", ".join(name for name, value in core.items() if not value)
            raise ValueError(f"partial S3 configuration; missing: {missing}")
        if (self.s3_public_endpoint_url.strip() or self.s3_ca_bundle.strip()) and not configured:
            raise ValueError("S3 public endpoint/CA bundle requires complete S3 configuration")
        return self

    @property
    def use_s3(self) -> bool:
        return bool(
            self.s3_endpoint_url.strip()
            and self.s3_bucket_name.strip()
            and self.s3_access_key.strip()
            and self.s3_secret_key.strip()
        )

    @property
    def s3_presign_endpoint_url(self) -> str:
        return self.s3_public_endpoint_url.strip() or self.s3_endpoint_url.strip()


@lru_cache
def get_settings() -> Settings:
    return Settings()
