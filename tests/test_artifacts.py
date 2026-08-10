"""Artifact storage configuration, failure semantics, and reference contract."""
from __future__ import annotations

import hashlib
from pathlib import Path
from types import SimpleNamespace

import pytest
from botocore.exceptions import ClientError, EndpointConnectionError
from pydantic import ValidationError

from server import artifacts
from server.config import Settings


@pytest.fixture(autouse=True)
def _clear_clients():
    artifacts.clear_client_cache()
    yield
    artifacts.clear_client_cache()


def _settings(tmp_path: Path, *, allow_fallback: bool = False):
    return SimpleNamespace(
        use_s3=True,
        artifacts_dir=str(tmp_path),
        artifact_url_base="",
        s3_endpoint_url="http://minio:9000",
        s3_presign_endpoint_url="http://localhost:9000",
        s3_access_key="artifact-writer",
        s3_secret_key="test-secret",
        s3_bucket_name="coscientist-artifacts",
        s3_region="us-east-1",
        s3_addressing_style="path",
        s3_ca_bundle="",
        s3_allow_local_fallback=allow_fallback,
        s3_url_expiration=3600,
    )


def test_partial_s3_configuration_is_rejected():
    with pytest.raises(ValidationError, match="partial S3 configuration"):
        Settings(
            _env_file=None,
            s3_endpoint_url="http://minio:9000",
            s3_access_key="writer",
            s3_secret_key="",
            s3_bucket_name="bucket",
        )


def test_internal_upload_and_public_presign_clients_are_separate(tmp_path, monkeypatch):
    settings = _settings(tmp_path)
    monkeypatch.setattr(artifacts, "get_settings", lambda: settings)
    calls: list[tuple[str, dict]] = []

    class FakeClient:
        def __init__(self, endpoint: str):
            self.endpoint = endpoint

        def put_object(self, **kwargs):
            calls.append(("put", kwargs))

        def generate_presigned_url(self, _method, *, Params, ExpiresIn):
            calls.append(("presign", {"Params": Params, "ExpiresIn": ExpiresIn}))
            return f"{self.endpoint}/{Params['Bucket']}/{Params['Key']}?signed=1"

    endpoints: list[str] = []

    def fake_client(_service, **kwargs):
        endpoints.append(kwargs["endpoint_url"])
        return FakeClient(kwargs["endpoint_url"])

    monkeypatch.setattr(artifacts.boto3, "client", fake_client)
    data = b"\x89PNG\r\n\x1a\nstrict-storage-test"
    result = artifacts.store_png(data, "roundtrip")

    assert endpoints == ["http://minio:9000", "http://localhost:9000"]
    assert result["kind"] == "s3"
    assert result["artifact"].startswith("http://localhost:9000/")
    assert result["key"].startswith("tox_antitargets/")
    assert result["sha256"] == hashlib.sha256(data).hexdigest()
    put = next(payload for kind, payload in calls if kind == "put")
    assert put["Body"] == data
    assert put["ContentType"] == "image/png"
    assert put["Metadata"]["sha256"] == result["sha256"]


def test_s3_upload_failure_is_strict_by_default(tmp_path, monkeypatch):
    settings = _settings(tmp_path)
    monkeypatch.setattr(artifacts, "get_settings", lambda: settings)

    class FailingClient:
        def put_object(self, **_kwargs):
            raise EndpointConnectionError(endpoint_url=settings.s3_endpoint_url)

    monkeypatch.setattr(artifacts, "_upload_client", lambda: FailingClient())
    with pytest.raises(artifacts.ArtifactStorageError, match="failed to store PNG"):
        artifacts.store_png(b"png", "strict")
    assert not list(tmp_path.iterdir())


def test_s3_local_fallback_requires_opt_in(tmp_path, monkeypatch):
    settings = _settings(tmp_path, allow_fallback=True)
    monkeypatch.setattr(artifacts, "get_settings", lambda: settings)

    class FailingClient:
        def put_object(self, **_kwargs):
            raise EndpointConnectionError(endpoint_url=settings.s3_endpoint_url)

    monkeypatch.setattr(artifacts, "_upload_client", lambda: FailingClient())
    result = artifacts.store_png(b"png", "fallback")
    assert result["kind"] == "local"
    assert Path(result["artifact"]).exists()
    assert "explicit local fallback" in result["warning"]


def test_missing_bucket_fails_startup(tmp_path, monkeypatch):
    settings = _settings(tmp_path)
    monkeypatch.setattr(artifacts, "get_settings", lambda: settings)

    class MissingBucketClient:
        def head_bucket(self, **_kwargs):
            raise ClientError({"Error": {"Code": "404", "Message": "Not Found"}}, "HeadBucket")

    monkeypatch.setattr(artifacts, "_upload_client", lambda: MissingBucketClient())
    with pytest.raises(artifacts.ArtifactStorageError, match="not ready"):
        artifacts.ensure_storage_ready()


def test_unwritable_local_dir_falls_back_with_explicit_warning(tmp_path, monkeypatch):
    blocked = tmp_path / "regular-file-not-directory"
    blocked.write_text("blocks mkdir", encoding="utf-8")
    settings = SimpleNamespace(
        use_s3=False,
        artifacts_dir=str(blocked),
        artifact_url_base="https://would-be-wrong-for-the-temp-file.invalid",
    )
    monkeypatch.setattr(artifacts, "get_settings", lambda: settings)
    monkeypatch.setattr(artifacts.tempfile, "gettempdir", lambda: str(tmp_path))

    result = artifacts.store_png(b"png", "local-fallback")
    artifact = Path(result["artifact"])
    assert artifact.exists()
    assert artifact.parent == tmp_path / "tox_antitargets"
    assert "temp storage" in result["warning"]
    assert not result["artifact"].startswith("https://")
