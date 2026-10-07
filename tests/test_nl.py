"""The intent parser's call to the Anthropic SDK.

Offline: the request shape is checked against the INSTALLED SDK's signature, so a parameter
the SDK drops (temperature went in 1.0) fails here instead of as a TypeError on a live key,
and the JSON schema sent as output_config.format names exactly the keys the prompt asks for.
Live (RUN_LIVE=1 and ANTHROPIC_API_KEY): one real call through parse().
"""
import inspect
import os
import sys
import types
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import app.nl as nl  # noqa: E402


class _FakeMessages:
    def __init__(self, store, text):
        self.store, self.text = store, text

    def create(self, **kwargs):
        self.store.update(kwargs)
        return types.SimpleNamespace(content=[types.SimpleNamespace(type="text", text=self.text)])


def _fake_client(store, text):
    def factory(**_kw):
        return types.SimpleNamespace(messages=_FakeMessages(store, text))
    return factory


def test_request_shape_is_accepted_by_the_installed_sdk(monkeypatch):
    import anthropic
    from anthropic.resources.messages import Messages

    sent = {}
    reply = '{"variant": "rs6025", "gene": "F5", "window_bp": null, "ancestry": null, "common_maf": null}'
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key-never-used")
    monkeypatch.setattr(anthropic, "Anthropic", _fake_client(sent, reply))

    assert nl._llm_intent("the factor V Leiden mutation")["variant"] == "rs6025"

    accepted = inspect.signature(Messages.create).parameters
    unknown = sorted(k for k in sent if k not in accepted)
    assert not unknown, f"messages.create() does not take {unknown} in anthropic {anthropic.__version__}"
    assert "temperature" not in sent
    assert sent["output_config"]["format"]["schema"] is nl._INTENT_SCHEMA


def test_schema_keys_are_the_prompt_keys():
    keys = set(nl._INTENT_SCHEMA["properties"])
    assert keys == set(nl._INTENT_SCHEMA["required"])
    for k in keys:
        assert f'"{k}"' in nl._SYSTEM, f"{k} is in the schema but not named in the prompt"


def test_fenced_reply_is_still_read(monkeypatch):
    import anthropic

    sent = {}
    reply = '```json\n{"variant": "rs6025", "gene": null, "window_bp": null, "ancestry": null, "common_maf": null}\n```'
    monkeypatch.setenv("ANTHROPIC_API_KEY", "test-key-never-used")
    monkeypatch.setattr(anthropic, "Anthropic", _fake_client(sent, reply))
    assert nl._llm_intent("x")["variant"] == "rs6025"


@pytest.mark.skipif(os.environ.get("RUN_LIVE") != "1" or not os.environ.get("ANTHROPIC_API_KEY"),
                    reason="live model call; set RUN_LIVE=1 and ANTHROPIC_API_KEY")
def test_live_free_text_names_a_variant():
    nl._cached_intent.cache_clear()
    r = nl.parse("the factor V Leiden mutation, 1 Mb window")
    assert r.used_llm
    assert nl._RSID.fullmatch(r.query.variant) or nl._HGVS.fullmatch(r.query.variant)
    nl._cached_intent.cache_clear()
