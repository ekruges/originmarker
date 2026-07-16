# Enabling free-text search

Free-text search ("markers near the ABCC8 splice mutation in Europeans") needs an
Anthropic API key. Everything else works without one: rsID and HGVS input, panel builds,
every export. This is the only feature the key buys.

## Where the key goes

```bash
cd /path/to/originmarker
cp .env.example .env
vi .env                      # paste the key on the ANTHROPIC_API_KEY= line
chmod 600 .env               # owner-only
docker compose up -d         # compose reads .env and injects it at run time
```

Confirm it took, without printing the key:

```bash
curl -s https://ezrakruger.cc/originmarker/api/health | grep -o '"nl_enabled":[a-z]*'
# "nl_enabled":true
```

## Why the key cannot leak

Not by discipline. By construction:

- **It cannot reach the browser.** The frontend is Vite, and Vite only inlines variables
  prefixed `VITE_`. No `VITE_*` variable exists anywhere in this project, so there is no
  mechanism by which a secret could be compiled into the bundle even by mistake. Verify:
  `grep -r "VITE_" web/src/` returns nothing.
- **It is read server-side only**, in `app/nl.py`, inside the container.
- **The API exposes a boolean, never the value.** `/api/health` returns
  `nl_enabled: true|false`. There is no endpoint that returns the key.
- **It cannot be committed.** `.env`, `.env.*`, `*.pem` and `credentials*.json` are in
  `.gitignore`. Verify with `git check-ignore -v .env`.
- **It cannot be baked into an image.** The same patterns are in `.dockerignore`, so the
  key is never part of a layer. It arrives at `docker compose up` time as process
  environment. The image stays clean and is safe to rebuild or push.

The same applies to `LDLINK_TOKEN` and `NCBI_API_KEY`.

## What it costs

Close to nothing, because the model is a last resort rather than the front door.

Any input containing an rsID or an HGVS expression is parsed by regex and **never reaches
the API**. That is the overwhelming majority of real queries. Window, MAF and ancestry
modifiers are also parsed locally, even on the LLM path, because a regex is free and
deterministic. The model is called only for prose that names no identifier at all: Haiku,
temperature 0, `max_tokens` capped, system prompt cached.

You can confirm the fast path is free from the response itself: `used_llm` is `false`
whenever the regex handled it.

## What the model is allowed to do

It resolves **intent only**. It picks which variant you meant; it never supplies a
position. `pb.StructuredQuery`, the only shape the parser can produce, has no `chrom`,
`pos`, `ref`, `alt` or `strand` field, so there is no channel through which a model could
pass a coordinate even if it hallucinated one. Every genomic fact still comes from the
same live lookup the manual path uses, and `app/nl.py` rejects coordinate-shaped output at
runtime rather than trusting the prompt to have forbidden it.

## The failure mode that remains

The model cannot invent a position, but it **can pick the wrong variant** from ambiguous
wording. If it does, every number downstream will be correct about the wrong locus, which
is a worse kind of wrong than an obvious error, because nothing looks broken.

That is why any model-parsed query renders a red warning in the interface telling the
reader to check the resolved variant before building a panel, and why the resolved-variant
card is shown first, with the rsID, coordinate and ClinVar link, before the expensive pull
runs. Read the card, not the search box.

If that risk is not acceptable for a given piece of work, leave the key unset and use
**Manual input**, which does no parsing at all.
