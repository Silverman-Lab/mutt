# Continue the MUTT release on a laptop

## If you already have the pre-migration repository locally

There is **no published download link** for the cloud handoff ZIP. Its existence
in the cloud filesystem does not make it downloadable from this chat.

If your existing local checkout contains the original processed study files, use
it as the source instead of transferring the cloud archive. In Codex desktop,
open that existing folder and ask the local agent to:

1. Inspect and preserve local changes and the pre-migration checkout at `01c648b`.
2. Fetch `codex/mutt-v1-release` and create a separate sibling worktree for it.
3. Retain/install the old API in an isolated R library, then use
   `scripts/migration/capture-current.R` from the new worktree against the old
   processed study inputs. This loads existing tables; it must not denoise,
   classify, profile, or invoke raw-data processing.
4. Compare captured default objects with the committed `scientific_sha256`
   values in `docs/migration/current-object-comparison.csv`. Investigate any
   mismatch instead of accepting changed scientific data.
5. Use `stage-current.R` and the permanent release builder locally. Archives
   built on another machine may have different ZIP bytes/timestamps, so generate
   the local candidate hashes rather than reusing the cloud asset hashes.
6. Verify old/new complete scientific objects in both forms, upload/download
   the draft release with local gh authentication, and stop for review.

If processed files are missing or are LFS pointers, stop and report which inputs
are missing. Do not silently reprocess studies or substitute historical objects.
This route still needs the old package's parsing dependencies; the separate
already-built-asset route below needs only digest and jsonlite.

## If the already-built cloud archive can be exported

The reviewed work is on branch `codex/mutt-v1-release`, draft PR
https://github.com/Silverman-Lab/mutt/pull/2. The large release assets are **not**
in Git. They are provided separately in `mutt-laptop-handoff.zip`, under
`mutt-laptop-handoff/release-v1/` (about 286 MB). Download that bundle from the
cloud workspace before discarding it. Cloning the branch alone does not transfer
the already-built study assets.

## Open a local working copy

```sh
git clone --branch codex/mutt-v1-release --single-branch --depth 1 https://github.com/Silverman-Lab/mutt.git mutt-release
cd mutt-release
```

Open this `mutt-release` folder as a local project in Codex/your editor. A source
snapshot is also included in the handoff for reference; use the Git checkout for
release commands so the release is tied to a real pushed source commit.

Required software: Git, R (>= 4.0), GitHub CLI (`gh`). For this handoff you need
only two R packages, not the old scientific processing environment:

```sh
Rscript -e 'install.packages(c("digest", "jsonlite"), repos="https://cloud.r-project.org")'
gh auth status
```

If GitHub CLI is not authenticated, run `gh auth login` locally with an account
that can create releases in Silverman-Lab/mutt. Do not transfer a cloud token.

## Restore the existing assets; do not rebuild them

Unzip the handoff. In the local repository create `.build/` and copy the included
`release-v1` folder into it. The result must contain:

    mutt-release/.build/release-v1/studies.candidate.json
    mutt-release/.build/release-v1/vandeputte2017.zip
    ... 39 other study ZIPs ...

Then run from the repository root:

```sh
Rscript scripts/restore-staging.R .build/release-v1 .build/data-v1.0.0
```

This checks each ZIP's pinned byte size/SHA-256, archive members, provenance,
RDS checksums and scientific object loading. It restores staging from the assets
without calling study processing or changing the ZIP bytes. It refuses to
overwrite a staging directory. The handoff also includes SHA256SUMS for its
payload and the original scientific-comparison/check records.

## Create and verify the draft data release

```sh
Rscript scripts/publish-draft.R .build/data-v1.0.0 .build/release-v1
gh release view data-v1.0.0 --repo Silverman-Lab/mutt --web
```

The first command creates/populates a **draft**, downloads every uploaded study,
and checks it against staging. It never overwrites a published release or
replaces an existing asset. Success writes
`.build/release-v1/draft-verification.json`. Stop at the verified draft for review;
this handoff does not authorize automatic publication or merging.

After explicit publication, generate the public manifest by verifying all public
URLs. This remains a separate step:

```sh
Rscript scripts/verify-public.R .build/data-v1.0.0 .build/release-v1 inst/extdata/studies.json
```

Commit/push the resulting available manifest only after that verification succeeds.
Other unavailable study records are retained automatically.

## Context for a local coding agent

Read AGENTS.md, docs/migration/FINAL.md and docs/RELEASING.md first. The baseline
is current mutt at `01c648b`, as clarified by the user. Forty working studies
match the old API exactly in both default and aligned forms. Seventeen existing
parsers failed in both forms; six other directories were unregistered source
records. Do not substitute historical objects for the current results.

All 40 default and 40 aligned objects passed exact end-to-end loader comparisons.
953 package fixture assertions and 5 release-tool assertions passed. R CMD check
returned Status OK on R 4.5.3. A Vandeputte 2021 recipe smoke check also reproduced
the old aligned object exactly. Scientific preprocessing was not rerun. The empty
functional-inference field was deliberately removed; other scientific contents
were retained. The 37-study historical object was investigated but is not the
selected baseline.

The cloud network proxy blocked api.github.com/uploads.github.com, preventing
data-release creation there. No data release upload has succeeded. The pinned
package manifest therefore says staged, and normal first downloads are not live.
Finish the draft upload/download verification locally. Preserve existing asset
bytes; do not invoke migration, DADA2, MetaPhlAn, or build-release to repeat work.
