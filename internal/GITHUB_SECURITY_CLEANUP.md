# GitHub security cleanup — IRB raw data purge

## What happened
The in-house orbital RNA-seq raw counts `data.txt` (15.3 MB, 46,428 rows,
per-sample read counts) were committed in the initial commit `31670fd` and
removed in `361ced6`. Deletion alone left the file recoverable from git history.

## What was done (✅ complete)
1. Safety backup of all refs taken before rewriting (local bundle, not pushed).
2. `git filter-repo --invert-paths --path data.txt` removed the blob
   (`1ead1f5…`) from **all** commits on both `main` and the working branch.
3. Verified: `data.txt` no longer present in any object in history.
4. Force-pushed rewritten history: `main` (88eabcb → 4abfb6d) and the working
   branch. Tracked file count unchanged (529); `data.txt` was never in the
   current tree, only in history.

## Verified status (2026-06-01)
- Branches on GitHub are **clean**: `main` → `4abfb6d`, `claude/brave-sagan-R84Lj`
  → `01653dd` (both on purged history; no `data.txt`).
- The old commit `88eabcb` is **still reachable by direct SHA** via the GitHub
  API/web — confirmed today. So `data.txt` remains retrievable server-side until
  GitHub garbage-collects it. **Action below is still required.**
- Note: the final v5 work lives on branch `claude/brave-sagan-R84Lj`; `main` holds
  only the purged original (no submission package). Decide the final `main` state
  (merge the branch) before or as part of the recreate.

## ⚠️ Remaining action — REQUIRED, only the repo owner can do this
A force-push does **not** immediately erase the old commit from GitHub's servers.
The dangling commit `88eabcb` (and therefore `data.txt`) can still be retrieved
by direct SHA on GitHub until GitHub garbage-collects it. There are **no open
PRs and no forks**, so no ref pins it — but to be safe and definitive, do **one**
of the following:

- **Option 1 (most reliable): delete and recreate the GitHub repo.**
  Private repo, owner-only, no PRs/forks → low cost. Steps:
  1. GitHub → repo → Settings → Danger Zone → *Delete this repository*.
  2. Recreate an empty private `TEDGWAS`.
  3. From a clean local clone of the purged history: `git push -u origin --all`.
  This guarantees the old blob is gone server-side.

- **Option 2: contact GitHub Support.** Ask them to purge cached views /
  garbage-collect unreachable objects after a force-push to remove sensitive
  data. Reference: repo `JungyulPark/TEDGWAS`, commit `88eabcb`, blob `1ead1f5`.

## Going forward
- `.gitignore` blocks `data.txt`, `*.gz`, FASTQ/BAM, patent files. Do not weaken it.
- Never `git add -f` IRB or licensed data.
