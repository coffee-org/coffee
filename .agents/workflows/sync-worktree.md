---
description: Git worktree management for coffee.
---

# Sync Worktree

To review PRs or work on isolated features without altering your active branch workspace:
1. Add a worktree:
   ```bash
   git worktree add ../coffee-worktree framework-dev
   ```
2. Perform work and test builds in `../coffee-worktree`.
3. Prune and remove when finished:
   ```bash
   git worktree remove ../coffee-worktree
   ```
