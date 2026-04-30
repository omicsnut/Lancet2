# /sync-cost-model — Refresh cost-model.md against current Anthropic docs

Refresh `cost-model.md` to reflect current Claude Code documentation, current community-converged practice, and the current state of this bundle. Run this on a quarterly cadence, or whenever you suspect the cost model has drifted (Anthropic released a major Claude Code update, the bundle's structure changed substantively, or numbers in the document feel obviously stale).

## Why this is a slash command and not an agent

Slash commands are user-initiated and cost essentially nothing baseline. This refresh is a maintenance ritual that should happen rarely and deliberately — that's exactly the slash-command pattern. An always-loaded agent for a task invoked maybe quarterly would pay session-baseline cost on every session for one rare benefit; the math doesn't work. If usage shows the slash command isn't enough (e.g., you find yourself wanting to refresh on demand and can't remember the command name), upgrading to an agent is a one-file change.

## Procedure

The refresh has five phases. Do not skip phase 4 — the bundle's own numbers are the most likely thing to be stale.

### Phase 1 — Re-fetch Anthropic guidance

Use `web_fetch` to read each of these in turn:

- `https://code.claude.com/docs/en/costs` — the costs and CLAUDE.md guidance
- `https://code.claude.com/docs/en/sub-agents` — subagent semantics and description-as-trigger
- `https://code.claude.com/docs/en/skills` — skills lifecycle, auto-compaction behavior
- `https://code.claude.com/docs/en/hooks` — hook events and exit-code semantics

For each, note: did the guidance change? Was a new mechanism added? Was a published number updated (token budgets, line targets, etc.)? Are there new caveats or recommendations?

### Phase 2 — Sample current community practice

Use `web_search` to look for substantive community writeups on Claude Code bundle structure published since the last refresh of cost-model.md. Useful queries:

- `Claude Code skills agents hooks best practices 2026` (or current year)
- `CLAUDE.md size optimization tokens`
- `Claude Code subagent vs skill when to use`

Note any consensus shifts. The community is the leading indicator for emerging patterns; Anthropic docs lag by a quarter or two.

### Phase 3 — Diff against cost-model.md

Read the current `cost-model.md`. For each section, ask:

- Does the "What costs what" section accurately describe how each mechanism is loaded? Did Anthropic change any of this?
- Does the "Anthropic-published guidance and citations" section still match what's at those URLs? Are the quoted phrases still verbatim?
- Are the numbers in the "Concrete numbers" table still in the right ballpark?
- Does the "How the bundle's choices map to this" section reflect the bundle's current contents (count of agents, skills, slash commands, hooks)?

Compute a list of edits needed. Be conservative — only edit what is genuinely stale or wrong; do not rewrite for style.

### Phase 4 — Re-measure the bundle's own baseline

Count the bundle's current contents:

```bash
ls .claude/agents/ | grep -v README | grep '\.md$' | wc -l       # subagent count
find .claude/skills -maxdepth 1 -mindepth 1 -type d | wc -l       # skill count (folders)
ls .claude/commands/*.md | grep -v README | wc -l                # slash command count
ls .claude/hooks/*.{py,sh} 2>/dev/null | wc -l                   # hook count
wc -l AGENTS.md                                                  # canonical content size
wc -l CLAUDE.md                                                  # wrapper size (should be ~3 lines)
```

Update the "How the bundle's choices map to this" section to reference these current counts. If the numbers have drifted significantly from what's in the doc, that itself is signal — either the bundle has grown without good reason (prune), or the doc is just stale (update).

For the rough-token-cost table, you generally do not need to re-measure precise tokens. The table is for planning, not budgeting. Update only if the relative ordering of mechanisms has changed (e.g., if a new mechanism was added that fits between two existing rows).

### Phase 5 — Apply edits and commit

Apply the edits identified in phase 3 and the bundle measurement from phase 4. Commit as a single small commit:

```
chore: refresh cost-model.md

Re-verified Anthropic costs/sub-agents/skills/hooks docs.
Updated bundle counts: <numbers>.
<one-line note about anything substantive that changed>
```

If nothing substantive changed (the doc is still accurate), say so in chat and skip the commit. The doc not changing for a quarter is a healthy signal, not a problem.

## When NOT to use this command

Do not use this command if you have not made bundle structural changes recently and have not seen Anthropic announce major Claude Code updates. There is no value in refreshing a document whose subject has not changed. Quarterly cadence is the upper bound; less often is fine.

Do not use this command as a way to reorganize cost-model.md for style. The command is for accuracy, not for prose improvement. Style edits go in their own commit.

## Maintenance

This command's procedure is itself subject to drift — if Anthropic moves docs, adds a new mechanism, or deprecates a section, this command's phase-1 URL list will need updating. When you find yourself working around stale instructions in this command, edit them.

The "useful queries" list in phase 2 should grow over time. As you find queries that consistently surface useful results, add them; as you find queries that turn out to be noise, remove them.
