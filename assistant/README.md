# gact AI assistant maintenance

The assistant follows the portable-prompt pattern in the public
[gsuite assistant](https://psoerensen.github.io/gtools/help.html).
`instructions.txt` is the single maintained instruction source. `build.py`
defines four short task starters plus a bug-report/feature-request starter and generates the root `help.html` with
complete ChatGPT prefill links and a readable prompt. `launch.js` copies the
same embedded instructions plus the selected task, with a manual fallback.
The page works without JavaScript through links, download and manual copying.
No network requests occur until the user follows an external link.

After editing instructions or page layout, run from the repository root:

```
python assistant/build.py
python assistant/build.py --check
python -m unittest discover -s assistant
node assistant/test_launch.cjs
```

Python and Node are maintenance checks, not R package or visitor dependencies.
No full website render, package install or data acquisition is needed.
Keep the README.Rmd introduction and generated README.md entry in sync.
Review public documentation links when capabilities change; linked tutorials
do not establish scientific validation or current release availability.

The analysis guidance uses gact data/annotations with the owning gsuite analysis
packages. Request preparation shares the same instructions and conversation,
and links to the existing gtools support service guidance rather than adding a
submission backend. Check that guidance for current hosting availability; drafting
never submits requests or launches development.

Checks cover generated-page drift, all five complete launch payloads, local
documentation links, copy actions and denied-clipboard fallback. They do not
validate model responses, scientific methods, provider sign-in or live prefill
delivery. Long URL prefills may be interrupted by providers or browsers; the
complete copy/download route remains available. Portable instructions guide
behavior; they cannot change permissions of tools in a user's AI account.
Publishing remains a separate, explicitly authorized action.
