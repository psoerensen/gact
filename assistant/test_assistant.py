"""Focused static acceptance checks: python -m unittest discover -s assistant."""
import unittest
from html.parser import HTMLParser
from urllib.parse import urlsplit, parse_qs
from build import ROOT, TASKS, render


class Page(HTMLParser):
    def __init__(self, source):
        super().__init__()
        self.links = []
        self.tasks = []
        self.feed(source)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if tag == 'a':
            self.links.append(attrs['href'])
        if tag == 'button':
            self.tasks.append(attrs['data-task'])


class AssistantTests(unittest.TestCase):
    def test_generated_page_and_launch_payloads(self):
        source = render()
        self.assertEqual((ROOT / 'help.html').read_text(encoding='utf-8'), source)
        page = Page(source)
        prompt = (ROOT / 'assistant/instructions.txt').read_text(encoding='utf-8')
        links = [url for url in page.links if url.startswith('https://chatgpt.com/')]
        self.assertEqual(len(links), 5)
        self.assertEqual(page.tasks, [task for _, task in TASKS.values()])
        for url, (_, task) in zip(links, TASKS.values()):
            self.assertEqual(parse_qs(urlsplit(url).query)['q'], [prompt + '\nInitial task: ' + task])
        self.assertNotIn('GACT', source)
        self.assertIn('https://psoerensen.github.io/gtools/support.html', page.links)
        for url in page.links:
            if not urlsplit(url).scheme:
                self.assertTrue((ROOT / url).exists(), url)

    def test_public_workflow_targets(self):
        import re
        prompt = (ROOT / 'assistant/instructions.txt').read_text(encoding='utf-8')
        for path in re.findall(r'https://psoerensen.github.io/gact/(Document/\S+)', prompt):
            self.assertTrue((ROOT / path).is_file(), path)
        for path in ['README.Rmd', 'README.md']:
            self.assertIn('(help.html)', (ROOT / path).read_text(encoding='utf-8'))
