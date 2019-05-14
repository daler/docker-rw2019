#!/usr/bin/env python
with open('template.html') as template, open('README.md') as readme, open('presentation.html', 'w') as pres:
    pres.write(
        template.read()
        .replace('PLACEHOLDERPLACEHOLDERPLACEHOLDER', readme.read())
    )
