pandoc -s -S --bibliography paper.bib --filter pandoc-citeproc -f markdown -o paper.pdf paper.md
