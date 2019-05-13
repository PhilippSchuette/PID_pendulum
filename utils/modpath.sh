#!/bin/bash
# replace references to style sheets and javascript sources
replacements=(alabaster.css basic.css custom.css pygments.css jquery.js underscore.js doctools.js language_data.js documentation_options.js searchindex.js)
html_files=(api/static/docs/index.html api/static/docs/search.html api/static/docs/genindex.html api/static/docs/py-modindex.html)

for file in "${html_files[@]}"; do
    for repl in "${replacements[@]}"; do
        sed -i "s%_static/$repl%static/docs/_static/$repl%g" "$file"
    done
done

# replace routes to html files
replacements=(index.html search.html genindex.html py-modindex.html)

for file in "${html_files[@]}"; do
    for repl in "${replacements[@]}"; do
        sed -i "s%href=\"$repl\"%href=\"/docs/$repl\"%g" "$file"
    done
done
