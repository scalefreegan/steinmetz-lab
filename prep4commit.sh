cp .gitignore_template .gitignore
find . -size +25M | sed -r 's/^.{2}//' >> .gitignore
