import re

with open("app.py", "r", encoding="utf-8") as f:
    lines = f.readlines()

new_lines = []
for line in lines:
    if "self._specificity_service().clear_cache()" in line:
        continue
    if "self._clear_specificity_results()" in line:
        continue
    if 'self.vars["ncbi_specificity_top_n"].get()' in line:
        new_lines.append(line.replace('self.vars["ncbi_specificity_top_n"].get()', '10')) # fallback value
        continue
    if "self.ncbi_specificity_status.set(" in line:
        continue
    new_lines.append(line)

with open("app.py", "w", encoding="utf-8") as f:
    f.writelines(new_lines)
