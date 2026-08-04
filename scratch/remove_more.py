with open("app.py", "r", encoding="utf-8") as f:
    lines = f.readlines()

new_lines = []
skip = False
indent = 0
for line in lines:
    if line.startswith("    def show_ncbi_primer_fasta"): skip = True
    elif line.startswith("    def export_ncbi_primer_fasta"): skip = True
    elif line.startswith("    def _ncbi_primer_fasta"): skip = True
    elif skip and line.startswith("    def "): skip = False

    if not skip:
        new_lines.append(line)

with open("app.py", "w", encoding="utf-8") as f:
    f.writelines(new_lines)
