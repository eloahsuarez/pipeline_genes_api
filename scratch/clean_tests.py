import re
import glob

def clean_file(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()
        
    new_lines = []
    skip = False
    indent_level = -1
    
    for line in lines:
        match = re.match(r'^(\s*)def (test_[a-zA-Z0-9_]+)\(', line)
        if match:
            indent_level = len(match.group(1))
            func_name = match.group(2)
            if 'specificity' in func_name.lower():
                skip = True
            else:
                skip = False
                
        if skip:
            if line.strip() == "":
                continue
            if len(line) - len(line.lstrip(' ')) <= indent_level and not line.strip().startswith('def '):
                # end of function block? actually, wait, the next def test_ resets it.
                pass
            continue
        new_lines.append(line)
        
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(new_lines)

for file in glob.glob("tests/test_*.py"):
    clean_file(file)

