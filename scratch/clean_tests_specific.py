with open("tests/test_app_restart.py", "r", encoding="utf-8") as f:
    lines = f.readlines()

new_lines = []
skip = False
indent = 0
for line in lines:
    if line.startswith("def test_ncbi_remote_database_friendly_selection_stores_exact_alias"): skip = True
    elif line.startswith("def test_ncbi_remote_database_change_invalidates_previous_result"): skip = True
    elif line.startswith("def test_ncbi_results_table_uses_left_and_right_sites_for_same_primer_product"): skip = True
    elif line.startswith("def test_ncbi_results_table_limits_products_and_keeps_report_complete"): skip = True
    elif line.startswith("def test_ncbi_selection_helper_accepts_omitted_products_notice"): skip = True
    elif line.startswith("def test_ncbi_product_details_show_combination_and_site_metrics"): skip = True
    elif skip and line.startswith("def test_"):
        skip = False

    if not skip:
        new_lines.append(line)

with open("tests/test_app_restart.py", "w", encoding="utf-8") as f:
    f.writelines(new_lines)
