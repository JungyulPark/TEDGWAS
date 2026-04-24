import docx

out = open(r'c:\ProjectTEDGWAS\docs_dump.txt', 'w', encoding='utf-8')

for fpath in [r'c:\ProjectTEDGWAS\TED_TRAP_Main_Tables_1_2.docx', r'c:\ProjectTEDGWAS\TED_TRAP_Tables_1_2_v2.docx']:
    out.write(f"\n================= {fpath} =================\n")
    try:
        doc = docx.Document(fpath)
        for t_idx, table in enumerate(doc.tables):
            out.write(f"--- Table {t_idx + 1} ---\n")
            for row in table.rows:
                cells = [c.text.strip().replace('\n', ' ') for c in row.cells]
                out.write(" | ".join(cells) + "\n")
    except Exception as e:
        out.write(f"ERROR: {e}\n")

out.close()
