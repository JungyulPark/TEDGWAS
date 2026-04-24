import pandas as pd

excel_path = 'c:/ProjectTEDGWAS/Final_Numbers_Frozen_20260421.xlsx'

print('=== A. Gene_Sets sheet ===')
gene_sets_df = pd.read_excel(excel_path, sheet_name='Gene_Sets')
print(gene_sets_df)

print('\n=== B. Intersections sheet ===')
intersections_df = pd.read_excel(excel_path, sheet_name='Intersections')
print(intersections_df)
