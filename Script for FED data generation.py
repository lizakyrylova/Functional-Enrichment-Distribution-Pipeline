import pandas as pd
from neo4j import GraphDatabase

#variables
filter_columns_output = r"PATH" #path to your functional annotation library (can be the one provided under filtered_combined_file.tsv)
gene_list = r"PATH" #your input gene list
script1_path = r"PATH" #enter your path for saving on your computer

uri = "bolt://localhost:XXXX" #enter number
username = #enter your username for neo4j
password = #enter your passowrd for neo4j
cypher_query_for_root_superpathway = """
MATCH (p:Pathway{{stId:{input_id}}})
OPTIONAL MATCH path = (p)<-[:hasEvent*]-(sp:Pathway)
WITH p, sp
ORDER BY LENGTH(path) DESC
LIMIT 1
RETURN
  p.stId AS Pathway,
  CASE
    WHEN sp IS NULL THEN p.stId
    ELSE COALESCE(sp.stId, p.stId)
  END AS SuperPathway,
  CASE
    WHEN sp IS NULL THEN p.displayName
    ELSE COALESCE(sp.displayName, p.displayName)
  END AS DisplayName
"""
#
cypher_query_for_all_root_superpathways = """
MATCH (p:Pathway{{stId:{input_id}}})
OPTIONAL MATCH path = (p)<-[:hasEvent*]-(sp:Pathway)
WITH p, sp
ORDER BY LENGTH(path) DESC
RETURN
  p.stId AS Pathway,
  CASE
    WHEN sp IS NULL THEN p.stId
    ELSE COALESCE(sp.stId, p.stId)
  END AS SuperPathway,
  CASE
    WHEN sp IS NULL THEN p.displayName
    ELSE COALESCE(sp.displayName, p.displayName)
  END AS DisplayName
"""
#to obtain the desired heirarchy level, edit "SKIP #", where level n = SKIP n-1, i.e. in case where root superpathway is level 1
cypher_query_for_sub_root_superpathway = """
MATCH (p:Pathway{{stId:{input_id}}})
OPTIONAL MATCH path = (p)<-[:hasEvent*]-(sp:Pathway)
WITH p, sp
ORDER BY LENGTH(path) DESC
SKIP 3
LIMIT 1
RETURN
  p.stId AS Pathway,
  CASE
    WHEN sp IS NULL THEN p.stId
    ELSE COALESCE(sp.stId, p.stId)
  END AS SuperPathway,
  CASE
    WHEN sp IS NULL THEN p.displayName
    ELSE COALESCE(sp.displayName, p.displayName)
  END AS DisplayName
"""
concat_final_file_output = r"PATH" #path to your output
concat_and_non_redundant_final_file_output = r"PATH" #path to your output (each gene is associated with non redundant functional annotations)
#load DFs
filter_columns_output_df = pd.read_csv(filter_columns_output, sep='\t')
gene_list_df = pd.read_csv(gene_list, sep='\t')

#Script 1: Find Reactome ID Indexes
def script1(filter_columns_output_df, gene_list_df):
    gene_identifier_to_rows = {}

    for idx, row in filter_columns_output_df.iterrows():
        gene_identifiers = row['genes_hugo'].split(';')
        for gene_id in gene_identifiers:
            if gene_id not in gene_identifier_to_rows:
                gene_identifier_to_rows[gene_id] = []
            gene_identifier_to_rows[gene_id].append(idx)

    result_data = {'Gene Identifier': [], 'Associated Rows': []}

    for idx, row in gene_list_df.iterrows():
        gene_id = row['Gene Identifier']
        associated_rows = gene_identifier_to_rows.get(gene_id, [])
        result_data['Gene Identifier'].append(gene_id)
        result_data['Associated Rows'].append(associated_rows)

    script1_result_df = pd.DataFrame(result_data)
    return script1_result_df

script1_result = script1(filter_columns_output_df, gene_list_df)
print(script1_result)
script1_result.to_csv(script1_path, sep='\t', index=False)

#Script 2: Change Indexes to associated REAC IDs
def script2(script1_path, filter_columns_output_df):
    script1_result_df = pd.read_csv(script1_path, sep='\t')
    row_index_to_REAC = dict(zip(filter_columns_output_df.index, filter_columns_output_df['term_id_REAC']))
    def replace_with_REAC(row_indexes):
        return [row_index_to_REAC[int(idx)] for idx in row_indexes if idx]
    script1_result_df['Associated Rows'] = script1_result_df['Associated Rows'].str.strip('[]').str.split(', ')
    script1_result_df['Associated Rows'] = script1_result_df['Associated Rows'].apply(replace_with_REAC)

    return script1_result_df

script2_result = script2(script1_path, filter_columns_output_df)
print(script2_result)
#script2_result.to_csv(r"PATH", sep='\t', index=False)

#Script 3: Remove genes without the Reactome IDs
def script3(script2_result):
    df = script2_result.copy()

    def process_associated_rows(row):
        values = row['Associated Rows']
        if not values:
            return pd.Series({'REAC ID': ''})
        cleaned_values = ["'R" + value.strip('REAC:') + "'" for value in values]
        result = ', '.join(cleaned_values)
        return pd.Series({'REAC ID': result})

    df[['REAC ID']] = df.apply(process_associated_rows, axis=1)
    df = df[df['REAC ID'] != '']
    df = df.drop(columns=['Associated Rows'])

    return df

script3_result = script3(script2_result)
print (script3_result)
#script3_result.to_csv(r"PATH", sep='\t', index=False)

#Script 4: Melting the Reactome ID dataframe
def script4(script3_result):
    script3_result_reset = script3_result.reset_index(drop=True)
    s = script3_result_reset['REAC ID'].str.split(',', expand=True).stack()
    s.index = s.index.droplevel(1)

    melted_df = pd.concat([script3_result_reset.drop('REAC ID', axis=1), s.rename('REAC ID')], axis=1)
    return melted_df

script4_result = script4(script3_result)
print(script4_result)
#script4_result.to_csv(r"PATH", sep='\t', index=False)

#Script 5: Get needed level of annotation resolution by querying Reactome
def script5(script4_result, uri, username, password):
    script4_result_df = pd.DataFrame(script4_result)
    input_ids = script4_result_df["REAC ID"].tolist()

    results = [{'Pathway': input_id, 'SuperPathway': None, 'DisplayName': None} for input_id in input_ids]

    cypher_query_template = cypher_query_for_root_superpathway #Choose Cypher query
    with GraphDatabase.driver(uri, auth=(username, password)) as driver:
        with driver.session() as session:
            # Iterate over the input_ids and update the results list with actual data if available
            for index, input_id in enumerate(input_ids):
                cypher_query = cypher_query_template.format(input_id=input_id)
                result = session.run(cypher_query).single()
                if result:
                    results[index].update(result.data())
                
    reactome_surf_output_df = pd.DataFrame(results)
    return reactome_surf_output_df

script5_result = script5(script4_result, uri, username, password)
print(script5_result)
#script5_result.to_csv(r"PATH", sep='\t', index=False)

#Script 6: Final Processing step
def script6(script5_result, script4_result):
    gene_identifier_column = script4_result[['Gene Identifier']]
    pathway_columns = script5_result[['Pathway', 'SuperPathway', 'DisplayName']]
    result_df = pd.concat([gene_identifier_column.reset_index(drop=True), pathway_columns.reset_index(drop=True)], axis=1)
    result_df = result_df.dropna(subset=['SuperPathway', 'DisplayName'], how='all')
    result_df_unique = result_df.drop_duplicates(subset=['Gene Identifier', 'SuperPathway'])
    return result_df, result_df_unique

script6_result, script6_result_unique = script6(script5_result, script4_result)
print (script6_result, script6_result_unique)
script6_result.to_csv(concat_final_file_output, sep='\t', index=False)
script6_result_unique.to_csv(concat_and_non_redundant_final_file_output, sep='\t', index=False)






    

