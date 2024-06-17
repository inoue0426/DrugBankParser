import argparse
import xml.etree.ElementTree as ET

import pandas as pd
from tqdm import tqdm


def parse_xml(xml_file_path):
    print("Parsing XML file...")
    tree = ET.parse(xml_file_path)
    root = tree.getroot()
    namespace = "{http://www.drugbank.ca}"
    drug_data_list = []

    for drug in tqdm(root.findall(f"{namespace}drug")):
        drug_info = {
            "Drug Name": (
                drug.find(f"{namespace}name").text
                if drug.find(f"{namespace}name") is not None
                else "N/A"
            ),
            "DrugBank ID": None,
            "PubChem CID": None,
            "PubChem SID": None,
            "SMILES": None,
            "Targets Name": [],
            "Targets": [],
        }
        _parse_drug_info(drug, namespace, drug_info)
        drug_data_list.append(drug_info)

    return drug_data_list


def _parse_drug_info(drug, namespace, drug_info):
    print("Parsing drug information...")
    _parse_drugbank_id(drug, namespace, drug_info)
    _parse_external_identifiers(drug, namespace, drug_info)
    _parse_calculated_properties(drug, namespace, drug_info)
    _parse_targets(drug, namespace, drug_info)


def _parse_drugbank_id(drug, namespace, drug_info):
    print("Parsing DrugBank ID...")
    drug_id_element = drug.find(f"{namespace}drugbank-id[@primary='true']")
    if drug_id_element is not None:
        drug_info["DrugBank ID"] = drug_id_element.text


def _parse_external_identifiers(drug, namespace, drug_info):
    print("Parsing external identifiers...")
    external_identifiers = drug.find(f"{namespace}external-identifiers")
    if external_identifiers is not None:
        for identifier in external_identifiers:
            resource = identifier.find(f"{namespace}resource")
            id_value = identifier.find(f"{namespace}identifier")
            if resource is not None and id_value is not None:
                if resource.text == "PubChem Compound":
                    drug_info["PubChem CID"] = id_value.text
                elif resource.text == "PubChem Substance":
                    drug_info["PubChem SID"] = id_value.text


def _parse_calculated_properties(drug, namespace, drug_info):
    calculated_properties = drug.find(f"{namespace}calculated-properties")
    if calculated_properties is not None:
        for prop in calculated_properties:
            kind = prop.find(f"{namespace}kind")
            value = prop.find(f"{namespace}value")
            if kind is not None and value is not None and kind.text == "SMILES":
                drug_info["SMILES"] = value.text


def _parse_targets(drug, namespace, drug_info):
    targets = drug.find(f"{namespace}targets")
    if targets is not None:
        for target in targets:
            target_name = target.find(f"{namespace}name")
            if target_name is not None:
                drug_info["Targets Name"].append(target_name.text)
            polypeptide = target.find(f"{namespace}polypeptide")
            if polypeptide is not None and "id" in polypeptide.attrib:
                drug_info["Targets"].append(polypeptide.attrib["id"])


def main():
    parser = argparse.ArgumentParser(description="DrugBank XML Parser")
    parser.add_argument("xml_file", help="Path to the DrugBank XML file")
    parser.add_argument("output_file", help="Path to the output CSV file")
    args = parser.parse_args()

    drug_data_list = parse_xml(args.xml_file)

    drug_df = pd.DataFrame(drug_data_list)
    drug_df.to_csv(args.output_file, index=False)


if __name__ == "__main__":
    main()
