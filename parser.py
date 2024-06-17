import argparse
import xml.etree.ElementTree as ET

import pandas as pd
from tqdm import tqdm


def load_xml(xml_file):
    tree = ET.parse(xml_file)
    root = tree.getroot()
    return root


def parse_drug(root, ns):
    res = pd.DataFrame()
    for drug in tqdm(root.findall("db:drug", ns)):
        drug_name = drug.find("db:name", ns).text
        drugbank_id = drug.find("db:drugbank-id", ns).text[2:]
        res = parse_targets(drug, ns, res, drug_name, drugbank_id)
    return res


def parse_targets(drug, ns, res, drug_name, drugbank_id):
    properties = {
        "SMILES": get_property_value(drug, ns, "SMILES"),
        "PubChem Compound": get_property_value(drug, ns, "PubChem Compound"),
        "PubChem Substance": get_property_value(drug, ns, "PubChem Substance"),
    }
    targets = drug.find("db:targets", ns)
    if targets is not None:
        for target in targets.findall("db:target", ns):
            target_info = {
                "ID": get_target_id(target, ns),
                "Name": get_target_name(target, ns),
                "UniProt ID": get_uniprot_id(target, ns),
                "Gene Name": get_gene_name(target, ns),
            }
            res = parse_articles(
                target, ns, res, drug_name, drugbank_id, **properties, **target_info
            )
    return res


def get_property_value(element, ns, property_name):
    return (
        next(
            (
                p.find("db:value", ns).text
                for p in element.findall(
                    f'db:calculated-properties/db:property[db:kind="{property_name}"]',
                    ns,
                )
            ),
            None,
        )
        or ""
    )


def get_target_id(target, ns):
    return (
        target.find("db:id", ns).text if target.find("db:id", ns) is not None else "N/A"
    )


def get_target_name(target, ns):
    return (
        target.find("db:name", ns).text
        if target.find("db:name", ns) is not None
        else "N/A"
    )


def get_uniprot_id(target, ns):
    polypeptide = target.find("db:polypeptide", ns)
    return polypeptide.get("id") if polypeptide is not None else "N/A"


def get_gene_name(target, ns):
    polypeptide = target.find("db:polypeptide", ns)
    if polypeptide is None:
        return "N/A"
    gene_name = polypeptide.find("db:gene-name", ns)
    return gene_name.text if gene_name is not None else "N/A"


def parse_articles(target, ns, res, drug_name, drugbank_id, **kwargs):
    articles = target.findall("db:references/db:articles/db:article", ns)
    pmids = [
        article.find("db:pubmed-id", ns).text
        for article in articles
        if article.find("db:pubmed-id", ns) is not None
        and article.find("db:pubmed-id", ns).text is not None
    ]
    pubmed_ids = ";".join(pmids) if pmids else "N/A"
    res = pd.concat(
        [
            res,
            pd.DataFrame(
                {
                    "Drug": [drug_name],
                    "DrugBank ID": [drugbank_id],
                    "SMILES": [kwargs["SMILES"]],
                    "PubChem CID": [kwargs["PubChem Compound"]],
                    "PubChem SID": [kwargs["PubChem Substance"]],
                    "Target ID": [kwargs["ID"]],
                    "Target Name": [kwargs["Name"]],
                    "UniProt ID": [kwargs["UniProt ID"]],
                    "Gene Name": [kwargs["Gene Name"]],
                    "PMIDs": [pubmed_ids],
                    "PMID Count": [len(pmids)],
                }
            ),
        ],
        ignore_index=True,
    )
    return res


def save_to_csv(res, output_file):
    res.to_csv(output_file, index=False)


if __name__ == "__main__":
    print("Starting to parse DrugBank XML file...")
    parser = argparse.ArgumentParser(
        description="Parse DrugBank XML file and save target information to CSV"
    )
    parser.add_argument(
        "--input", default="drugbank.xml", help="Path to the DrugBank XML file to parse"
    )
    parser.add_argument(
        "--output", default="targets.csv", help="Path to the output CSV file"
    )
    args = parser.parse_args()
    root = load_xml(args.input)
    ns = {"db": "http://www.drugbank.ca"}
    res = parse_drug(root, ns)
    save_to_csv(res, args.output)
    print("Parsing completed. Results saved to", args.output)
