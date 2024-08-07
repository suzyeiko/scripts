import xml.etree.ElementTree as ET
import sys
import requests
import time
import csv

SEQ_METADATA_FILE = str(sys.argv[1])

#ACCESSION = str(sys.argv[1])
#SEQ_ACCESSION = str(sys.argv[2])

def parse_assembly_statistics(xml_data):
    root = ET.fromstring(xml_data)
    metadata = {
        "accession": "",
        "assembly_title": "",
        "assembly_level": "",
        "genome_representation": "",
        "taxon_id": "",
        "scientific_name": "",
        "strain": "",
        "sample_id": "",
        "study_id": "",
        "total_length": "",
        "ungapped_length": "",
        "number_of_contigs": "",
        "number_of_scaffolds": "",
        "number_of_replicon": "",
        "number_of_non_chromosome_replicons": "",
        "contig_n50": "",
        "contig_l50": "",
        "contig_n75": "",
        "contig_n90": "",
        "scaffold_n50": "",
        "scaffold_l50": "",
        "scaffold_n75": "",
        "scaffold_n90": "",
        "chromosomes": {},
    }
    # Extract the assembly accession
    identifiers = root.find(".//IDENTIFIERS")
    if identifiers is not None:
        metadata["accession"] = identifiers.find("PRIMARY_ID").text
    # Extract single-value attributes
    for element, field in [
        (".//TITLE", "assembly_title"),
        (".//ASSEMBLY_LEVEL", "assembly_level"),
        (".//GENOME_REPRESENTATION", "genome_representation"),
    ]:
        assembly_info = root.find(element)
        if assembly_info is not None:
            metadata[field] = assembly_info.text
    # Extract taxon information
    taxon = root.find(".//TAXON")
    if taxon is not None:
        metadata["taxon_id"] = taxon.find("TAXON_ID").text
        metadata["scientific_name"] = taxon.find("SCIENTIFIC_NAME").text
        metadata["strain"] = taxon.find("STRAIN").text
    # Extract sample and study IDs
    sample_ref = root.find(".//SAMPLE_REF")
    if sample_ref is not None:
        metadata["sample_id"] = sample_ref.find("IDENTIFIERS/PRIMARY_ID").text
    study_ref = root.find(".//STUDY_REF")
    if study_ref is not None:
        metadata["study_id"] = study_ref.find("IDENTIFIERS/PRIMARY_ID").text
    # Extract assembly attributes
    for attribute in root.findall(".//ASSEMBLY_ATTRIBUTE"):
        tag = attribute.find("TAG").text
        value = attribute.find("VALUE").text
        mapping = {
            "total-length": "total_length",
            "ungapped-length": "ungapped_length",
            "count-contig": "number_of_contigs",
            "scaffold-count": "number_of_scaffolds",
            "replicon-count": "number_of_replicon",
            "count-non-chromosome-replicon": "number_of_non_chromosome_replicons",
            "contig-n50": "contig_n50",
            "contig-L50": "contig_l50",
            "contig-n75": "contig_n75",
            "contig-n90": "contig_n90",
            "n50": "scaffold_n50",
            "scaf-L50": "scaffold_l50",
            "scaf-n75": "scaffold_n75",
            "scaf-n90": "scaffold_n90",
        }
        if tag in mapping:
            metadata[mapping[tag]] = value
    # Extract chromosome information
    chromosomes = root.findall(".//CHROMOSOME")
    metadata["chromosomes"] = {
        chromosome.get("accession"): (
            chromosome.find("NAME").text if chromosome.find("NAME") is not None else "",
            chromosome.find("TYPE").text if chromosome.find("TYPE") is not None else "",
        )
        for chromosome in chromosomes
    }
    return metadata

def get_assembly_metadata(accession, n_tries=3):
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{accession}"
    response = requests.get(url)
    if response.status_code == 200:
        return parse_assembly_statistics(response.text)
    elif n_tries > 0:
        return get_assembly_metadata(accession, n_tries - 1)
    else:
        response.raise_for_status()

def tabulate_output(metadata,SEQ_ACCESSION):
    output = SEQ_ACCESSION+"\t"
    for key, value in metadata.items():
        if key == "sequence_version":
            output = output+value
        elif key == "chromosomes":
            continue
        else:
            output = output+value+"\t"
    return output

def test_script(accession="GCA_040084815", seq_accession="CP157586"):
    assembly_metadata = get_assembly_metadata(accession)
    print(assembly_metadata)
    #print(result)
    #print(assembly_metadata.keys())
    #print(f"Accession: {assembly_metadata['accession']}")
    #print(f"Total length: {assembly_metadata['total_length']}")
    #print(f"Assembly level: {assembly_metadata['assembly_level']}")
    #print(f"Number of scaffolds: {assembly_metadata['number_of_scaffolds']}")
    #print(f"Scaffold N50: {assembly_metadata['scaffold_n50']}")
    #print(f"Scaffold L50: {assembly_metadata['scaffold_l50']}")

#test_script()


with open(SEQ_METADATA_FILE) as file:
    tsv_file = csv.reader(file, delimiter="\t")
    for line in tsv_file:
        try:
            ASSEMBLY_ACCESSION = line[55]
            SEQ_ACCESSION = line[0]
            #print(SEQ_ACCESSION, ASSEMBLY_ACCESSION)
            if SEQ_ACCESSION != "accession":
                assembly_metadata = get_assembly_metadata(ASSEMBLY_ACCESSION)
                result = tabulate_output(assembly_metadata, SEQ_ACCESSION)
                print(result)
            elif SEQ_ACCESSION == "accession":
                header = "sequence_accession\tassembly_accession\tassembly_title\tassembly_level\tgenome_representation\ttaxon_id\tscientific_name\tstrain\tsample_id\tstudy_id\ttotal_length\tungapped_length\tnumber_of_contigs\tnumber_of_scaffolds\tnumber_of_replicon\tnumber_of_non_chromosome_replicons\tcontig_n50\tcontig_l50\tcontig_n75\tcontig_n90\tscaffold_n50\tscaffold_l50\tscaffold_n75\tscaffold_n90"
                print(header)
            else:
                continue
        except:
            print(SEQ_ACCESSION+"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t")
            #print(SEQ_ACCESSION)
 
