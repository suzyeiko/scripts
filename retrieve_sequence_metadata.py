import sys
import requests
import time

ACCESSION = str(sys.argv[1])

def get_sequence_metadata(accession, n_tries=3):
    # Define the base URL for the ENA API
    base_url = "https://www.ebi.ac.uk/ena/portal/api/filereport"
    parameters = {
        "accession": accession,  # Sequence accession
        "result": "sequence",  # Specify the type of data you want
        "format": "json",  # Request the data in JSON format
        "fields": "all",  # Request all fields
    }
    # Make the GET request to the ENA API
    response = requests.get(base_url, params=parameters)
    # Check if the request was successful (status code 200)
    if response.status_code == 200:
        metadata = response.json()[0]
        return metadata
    elif n_tries > 0:
        # If the request was not successful, wait for 10 seconds and try again
        time.sleep(10)
        return get_sequence_metadata(accession, n_tries - 1)
    else:
        response.raise_for_status()
        print(accession)

def tabulate_output(metadata):
    output = ""
    for key, value in metadata.items():
        if key == "sequence_version":
            output = output+value
        else:
            output = output+value+"\t"
    return output

#def test_script(ACCESSION="CP059301"):
#    seq_metadata = get_sequence_metadata(ACCESSION)
#    print(seq_metadata)
#    return seq_metadata

#test_script("CP059301")

try:
    seq_metadata = get_sequence_metadata(ACCESSION)
    result = tabulate_output(seq_metadata)
    print(result)
except:
    print(ACCESSION)

#result = tabulate_output(seq_metadata)
#print(result)


