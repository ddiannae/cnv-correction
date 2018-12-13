import requests
import json
import time

fields = [
    "submitter_id",
    "case_id",
    "aliquot_ids"
    ]

fields = ",".join(fields)

cases_endpt = "https://api.gdc.cancer.gov/cases"
aliquot_id = "1dd196c1-6198-4edf-b8a1-66fed90a68b3"
all_caseid_submitterid_aliquotid = ""

f = open("cancer-aliquotes.tsv", "r")
for x in f:
    aliquot_id = x.strip()
    filters = {
        "op": "=",
        "content":{
            "field": "aliquot_ids",
            "value": [aliquot_id]
            }
        }
    print "Request for: " + aliquot_id + "\n"

    # With a GET request, the filters parameter needs to be converted
    # from a dictionary to JSON-formatted string
    params = {
        "filters": json.dumps(filters),
        "fields": fields,
        "format": "JSON",
        "size": "100"
        }
    response = requests.get(cases_endpt, params = params)
    content = json.loads(response.content)
    if len(content["data"]["hits"]) > 0:
        first_hit = content["data"]["hits"][0]
        matching_aliquot_id = next(x for x in first_hit["aliquot_ids"] if x == aliquot_id)
        caseid_submitterid_aliquotid = first_hit["case_id"] + "\t" + first_hit["submitter_id"] + "\t" + matching_aliquot_id + "\n"
        all_caseid_submitterid_aliquotid = all_caseid_submitterid_aliquotid + caseid_submitterid_aliquotid
    else:
        print "No matching entries for aliquote: " + aliquot_id
print(all_caseid_submitterid_aliquotid)
f = open("cancer-caseid-aliquotes.tsv", "w")
f.write(all_caseid_submitterid_aliquotid)
