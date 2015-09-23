import urllib
# If you are using Python 3+, import urllib instead of urllib2
import json


data =  {

        "Inputs": {

                "input1":
                {
                    "ColumnNames": ["sequence", "cutsite", "percentpeptide"],
                    "Values": [ [ "ATGAAATCTAACAATGCGCTCATCGTCATCCTCGGCACCGTCACCCTGGATGCTGTAGGCATAGGCTTGGTTATGCCGGTACTGCCGGGCCTCTTGCGGGATAT",
                    "10", "100" ], ]
                },        },
            "GlobalParameters": {
}
    }

body = str.encode(json.dumps(data))

url = 'https://ussouthcentral.services.azureml.net/workspaces/ee5485c1d9814b8d8c647a89db12d4df/services/c24d128abfaf4832abf1e7ef45db4b54/execute?api-version=2.0&details=true'
api_key = '4dLjIIiPRWRYXwW4ZTwgCXYWU78H9jF1W+F4T2ci7Nk7kpYjlvKw/wCfOHZwQAX/uMugAl/o7T2O0sbU7xSrtQ==' # Replace this with the API key for the web service
headers = {'Content-Type':'application/json', 'Authorization':('Bearer '+ api_key)}

req = urllib.request.Request(url, body, headers)


try:
    response = urllib.request.urlopen(req)

    # If you are using Python 3+, replace urllib2 with urllib.request in the above code:
    # req = urllib.request.Request(url, body, headers)
    # response = urllib.request.urlopen(req)

    result = response.read()
    print(result)
except urllib.error.HTTPError as error:
    print("The request failed with status code: " + str(error.code))

    # Print the headers - they include the requert ID and the timestamp, which are useful for debugging the failure
    print(error.info())

    print(json.loads(error.read().decode()))
