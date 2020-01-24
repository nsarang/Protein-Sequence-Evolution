# python version
from urllib import request, parse


SERVICELOCATION = "https://www.rcsb.org/pdb/rest/postBLAST/"

params = {
    "sequence": "VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTAVAHVDDMPNAL",
    "eCutOff": "10.0",
    "matrix": "BLOSUM62",
    "outputFormat": "XML",
}

data = parse.urlencode(params).encode()
req = request.Request(SERVICELOCATION, data=data)  # this will make the method "POST"
resp = request.urlopen(req)

print(resp.read().decode('utf-8'))
