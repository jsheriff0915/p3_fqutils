import os
import sys
import re
import requests
import urllib.request, urllib.parse, urllib.error
import json
from requests.packages.urllib3.exceptions import InsecureRequestWarning
requests.packages.urllib3.disable_warnings(InsecureRequestWarning)

Debug = False #shared across functions defined here
LOG = sys.stderr
Base_url="https://www.bv-brc.org/api/"
PatricUser = None

def createTSVGet(api_url=None):
    if api_url == None:
        api_url=Base_url
    Session = requests.Session()
    Session.headers.update({ 'accept': "text/tsv" })
    Session.headers.update({ "Content-Type": "application/rqlquery+x-www-form-urlencoded" })
    if not authenticateByEnv(Session):
        authenticateByFile(None, Session)
    return Session


def getHostManifest(api_url=None):
    if api_url is None:
        api_url = Base_url
    r = requests.get(api_url + "content/host/patric_host_summary.json")
    manifest = r.json()
    return {item["species_taxid"]: item for item in manifest["genomes"]}


def authenticateByFile(tokenFile=None, Session=None):
    if not tokenFile:
        tokenFile = os.path.join(os.environ.get('HOME'), ".patric_token")
    if os.path.exists(tokenFile):
        LOG.write("reading auth key from file %s\n"%tokenFile)
        with open(tokenFile) as F:
            tokenString = F.read().rstrip()
            authenticateByString(tokenString, Session)
        return True
    return False


def authenticateByEnv(Session):
    if "KB_AUTH_TOKEN" in os.environ:
        LOG.write("reading auth key from environment\n")
        authenticateByString(os.environ.get('KB_AUTH_TOKEN'), Session)
        return True
    else:
        return authenticateByFile(None, Session)

def authenticateByString(tokenString, Session):
    Session.headers.update({ 'Authorization' : tokenString })
    if "Authorization" in Session.headers:
        global PatricUser
        PatricUser = Session.headers["Authorization"].split(r"|")[3].split("=")[1]
        LOG.write("Patric user = %s\n"%PatricUser)

def getGenomeIdsNamesByName(name, limit='10', Session=None):
    query = "eq(genome_name,%s)"%name
    query += "&select(genome_id,genome_name)"
    query += "&limit(%s)"%limit
    ret = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write(ret.url+"\n")
    return(ret.text.replace('"', ''))

def getGenomeGroupIds(genomeGroupName, Session):
    LOG.write("getGenomeGroupIds(%s), PatricUser=%s\n"%(genomeGroupName, PatricUser))
    genomeGroupSpecifier = PatricUser+"/home/Genome Groups/"+genomeGroupName
    genomeGroupSpecifier = "/"+urllib.parse.quote(genomeGroupSpecifier)
    genomeGroupSpecifier = genomeGroupSpecifier.replace("/", "%2f")
    query = "in(genome_id,GenomeGroup("+genomeGroupSpecifier+"))"
    query += "&select(genome_id)"
    query += "&limit(10000)"
    if Debug:
        LOG.write("requesting group %s for user %s\n"%(genomeGroupName, PatricUser))
        LOG.write("query =  %s\n"%(query))
    ret = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write(ret.url+"\n")
    return(ret.text.replace('"', '').split("\n"))[1:-1]

def getNamesForGenomeIds(genomeIds, Session):
#    return getDataForGenomes(genomeIdSet, ["genome_id", "genome_name"])
    retval = {}
    for genome in genomeIds:
        retval[genome] = ""
    query="in(genome_id,("+",".join(genomeIds)+"))&select(genome_id,genome_name)"
    response = Session.get(Base_url+"genome/", params=query) #, 
    if Debug:
        LOG.write("    response URL: %s\n"%response.url)
        LOG.write("    len(response.text)= %d\n"%len(response.text))
    if not response.ok:
        LOG.write("Error code %d returned by %s in getNamesForGenomeIds\n"%(response.status_code, response.url))
    for line in response.text.split("\n"):
        line = line.replace('"','')
        row = line.split("\t", 1)
        if len(row) >= 2:
            genome, name = row
            retval[genome] = name
    return retval

def getNamesForGenomeIdsByN(genomeIds, n=5):
    """ For some reason, grabbing them in bulk misses some, so grab N at a time.
    """
    retval = {}
    i = 0
    genomeIds = list(genomeIds)
    while i < len(genomeIds):
        subset = genomeIds[i:i+n]
        retval.update(getNamesForGenomeIds(subset))
        i += n
    return retval


def getGenomeIdByFieldValue(queryField, queryValue, Session):
    query = "eq(%s,%s)"%(queryField, queryValue)
    query += "&select(genome_id)"
    req = Session.get(Base_url+"genome/", params=query) 
    if Debug:
        LOG.write("getGenomeIdsByQuery: "+req.url+"\n")
        LOG.write(req.text+"\n")
    data = req.text.split("\n")
    genomeId = ""
    if len(data) > 1:
        genomeId = data[1]
        genomeId = genomeId.replace('\"', '')
    return genomeId

def getDataForGenomes(genomeIdSet, fieldNames, Session):
    query = "in(genome_id,(%s))"%",".join(genomeIdSet)
    if fieldNames:
        query += "&select(%s)"%",".join(fieldNames)
    query += "&limit(%s)"%len(genomeIdSet)

    response = Session.get(Base_url+"genome/", params=query)
    if Debug:
        LOG.write("getDataForGenomes:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getDataForGenomes\n"%(response.status_code, response.url))
        LOG.write("length of query was %d\n"%len(query))
        LOG.write("url="+req.url+"\nquery="+query+"\n")
        raise Exception(errorMessage)
    data = response.text.replace('"','') #get rid of quotes
    rows = data.split("\n")[:-1] # leave off empty last element
    retval = []
    for row in rows:
        fields = row.split("\t")
        #if len(fields) != len(fieldNames):
         #   continue
        retval.append(fields)
    return(retval)

def getProteinFastaForPatricIds(patricIds, Session):
    query="in(patric_id,("+",".join(map(urllib.parse.quote, patricIds))+"))"
    query += "&limit(%d)"%len(patricIds)
    response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/protein+fasta'})
    if False and Debug:
        LOG.write("getProteinFastaForByPatricIds:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getProteinFastaForPatricIds\n"%(response.status_code, Base_url))
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        LOG.write(errorMessage)
        LOG.flush()
        raise Exception(errorMessage)
    idsFixedFasta=""
    for line in response.text.split("\n"):
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) > 2:
                line = "|".join(parts[:2])
        idsFixedFasta += line+"\n"
    return idsFixedFasta
    
def getDnaFastaForPatricIds(patricIds, Session):
    query="in(patric_id,("+",".join(map(urllib.parse.quote, patricIds))+"))"
    query += "&limit(%d)"%len(patricIds)
    response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/dna+fasta'})
    if False and Debug:
        LOG.write("getDnaFastaForByPatricIds:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getDnaFastaForPatricIds\n"%(response.status_code, Base_url))
        errorMessage= "Error code %d returned by %s in getGenomeFeaturesByPatricIds\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        LOG.write(errorMessage)
        LOG.flush()
        raise Exception(errorMessage)
    idsFixedFasta=""
    for line in response.text.split("\n"):
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) > 2:
                line = "|".join(parts[:2])
        idsFixedFasta += line+"\n"
    return idsFixedFasta
    
def getProteinsFastaForGenomeId(genomeId, Session):
    query="in(genome_id,("+genomeId+"))"
    query += "&limit(25000)"
    response=Session.get(Base_url+"genome_feature/", params=query, headers={'Accept': 'application/protein+fasta'})
    if Debug:
        LOG.write("getProteinsFastaForGenomeId:\nurl="+response.url+"\nquery="+query+"\n")
    if not response.ok:
        LOG.write("Error code %d returned by %s in getProteinsFastaForGenomeId\n"%(response.status_code, Base_url))
        errorMessage= "Error code %d returned by %s in getProteinsFastaForGenomeId\nlength of query was %d\n"%(response.status_code, Base_url, len(query))
        LOG.write(errorMessage)
        LOG.flush()
        raise Exception(errorMessage)
    idsFixedFasta=""
    for line in response.text.split("\n"):
        if line.startswith(">"):
            parts = line.split("|")
            if len(parts) > 2:
                line = "|".join(parts[:2])+"\n"
        idsFixedFasta += line
    return idsFixedFasta

def getProductsForPgfams(pgfams, Session):
    retval = {}
    for pgfam in pgfams:
        retval[pgfam] = ""
    query="in(family_id,("+",".join(pgfams)+"))&select(family_id,family_product)"
    response = Session.get(Base_url+"protein_family_ref/", params=query) #, 
    if Debug:
        LOG.write("    response URL: %s\n"%response.url)
        LOG.write("    len(response.text)= %d\n"%len(response.text))
    if not response.ok:
        LOG.write("Error code %d returned by %s in getProductsForPgfams\n"%(response.status_code, response.url))
    for line in response.text.split("\n"):
        line = line.replace('"','')
        row = line.split("\t", 1)
        if len(row) >= 2:
            pgfam, product = row
            retval[pgfam] = product
    return retval

def getProductsForPgfamsByN(pgfams, n=5, Session=None):
    """ For some reason, grabbing them in bulk misses some, so grab N at a time.
    """
    retval = {}
    i = 0
    pgfams = list(pgfams)
    while i < len(pgfams):
        subset = pgfams[i:i+n]
        retval.update(getProductsForPgfams(subset))
        i += n
    return retval

def getPatricGenesPgfamsForGenomeSet(genomeIdSet, Session):
    if Debug:
        LOG.write("getPatricGenesPgfamsForGenomeSet() called for %d genomes\n"%len(genomeIdSet))
        LOG.write("    Session headers=\n"+str(Session.headers)+"\n")
    retval = []
    # one genome at a time, so using 'get' should be fine
    for genomeId in genomeIdSet:
        query="and(%s,%s,%s)"%("eq(genome_id,(%s))"%genomeId, "eq(feature_type,CDS)", "eq(pgfam_id,PGF*)")
        query += "&select(genome_id,patric_id,pgfam_id)"
        query += "&limit(25000)"
        response = Session.get(Base_url+"genome_feature/", params=query) #, 
        """
        req = requests.Request('POST', Base_url+"genome_feature/", data=query)
        prepared = Session.prepare_request(req) #req.prepare()
        response=Session.send(prepared, verify=False)
        """
        if Debug:
            LOG.write("    response URL: %s\n"%response.url)
            LOG.write("    len(response.text)= %d\n"%len(response.text))
        curLen = len(retval)
        for line in response.text.split("\n"):
            line = line.replace('"','')
            row = line.split("\t")
            if len(row) != 3:
                continue
            if not row[2].startswith("PGF"):
                continue
            retval.append(row)
        if Debug:
            LOG.write("    got %d pgfams for that genome\n"%(len(retval)-curLen))
    return(retval)

def getPgfamGenomeMatrix(genomeIdSet, ggpMat = None):
    """ Given list of genome ids: 
        tabulate genes per genome per pgfam 
        (formats data from getPatricGenesPgfamsForGenomeSet as table)
    """
    genomeGenePgfamList = getPatricGenesPgfamsForGenomeSet(genomeIdSet)
    if not ggpMat: # if a real value was passed, extend it
        ggpMat = {} # genome-gene-pgfam matrix (really just a dictionary)
    for row in genomeGenePgfamList:
        genome, gene, pgfam = row
        if pgfam not in ggpMat:
            ggpMat[pgfam] = {}
        if genome not in ggpMat[pgfam]:
            ggpMat[pgfam][genome] = []
        ggpMat[pgfam][genome].append(gene)
    return ggpMat

def getPgfamCountMatrix(genomeIdSet, ggpMat = None):
    """ Given list of genome ids: 
        tabulate counts per genome per pgfam 
        (formats data from getPatricGenesPgfamsForGenomeSet as table)
    """
    genomeGenePgfamList = getPatricGenesPgfamsForGenomeSet(genomeIdSet)
    if not ggpMat: # if a real value was passed, extend it
        ggpMat = {} # genome-gene-pgfam matrix (really just a dictionary)
    for row in genomeGenePgfamList:
        genome, gene, pgfam = row
        if pgfam not in ggpMat:
            ggpMat[pgfam] = {}
        if genome not in ggpMat[pgfam]:
            ggpMat[pgfam][genome] = 0
        ggpMat[pgfam][genome] += 1
    return ggpMat

def writePgfamGeneMatrix(ggpMat, fileHandle):
    """ write out pgfamGeneMatrix to file handle 
    data is list of genes per pgfam per genome
    rows are pgfams
    cols are genomes
    column headers identify genomes
    genes are comma-separated
    """
    # first collect set of all genomes
    genomeSet = set()
    for pgfam in ggpMat:
        genomeSet.update(set(ggpMat[pgfam].keys()))
    genomes = sorted(genomeSet)
    fileHandle.write("PGFam\t"+"\t".join(genomes)+"\n")
    for pgfam in ggpMat:
        fileHandle.write(pgfam)
        for genome in genomes:
            gene = ""
            if genome in ggpMat[pgfam]:
                gene = ",".join(ggpMat[pgfam][genome])
            fileHandle.write("\t"+gene)
        fileHandle.write("\n")

def writePgfamCountMatrix(ggpMat, fileHandle):
    """ write out matrix of counts per pgfam per genome to file handle 
    data is count of genes per pgfam per genome (integers)
    rows are pgfams
    cols are genomes
    column headers identify genomes
    """
    # first collect set of all genomes
    genomeSet = set()
    for pgfam in ggpMat:
        genomeSet.update(set(ggpMat[pgfam].keys()))
    genomes = sorted(genomeSet)
    fileHandle.write("PGFam\t"+"\t".join(genomes)+"\n")
    for pgfam in ggpMat:
        fileHandle.write(pgfam)
        for genome in genomes:
            count = 0
            if genome in ggpMat[pgfam]:
                count = len(ggpMat[pgfam][genome])
            fileHandle.write("\t%d"%count)
        fileHandle.write("\n")

def readPgfamGeneMatrix(fileHandle):
    """ read pgfamGeneMatrix from file handle
    Data are list of genes (comma-delimited) per genome per pgfam
    rows are pgfams, cols are genomes, column headers identify genomes
    """
    # genome ids are headers in first line
    header = fileHandle.readline().rstrip()
    genomes = header.split("\t")[1:] # first entry is placeholder for pgfam rownames
    pgMat = {} # genome-gene-pgfam matrix (really just a dictionary)
    for row in fileHandle:
        fields = row.rstrip().split("\t")
        pgfam = fields[0]
        pgMat[pgfam] = {}
        data = fields[1:]
        for i, genome in enumerate(genomes):
            if len(data[i]):
                pgMat[pgfam][genome] = data[i]
    return pgMat

def readPgfamCountMatrix(fileHandle):
    """ read pgfamCountMatrix from file handle
    rows are pgfams
    cols are genomes
    data are integer counts of that pgfam in that genome
    column headers identify genomes
    """
    # genome ids are headers in first line
    header = fileHandle.readline().rstrip()
    genomes = header.split("\t")[1:] # first entry is placeholder for pgfam rownames
    pcMat = {} # pgfam count matrix (really just a dictionary)
    for row in fileHandle:
        fields = row.rstrip().split("\t")
        pgfam = fields[0]
        pcMat[pgfam] = {}
        data = fields[1:]
        for i, genome in enumerate(genomes):
            pcMat[pgfam][genome] = int(float(data[i]))
    return pcMat

def getPatricGenesPgfamsForGenomeObject(genomeObject):
# parse a PATRIC genome object (read from json format) for PGFams
    retval = [] # a list of tupples of (genomeId, Pgfam, geneId)
    genomeId = genomeObject['id']
    for feature in genomeObject['features']:
        if 'family_assignments' in feature:
            for familyAssignment in feature['family_assignments']:
                if familyAssignment[0] == 'PGFAM':
                    retval.append((genomeId, feature['id'], familyAssignment[1]))
    return retval

def getGenomeObjectProteins(genomeObject):
# return dictionary of patricId -> BioPython.SeqRecord
    genomeId = genomeObject['id']
    retval = {}
    for feature in genomeObject['features']:
        patricId, product, genomeId, aa_sequence = '', '', '', ''
        patricId = feature['id']
        if "protein_translation" in feature:
            aa_sequence = feature["protein_translation"]
        if 'function' in feature:
            product = feature['function']
        simpleSeq = Seq(aa_sequence, IUPAC.extended_protein)
        seqRecord = SeqRecord(simpleSeq, id=patricId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval[patricId] = seqRecord
    return retval

def getGenomeObjectGeneDna(genomeObject):
# return dictionary of patricId -> BioPython.SeqRecord
    genomeId = genomeObject['id']
    contigSeq = {}
    for contig in genomeObject['contigs']:
        contigSeq[contig['id']] = contig['dna']
    retval = {} # dict of SeqRecords
    for feature in genomeObject['features']:
        geneId = feature['id']
        if geneId not in patricIds:
            continue
        product = ''
        if 'product' in feature:
            product = feature['function']
        if not 'location' in feature:
            continue
        contig, start, ori, length = feature['location'][0] # this should be an array of (contig, start, orientation, length)
        start = int(float(start))
        length = int(float(length))
        if ori == '+':
            start -= 1
            simpleSeq = Seq(contigSeq[contig][start:start+length], IUPAC.ambiguous_dna)
        if ori == '-':
            simpleSeq = Seq(contigSeq[contig][start-length:start], IUPAC.ambiguous_dna)
            simpleSeq = simpleSeq.reverse_complement()

        seqRecord = SeqRecord(simpleSeq, id=geneId, description=product)
        seqRecord.annotations["genome_id"] = genomeId
        retval[geneId] = seqRecord
    return retval

