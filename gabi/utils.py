################################################################################
################################# GABI UTILS ###################################
################################################################################
import sys
import requests
import aiohttp
import asyncio
import random
import numpy as np

# Default ensembl REST server URL
default_server='https://rest.ensembl.org'

def get_session(server=default_server):
    """return an aiohttp client session"""
    return aiohttp.ClientSession(base_url=server)

# Get an Ensembl region string from a region dict
get_region_str = lambda r: f"{r['seq_region']}:{r['start']}..{r['end']}:{r['strand']}"

def get_random_region(seq_region, sample_len, seq_region_length):
    """Create a random (uniform) sequence region within the range
    [1, seq_region_length] as a dict.

    NOTE: regions use 1-based indexing,
    as per the general and Ensembl standard

    :param        seq_region: the seq region name, usually the chromosome name
    :type         seq_region: str
    :param        sample_len: length of region sample to create
    :type         sample_len: int
    :param seq_region_length: the total length of the seq region (chromosome)
    :type  seq_region_length: int

    :return: region feature dictionary
    :rtype : dict {'seq_region': str, 'start': int, 'end': int, 'strand': int, 'len': int}
    """
    # get random end between min/max end point
    end = random.randint(sample_len, seq_region_length)
    # calculate start point,
    start = 1 + end - sample_len
    # Randomm strand
    strand = random.choice([1,-1])
    return {
        'seq_region': seq_region,
        'start': start,
        'end': end,
        'strand': strand,
        'len': end - start + 1
    }

def get_chromosomes(species, server=default_server):
    """get a dict of chromosome names and lengths for a species
    (blocking get request)

    :param species: ensembl species name
    :type  species: str
    :param  server: ensembl REST server URL, default https://rest.ensembl.org
    :type   server: str

    :return: dict of chromosome names, lengths as keys, values
    :rtype : dict['nae
    """
    ext = f"/info/assembly/{species}?"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
    res = r.json()
    chromosomes = {}
    for r in res['top_level_region']:
        if r['coord_system'] == 'chromosome':
            chromosomes[r['name']] = r['length']
    # return dict of chromosome names and lengths
    return chromosomes

async def region_seq(species, loc, session=None, server=default_server):
    """region_seq
    (async get request)

    Return the (DNA) sequence (str) for a species location

    :param species: ensembl species name
    :type  species: str
    :param     loc: ensembl location string
    :type      loc: str
    :param session: aiohttp client session, default None
    :type  session: aiohttp.ClientSession
    :param  server: ensembl REST server URL, default https://rest.ensembl.org
    :type   server: str

    :return: region DNA sequence
    :rtype : str (asyncio awaitable)
    """
    ext = f"/sequence/region/{species}/{loc}"
    async with session.get(ext, headers={"Content-Type":"text/plain"}) as r:
        if not r.ok:
            print(f"{r.status}\n{r.headers}\n", file=sys.stderr)
            r.raise_for_status()

        return await r.text()

async def region_labels(
        species,
        region,
        feature_types='exon',
        biotype='protein_coding',
        session=None,
        server=default_server
        ):
    """
    region_labels
    (async get request)

    This should return a numpy ndarray of zeros and ones
    with shape (region_length, len(feature_types)). That is
    a column of exon labels, and a column of CDS labels,
    by default.

    Transcript splice variants are handled by selecting only cannonical
    transcript sub-features.

    TODO (maybe): Add 'merge' option to instead combine overlapping features of
    the same type, if in the same gene, e.g. all exons from all transcripts of
    a particular gene.

    :param       species: ensembl species name
    :type        species: str
    :param        region: region dict, e.g from gabi.utils.get_random_region
    :type         region: dict
    :param feature_types: sub-feature tag(s) for labels, default exon
    :type  feature_types: str or [str]
    :param       biotype: limit the transcripts to biotype, default protein_coding
    :types       biotype: str
    :param       session: aiohttp client session, default None
    :type        session: aiohttp.ClientSession
    :param        server: ensembl REST server URL, default https://rest.ensembl.org
    :type         server: str

    :return: array of labels with shape (region length, num feature_types)
    :rtype : numpy.ndarray
    """
    # listify feature_types
    if isinstance(feature_types, str):
        feature_types = [feature_types]
    # make feature query string
    fstr = ';'.join([f"feature={f.lower()}" for f in feature_types])
    # get transcripts and target features in region
    loc = get_region_str(region)
    ext = f"/overlap/region/{species}/{loc}?feature=transcript;{fstr};biotype={biotype}"
    if not session: session = get_session(server)
    async with session.get(ext, headers={ "Content-Type" : "application/json"}) as r:
        if not r.ok:
            print(f"{r.status}\n{r.headers}\n", file=sys.stderr)
            r.raise_for_status()
        decoded = await r.json()

    # init label arrays of zeros for each feature_types
    labels = { k : np.zeros((region['len']), dtype=np.int8) for k in feature_types }

    if not decoded:
        return labels # no features, return labels array of zeros

    # Sort by length ascending
    data = sorted(decoded, key=lambda x: x['end'] - x['start'])
    # transform results into dict on ids
    # take longest where IDs clash since sorted by increasing length
    features = {f['id']:f for f in data}
    # filter
    label_features = {}
    for fid, f in features.items():
        # Do we want this feature type?
        if f['feature_type'].lower() not in map(str.lower, feature_types):
            continue
        # Skip only if feature is specifically non-cannonical
        if features[ f['Parent'] ].get('is_canonical', 1) == 0:
            continue
        # Check extents
        if f['end'] < region['start'] or f['start'] > region['end']:
            continue
        label_features[fid] = f
    # No desired feature types
    if not label_features:
        return labels

    # Transfer features to labels
    for fid, f in label_features.items():
        # Get ZERO-BASED indices for feature, relative to region
        # truncate at start of region
        start = max(0, f['start'] - region['start'])
        # truncate at end of region
        end   = min(region['len'] - 1, f['end'] - region['start'])
        # assign 1 to feature's range in feature_type label array
        labels[f['feature_type']][start:end] = 1
    return labels

async def sample_generator(
        species='homo_sapiens',
        feature_types='exon',
        biotype='protein_coding',
        n=None,
        sample_len=1000,
        server=default_server,
        seed=None
        ):
    """
    Generate sample regions for a species
    (async generator)

    USAGE
    -----

    import gabi.utils
    # create the default generator
    sgen = gabi.utils.sample_generator()

    # Get one sample
    r = await sgen.__anext__()
    print(r['seq_region']) # region name
    print(r['start'])
    print(r['end'])
    print(r['strand'])
    print(r['loc'])
    print(r['len'])
    print(sorted(set(r['seq']))) # check all nucleotide types in seq
    print(r['labels'].shape) # check shape of labels array

    # get samples until generator is empty
    async for r in sgen:
        print(r['loc'])

    :param       species: ensembl species name
    :type        species: str
    :param feature_types: sub-feature tag(s) for labels, default exon
    :type  feature_types: str or [str]
    :param       biotype: limit the transcripts to biotype, default protein_coding
    :types       biotype: str
    :param             n: max samples, default None yields infinite samples
    :type              n: int || None
    :param    sample_len: length of sample regions, default 1000
    :type     sample_len: int
    :param        server: ensembl REST service URL, default https://rest.ensembl.org
    :type         server: str
    :param          seed: set the random seed for sampling, default None (don't)
    :type           seed: int

    :return: generator of region dicts, with sequence and labels
    :rtype : async_generator
    """
    # set random seed?
    if seed:
        random.seed(seed)
        np.random.seed(seed)

    async with aiohttp.ClientSession(base_url=server) as session:
        # Get the chromosomes and their lengths
        chr_dict = get_chromosomes(species, server=server)
        chr_names, chr_lengths = zip(*chr_dict.items())
        count = 0
        # Yield forever if n is None
        while (n is None or count < n):
            # Weighted (by size) random choice of chromosome
            seq_id = random.choices(chr_names, chr_lengths, k=1)[0]
            #print(seq_id, sample_len, chr_dict)
            # Get random region on chosen chromosome
            region = get_random_region(seq_id, sample_len, chr_dict[seq_id])
            region['species'] = species # attach species name

            # Get and store a region location string
            region_loc = get_region_str(region)
            region['loc'] = region_loc
            
            # get sample region componenets asynchronously,
            # note result is ordered by input order
            region['seq'], region['labels'] = await asyncio.gather(
                # Get region sequence
                region_seq(species, region_loc, session=session),
                # get region labels
                region_labels(species, region, feature_types=feature_types, biotype=biotype, session=session),
            )

            # Increment generator count
            count += 1
            # yield the region dict
            yield region

def encode(seq: str, enc: dict) -> dict:
    """Encode (into numpy array uint8) a sequence str with and enc(oding) dict"""
    encoded = None
    try:
        encoded = np.array(
        [c for s in seq.upper() for c in enc[s]], # flat list
        dtype="uint8",
        )
    except KeyError as e:
        print(f"Cannot encode nucleotide {e} with this encoder! {enc}", file=sys.stderr)
        return None
    return {
        'seq': encoded,
        'bits': enc['bits'],
        'id': enc['id'],
        }

