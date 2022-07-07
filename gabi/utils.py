################################################################################
################################# SAMPLE DATA ##################################
################################################################################
import sys
import requests
import aiohttp
import asyncio
#import json
import numpy as np
import random
#from collections import Counter

# Default Server URL
default_server='https://rest.ensembl.org'
# get an aiohttp session
def get_session(server=default_server):
    return aiohttp.ClientSession(base_url=server)

# Get an Ensembl region string from a region dict
get_region_str = lambda r: f"{r['seq_region']}:{r['start']}..{r['end']}:{r['strand']}"

# BLOCKING
def get_random_region(seq_region, sample_len, seq_region_length):
    """
    NOTE: use 1-based indexing for Ensembl region query
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

# BLOCKING
def get_chromosomes(species, server=default_server):
    """ species info """
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

# ANSYNC
async def region_seq(
        species,
        loc,
        session=None,
        server=default_server
        ):
    ext = f"/sequence/region/{species}/{loc}"
    async with session.get(ext, headers={ "Content-Type" : "text/plain"}) as resp:
        if not resp.ok:
            print(f"{resp.status}\n{resp.headers}\n", file=sys.stderr)
            resp.raise_for_status()

        return await resp.text()

# ANSYNC
async def region_labels(
        species,
        region,
        feature_types='exon',
        transcript_strategy='cannonical',
        biotype='protein_coding',
        session=None,
        server=default_server
        ):
    """
    This should return a numpy ndarray of zeros and ones
    with shape (region_length, len(feature_types)). That is
    a column of exon labels, and a column of CDS labels,
    by default.
    Transcript splice variants are handled by transcript_strategy.
    'cannonical' selects only cannonical transcript sub-features.
    'merge' combines overlapping features of the same type,
    if in the same gene, e.g. all exons from all transcripts
    of a particular gene.
    """
    # listify feature_types
    if isinstance(feature_types, str):
        feature_types = [feature_types]
    # make feature query string
    fstr = ';'.join([f"feature={f}" for f in feature_types])
    # get transcripts and target features in region
    loc = get_region_str(region)
    ext = f"/overlap/region/{species}/{loc}?feature=transcript;{fstr};biotype={biotype}"
    if not session: session = get_session(server)
    async with session.get(ext, headers={ "Content-Type" : "application/json"}) as resp:
        if not resp.ok:
            print(f"{resp.status}\n{resp.headers}\n", file=sys.stderr)
            resp.raise_for_status()
        # init labels of zeros
        decoded = await resp.json()
    labels = np.zeros((region['len'], len(feature_types)), dtype=np.int8)
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
        # Check for target types
        if f['feature_type'] not in feature_types:
            continue
        # Check parent if canonical mode
        if features[ f['Parent'] ]['is_canonical'] == 0:
            continue
        # Check extents
        if f['end'] < region['start'] or f['start'] > region['end']:
            continue
        label_features[fid] = f
    # No desired feature types
    if not label_features:
        return labels

    # We'll need to "randomly" access the labels array,
    # columns (feature types), consistently by index
    feature_cols = {k:v for v, k in enumerate(feature_types)}
    # Transfer features to labels
    for fid, f in label_features.items():
        # Get ZERO-BASED indices for feature, relative to region
        start = max(0, f['start'] - region['start']) # truncate at start of region
        end   = min(region['len'] - 1, f['end'] - region['start']) # truncate at end of region
        # get the column index to transfer to
        col = feature_cols[ f['feature_type'] ]
        # assign label 1 to feature's position range on it's label list (array column)
        labels[start:end, col] = 1
    return labels

async def sample_generator(species='homo_sapiens', n=10, sample_len=1000, server=default_server, seed=None):
    """ Generator """
    if seed: random.seed(seed)

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

            # Get and store a region location string
            region_loc = get_region_str(region)
            region['loc'] = region_loc
            
            # Get region sequence
            region['seq'] = await region_seq(species, region_loc, session=session)
            # get region labels
            region['labels'] = await region_labels(species, region, session=session)

            # Increment generator count
            count += 1
            # yield the region dict
            yield region
