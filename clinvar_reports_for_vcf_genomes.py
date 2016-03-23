import bz2
import csv
import datetime
import gzip
import os
import re

import vcf2clinvar
from vcf2clinvar import clinvar_update
import myvariant

FILESDIR = 'source_files'


def guess_getevidence_url(clinvar_name):
    """
    Guess where GET-Evidence would have a page, if it has one.

    This URL doesn't mean a page necessarily exists. It's our best guess
    where to find one. GET-Evidence data is stored by gene and amino acid
    change (with some nonstandard format, implemented below).
    """
    base_url = 'http://evidence.pgp-hms.org/'
    reFS = r'\((.*?)\):.*?\(p\.(.*?[0-9]+)[A-Z][a-z][a-z]fs\)'
    reAASub = r'\((.*?)\):.*?\(p\.([A-Za-z]+[0-9]+[A-Za-z]+)\)'
    variant = ''
    if re.search(reFS, clinvar_name):
        gene, pos = re.search(reFS, clinvar_name).groups()
        variant = '{}-{}Shift'.format(gene, pos)
    elif re.search(reAASub, clinvar_name):
        gene, aachange = re.search(reAASub, clinvar_name).groups()
        aachange = re.sub(r'Ter', 'Stop', aachange)
        variant = '{}-{}'.format(gene, aachange)
    if variant:
        return '{}{}'.format(base_url, variant)
    else:
        return ''


def match_genome(inputfile, outputfile, inputfilename):
    """
    Produce a CSV genome report at outputfile for a given VCF inputfile.
    """
    data = dict()

    # Set up ClinVar data.
    clinvar_filepath = clinvar_update.get_latest_vcf_file(FILESDIR, 'b37')
    if clinvar_filepath.endswith('.vcf'):
            input_clinvar_file = open(clinvar_filepath)
    elif clinvar_filepath.endswith('.vcf.gz'):
        input_clinvar_file = gzip.open(clinvar_filepath)
    elif clinvar_filepath.endswith('.vcf.bz2'):
        input_clinvar_file = bz2.BZ2File(clinvar_filepath)
    else:
        raise IOError("ClinVar filename expected to end with '.vcf'," +
                      " '.vcf.gz', or '.vcf.bz2'.")

    # Run vcf2clinvar on genome data.
    clinvar_matches = vcf2clinvar.match_to_clinvar(
        inputfile, input_clinvar_file)
    # Set up to get myvariant.info data (mainly for ExAC data.)
    mv = myvariant.MyVariantInfo()

    # iterate through all ClinVar matches.
    for genome_vcf_line, allele, zygosity in clinvar_matches:
        # Discard low quality data.
        if genome_vcf_line.filters and 'PASS' not in genome_vcf_line.filters:
            continue
        # Check significance. Only keep this as a notable variant if one of the
        # submissions has reported "pathogenic" and "likely pathogenic" effect.
        sigs = [rec.sig for rec in allele.records]
        if not ('4' in sigs or '5' in sigs):
            continue
        # Store data in a dict according to HGVS position.
        poskey = myvariant.format_hgvs(
            genome_vcf_line.chrom,
            genome_vcf_line.start,
            genome_vcf_line.ref_allele,
            allele.sequence)
        data[poskey] = {'genome_vcf_line': genome_vcf_line,
                        'clinvar_allele': allele,
                        'zygosity': zygosity}

    # Add data from myvariant.info using the HGVS positions.
    variants = data.keys()
    mv_output = mv.getvariants(variants, fields=['clinvar', 'exac'])
    for i in range(len(variants)):
        if 'clinvar' in mv_output[i]:
            data[variants[i]]['mv_clinvar'] = mv_output[i]['clinvar']
        if 'exac' in mv_output[i]:
            data[variants[i]]['mv_exac'] = mv_output[i]['exac']

    # Write report as CSV.
    with open(outputfile, 'w') as f:
        csv_out = csv.writer(f)
        for var in variants:
            # Clinvar URL for variant.
            cv_url = 'http://www.ncbi.nlm.nih.gov/clinvar/{}/'.format(
                data[var]['clinvar_allele'].records[0].acc)
            disease_name = ''
            preferred_name = ''
            getev_url = ''
            # Disease name, preferred name, and GET-Evidence URL if we have
            # myvariant.info information with ClinVar data.
            if 'mv_clinvar' in data[var]:
                cv_url = 'http://www.ncbi.nlm.nih.gov/clinvar/variation/{}/'.format(
                    data[var]['mv_clinvar']['variant_id'])
                try:
                    disease_name = data[var]['mv_clinvar']['rcv']['conditions']['name']
                    preferred_name = data[var]['mv_clinvar']['rcv']['preferred_name']
                except TypeError:
                    disease_name = ', '.join(
                        set([rcv['conditions']['name'] for rcv in
                            data[var]['mv_clinvar']['rcv']]))
                    preferred_name = data[var]['mv_clinvar']['rcv'][0]['preferred_name']
                getev_url = guess_getevidence_url(preferred_name)
            exac_url = 'http://exac.broadinstitute.org/variant/{}-{}-{}-{}'.format(
                data[var]['genome_vcf_line'].chrom[3:],
                data[var]['genome_vcf_line'].start,
                data[var]['genome_vcf_line'].ref_allele,
                data[var]['clinvar_allele'].sequence)
            # Allele frequency using ExAC data, if myvariant.info had that.
            if 'mv_exac' in data[var]:
                total_freq = data[var]['mv_exac']['ac']['ac'] * 1.0 / data[var]['mv_exac']['an']['an']
                total_freq = str(total_freq)
                freq_source = 'ExAC'
            else:
                # If not, try to get it from our ClinVar data.
                try:
                    total_freq = str(data[var]['clinvar_allele'].frequency)
                    freq_source = 'ClinVar'
                except KeyError:
                    # If that fails, give up on frequency.
                    total_freq = ''
                    freq_source = 'Unknown'
            data_row = [
                inputfilename, var, preferred_name, disease_name, cv_url,
                exac_url, total_freq, freq_source, getev_url]
            csv_out.writerow(data_row)
    return


if __name__ == "__main__":
    """
    Produce a series of CSV genome reports for a directory of VCF genomes

    Run on the command line with 1st argument as VCF directory, 2nd argument
    the directory to produce reports. This script will skip reproducing
    a CSV report if the planned filename is already present in that output
    directory.
    """
    genomes_dir = 'vcf_directory'
    output_dir = 'report_directory'
    for filename in os.listdir(genomes_dir):
        if filename == '.placeholder':
            continue
        print "{}  Processing {}...".format(str(datetime.datetime.now()), filename)
        filepath = os.path.join(genomes_dir, filename)
        if filename.endswith('.gz'):
            inputfile = gzip.open(filepath)
            filename_stripped = filename[0:-3]
        elif filename.endswith('.bz2'):
            inputfile = bz2.BZ2File(filepath)
            filename_stripped = filename[0:-4]
        else:
            inputfile = open(filepath)
            filename_stripped = filename
        outputpath = os.path.join(output_dir, filename_stripped + '.csv')
        if os.path.exists(outputpath):
            print "Skipping {}: report already generated.".format(outputpath)
            continue
        match_genome(inputfile=inputfile, outputfile=outputpath,
                     inputfilename=filename)
